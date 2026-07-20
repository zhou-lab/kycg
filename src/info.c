// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of kycg.
 *
 * Copyright (C) 2026-present Wanding Zhou
 *
 * kycg is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * kycg is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with kycg.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * `kycg info` — describe the records in a .cg or .cm file.
 *
 * GOAL
 *   Answer the two questions a user has before running `kycg test`: how many
 *   rows does this file index, and how many records does it carry. A row
 *   count mismatch between a query and a knowledgebase is the one failure mode
 *   kycg can detect, and this is how a user checks for it in advance rather
 *   than discovering it mid-run.
 *
 * Output is TSV so it composes with the rest of the toolchain.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>

#include "kycg.h"

/* YAME (submodule) */
#include "cfile.h"
#include "cdata.h"
#include "summary.h"
#include "wzmisc.h"

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg info [options] <file.cg|file.cm> [...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Report the format, row count, and set size of each record.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -H        suppress the header line\n");
  fprintf(stderr, "    -h        this help\n");
  fprintf(stderr, "\n");
  return 1;
}

/* Human-readable name for a YAME format code. */
static const char *fmt_name(char fmt) {
  switch (fmt) {
  case '0': return "set";           /* 1 bit per CpG                       */
  case '1': return "set_rle";       /* run-length encoded set              */
  case '2': return "categorical";   /* key table + states                  */
  case '3': return "mu";            /* M/U counts, mu==0 means missing     */
  case '4': return "beta";          /* float, NA-capable                   */
  case '5': return "value";
  case '6': return "set_universe";  /* 2 bits: hit / background / NA       */
  case '7': return "coordinates";   /* the .cr reference row list          */
  default:  return "unknown";
  }
}

int main_info(int argc, char *argv[]) {
  int no_header = 0;
  int c;
  while ((c = getopt(argc, argv, "Hh")) >= 0) {
    switch (c) {
    case 'H': no_header = 1; break;
    case 'h': return usage();
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind >= argc) {
    usage();
    wzfatal("Please supply an input file.\n");
  }

  if (!no_header) fputs("file\trecord\tname\tformat\tn_rows\tn_set\n", stdout);

  for (int j = optind; j < argc; ++j) {
    char *fname = argv[j];
    cfile_t cf = open_cfile(fname);
    snames_t snames = loadSampleNamesFromIndex(fname);
    const char *disp = get_basename(fname);

    for (uint64_t k = 0;; ++k) {
      cdata_t cd = read_cdata1(&cf);
      if (cd.n == 0) break;

      char fmt = cd.fmt;
      prepare_mask(&cd);   /* normalize so cd.n is a row count, not bytes */

      /* n_set is only meaningful for a bitset; other formats report NA. */
      char n_set[32];
      if (cd.fmt == '0') snprintf(n_set, sizeof(n_set), "%zu", bit_count(cd));
      else snprintf(n_set, sizeof(n_set), "NA");

      const char *name = (snames.n && k < (unsigned)snames.n)
        ? snames.s[k] : NULL;

      if (name) fprintf(stdout, "%s\t%" PRIu64 "\t%s\t", disp, k + 1, name);
      else      fprintf(stdout, "%s\t%" PRIu64 "\tNA\t", disp, k + 1);

      fprintf(stdout, "%c:%s\t%" PRIu64 "\t%s\n",
              fmt, fmt_name(fmt), cd.n, n_set);

      free_cdata(&cd);
    }

    bgzf_close(cf.fh);
    cleanSampleNames2(snames);
  }

  return 0;
}
