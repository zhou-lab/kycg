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
 * kycg — functional analysis of DNA methylation at CpG resolution.
 *
 * Subcommand dispatcher. Every subcommand lives in its own translation unit
 * and exposes int main_<cmd>(int, char**); this file holds the only main().
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "kycg.h"
#include "registry.h"
#include "store.h"
#include "ui.h"
#include "yame_version.h"

/**
 * What this build is bound to.
 *
 * A build can verify only the tags whose digests were compiled into it, so
 * "which generation of the data does this kycg speak to" is a real question
 * with a fixed answer. It belongs in the help text rather than behind a flag
 * nobody knows to type.
 */
/*
 * The help text is styled, but only when someone is looking: every colour
 * below comes from ui.c, which returns an empty string off a TTY, under
 * NO_COLOR, or on a dumb terminal. So `kycg 2>&1 | less` stays readable and
 * the bytes are identical to what they were before any of this.
 */
#define H_TITLE  kycg_ui_bold()
#define H_KEY    kycg_ui_cyan()
#define H_NOTE   kycg_ui_dim()
#define H_OFF    kycg_ui_reset()

/**
 * The two facts worth stating right under the title: the backend this build is
 * coupled to, and where its store is. kycg links libyame.a from a pinned
 * submodule, so the YAME version is fixed at build time -- the store layout,
 * the fetch/verify engine and the `.cx` formats all come from it; every fetch
 * and every -m path resolves against the store directory.
 */
static void print_build_info(FILE *out) {
  const char *env = getenv("KYCG_DATA_DIR");

  fprintf(out, "    %sbuilt against%s  YAME %s\n", H_NOTE, H_OFF, YAME_VERSION);
  fprintf(out, "    %sstore%s          %s   %s%s%s\n",
          H_NOTE, H_OFF, kycg_store_root(NULL),
          H_NOTE, env && *env ? "(from $KYCG_DATA_DIR)"
                              : "($KYCG_DATA_DIR unset; -d overrides)", H_OFF);
}

static void cmd(FILE *out, const char *name, const char *what) {
  fprintf(out, "    %s%-8s%s %s\n", H_KEY, name, H_OFF, what);
}

static int usage(void) {
  FILE *o = stderr;
  fprintf(o, "\n");
  fprintf(o, "  %skycg%s %s%s%s  %s— functional analysis of DNA methylation "
             "at CpG resolution%s\n\n",
          H_TITLE, H_OFF, H_KEY, KYCG_VERSION, H_OFF, H_NOTE, H_OFF);

  print_build_info(o);
  fprintf(o, "\n");

  fprintf(o, "%sUsage%s\n", H_TITLE, H_OFF);
  fprintf(o, "    kycg %s<command>%s [options]\n\n", H_KEY, H_OFF);

  fprintf(o, "%sCommands%s\n", H_TITLE, H_OFF);
  cmd(o, "fetch", "browse, choose and download knowledgebases");
  cmd(o, "test",  "set enrichment against a knowledgebase");
  cmd(o, "annotate", "label a TSV of probe IDs by set membership");
  cmd(o, "info",  "describe the records in a .cg or .cm file");
  fprintf(o, "\n");

  return 1;
}

int main(int argc, char *argv[]) {
  int ret;

  if (argc < 2) return usage();

  if      (strcmp(argv[1], "test") == 0)  ret = main_test(argc - 1, argv + 1);
  else if (strcmp(argv[1], "info") == 0)  ret = kycg_main_info(argc - 1, argv + 1);
  else if (strcmp(argv[1], "fetch") == 0) ret = kycg_main_fetch(argc - 1, argv + 1);
  else if (strcmp(argv[1], "annotate") == 0) ret = main_annotate(argc - 1, argv + 1);
  else if (strcmp(argv[1], "-h") == 0 ||
           strcmp(argv[1], "--help") == 0) return usage();
  else if (strcmp(argv[1], "--version") == 0) {
    /* Title, plus what the build is coupled to and points at: the YAME version
     * (fixed by the pinned submodule) and the store directory. */
    printf("kycg %s\n", KYCG_VERSION);
    print_build_info(stdout);
    return 0;
  } else {
    fprintf(stderr, "[main] Unrecognized command '%s'.\n", argv[1]);
    return usage();
  }

  /* YAME's main() flushes and closes stdout explicitly before returning,
   * noting that this is "not enough for remote file systems" but is what
   * catches the common case of a full or failing output device. Preserve the
   * behavior here so kycg fails the same way under the same conditions. */
  fflush(stdout);
  fclose(stdout);

  return ret;
}
