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

#include "kycg.h"

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: kycg (functional analysis of DNA methylation)\n");
  fprintf(stderr, "Version: %s\n", KYCG_VERSION);
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   kycg <command> [options]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Commands:\n");
  fprintf(stderr, "    test      set enrichment against a knowledgebase\n");
  fprintf(stderr, "    info      describe the records in a .cg or .cm file\n");
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char *argv[]) {
  int ret;

  if (argc < 2) return usage();

  if      (strcmp(argv[1], "test") == 0) ret = main_test(argc - 1, argv + 1);
  else if (strcmp(argv[1], "info") == 0) ret = main_info(argc - 1, argv + 1);
  else if (strcmp(argv[1], "-h") == 0 ||
           strcmp(argv[1], "--help") == 0) return usage();
  else if (strcmp(argv[1], "--version") == 0) {
    printf("%s\n", KYCG_VERSION);
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
