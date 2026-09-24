// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
/**
 * This file is part of kycg.
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia.
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes subject to the terms of the 2-Clause BSD License.
 * For use or transfers of the software to commercial entities, please inquire
 * with Dr. Wanding Zhou at zhouw3@chop.edu.
 *
 * See the LICENSE file at the root of this repository for the full terms and
 * the warranty disclaimer.
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
#include "yame_ui.h"
#include "yame_version.h"
#include "assets.h"    /* yame_assets_have_curl */

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
 * below comes from YAME's ui, which returns an empty string off a TTY or on a
 * dumb terminal. So `kycg 2>&1 | less` stays readable and the bytes are
 * identical to what they were before any of this. NO_COLOR is not consulted --
 * see the note at the end of tests/t_browser.sh.
 */
#define H_TITLE  yame_ui_bold()
#define H_KEY    yame_ui_cyan()
#define H_NOTE   yame_ui_dim()
#define H_OFF    yame_ui_reset()

/**
 * The two facts worth stating right under the title: the backend this build is
 * coupled to, and where its store is. kycg links libyame.a from a pinned
 * submodule, so the YAME version is fixed at build time -- the store layout,
 * the fetch/verify engine and the `.cx` formats all come from it; every fetch
 * and every -m path resolves against the store directory.
 */
static void print_build_info(FILE *out) {
  const char *env = getenv("YAME_DATA_HOME");

  fprintf(out, "    %sbuilt against%s  YAME %s\n", H_NOTE, H_OFF, YAME_VERSION);
  fprintf(out, "    %sstore%s          %s   %s%s%s\n",
          H_NOTE, H_OFF, kycg_store_root(NULL),
          H_NOTE, env && *env ? "(from $YAME_DATA_HOME)"
                              : "($YAME_DATA_HOME unset; -d overrides)", H_OFF);
}

static void cmd(FILE *out, const char *name, const char *what) {
  fprintf(out, "    %s%-8s%s %s\n", H_KEY, name, H_OFF, what);
}

/* What data this build can verify, and whether it can fetch it. A build trusts
 * only the tags whose digests are compiled in, so --version states them; and
 * fetching needs libcurl, which lives in libyame and is a runtime fact. Both
 * are things the README says --version reports, so it does. */
static void print_registry_info(FILE *out) {
  size_t n_arr = 0, n_seq = 0;
  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r) ++n_arr;
  for (const kycg_seq_reg_t   *r = KYCG_SEQ_REGISTRY;   r->genome;   ++r) ++n_seq;
  /* Tags are per file now, so there is no single pinned tag to name. State
   * the sources this build carries rows from, and what it can verify. */
  fprintf(out, "    %sregistry%s       %s  (%zu arrays, %zu genomes)\n",
          H_NOTE, H_OFF, KYCG_REGISTRY_TAGS, n_arr, n_seq);
  fprintf(out, "    %snetwork%s        libcurl %s\n",
          H_NOTE, H_OFF,
          yame_assets_have_curl() ? "available"
                                  : "absent (built without it; fetch disabled)");
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
  fprintf(o, "\n");

  return 1;
}

int main(int argc, char *argv[]) {
  int ret;

  if (argc < 2) return usage();

  if      (strcmp(argv[1], "test") == 0)  ret = main_test(argc - 1, argv + 1);
  else if (strcmp(argv[1], "fetch") == 0) ret = kycg_main_fetch(argc - 1, argv + 1);
  else if (strcmp(argv[1], "annotate") == 0) ret = main_annotate(argc - 1, argv + 1);
  else if (strcmp(argv[1], "-h") == 0 ||
           strcmp(argv[1], "--help") == 0) { usage(); return 0; }
  else if (strcmp(argv[1], "--version") == 0 ||
           strcmp(argv[1], "-v") == 0) {
    /* Title, plus what the build is coupled to: the YAME version (fixed by the
     * pinned submodule), the store, the registry tag it can verify, and
     * whether it can fetch at all. */
    printf("kycg %s\n", KYCG_VERSION);
    print_build_info(stdout);
    print_registry_info(stdout);
    return 0;
  } else {
    fprintf(stderr, "kycg: unrecognized command '%s'.\n", argv[1]);
    usage();
    return 1;
  }

  /* YAME's main() flushes and closes stdout explicitly before returning,
   * noting that this is "not enough for remote file systems" but is what
   * catches the common case of a full or failing output device. Preserve the
   * behavior here so kycg fails the same way under the same conditions. */
  fflush(stdout);
  fclose(stdout);

  return ret;
}
