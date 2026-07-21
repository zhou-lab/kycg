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
 * Tests for argv permutation.
 *
 * These exist because the bug they guard against is silent and
 * platform-specific: on macOS, `kycg fetch -f mm10:CGI -d /tmp/store` read -d
 * as a target name, fetched into the wrong store, and only then complained.
 * A test that runs everywhere is the only way that stays fixed, since on
 * glibc the unpermuted code passes by accident.
 */

#include "args.h"

#include <stdio.h>
#include <string.h>

static int failures = 0;

/* Compare the permuted argv against an expected space-joined ordering. */
static void check(const char *name, const char *optstring,
                  const char *in, const char *want) {
  char buf[512];
  char *argv[64];
  int argc = 0;

  snprintf(buf, sizeof(buf), "%s", in);
  for (char *t = strtok(buf, " "); t && argc < 63; t = strtok(NULL, " "))
    argv[argc++] = t;
  argv[argc] = NULL;

  kycg_permute_args(argc, argv, optstring);

  char got[512] = {0};
  for (int i = 0; i < argc; ++i) {
    if (i) strncat(got, " ", sizeof(got) - strlen(got) - 1);
    strncat(got, argv[i], sizeof(got) - strlen(got) - 1);
  }

  if (strcmp(got, want) != 0) {
    printf("  FAIL %s\n    in   %s\n    got  %s\n    want %s\n",
           name, in, got, want);
    ++failures;
  }
}

int main(void) {
  const char *fetch = "d:o:t:frh";
  const char *test  = "m:a:Gs:MFHo:h";

  /* The reported bug, verbatim. */
  check("option after operand", fetch,
        "fetch -f mm10:CGI -d /tmp/store",
        "fetch -f -d /tmp/store mm10:CGI");

  /* Already correct: must not be disturbed. */
  check("already ordered", fetch,
        "fetch -f -d /tmp/store mm10:CGI",
        "fetch -f -d /tmp/store mm10:CGI");

  /* Nothing to do. */
  check("operands only", fetch, "fetch hg38 mm10", "fetch hg38 mm10");
  check("options only", fetch, "fetch -f -r", "fetch -f -r");

  /* Attached argument: the next element is an operand, not the option's arg. */
  check("attached arg", fetch,
        "fetch -d/tmp/store mm10",
        "fetch -d/tmp/store mm10");

  /* Cluster ending in an argument-taking option consumes the next element. */
  check("cluster then arg", fetch,
        "fetch mm10 -fd /tmp/store",
        "fetch -fd /tmp/store mm10");

  /* Cluster of flags only: the following element stays an operand. */
  check("flag cluster", fetch,
        "fetch mm10 -fr",
        "fetch -fr mm10");

  /* Relative order is preserved within each class. */
  check("order preserved", test,
        "test q.cg -m a.cm -G -m b.cm",
        "test -m a.cm -G -m b.cm q.cg");

  /* "--" stops permutation; the tail is untouched so a file named "-f" is
   * still reachable. */
  check("double dash", fetch,
        "fetch mm10 -- -f",
        "fetch mm10 -- -f");
  check("double dash with earlier option", fetch,
        "fetch mm10 -r -- -f",
        "fetch -r mm10 -- -f");

  /* A lone "-" is an operand, not an option. */
  check("lone dash", test, "test - -G", "test -G -");

  /* An unknown option is left for getopt to reject, not swallowed. */
  check("unknown option", fetch, "fetch mm10 -z", "fetch -z mm10");

  /* Missing argument at end of argv: nothing to consume, no read past argc. */
  check("dangling arg option", fetch, "fetch mm10 -d", "fetch -d mm10");

  if (failures) {
    printf("  %d check(s) failed\n", failures);
    return 1;
  }
  printf("  all checks passed\n");
  return 0;
}
