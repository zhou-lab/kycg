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
 * Tests for the store's member-name guard.
 *
 * This is a security control, so it gets a test that states the attack rather
 * than the implementation. Names in a downloaded SHA256SUMS were previously
 * joined straight onto the store directory, so an entry reading
 * "../../../../.ssh/authorized_keys" wrote outside the store -- and because
 * whoever writes the manifest also serves the bytes, the digest check passed.
 *
 * The accepted cases matter as much as the rejected ones: this runs on every
 * line of every manifest, and a guard that rejected a legitimate set name
 * would break fetching entirely.
 */

#include "store.h"

#include <stdio.h>

static int failures = 0;

static void want(const char *name, int expect, const char *why) {
  int got = kycg_store_safe_name(name);
  if (got != expect) {
    printf("  FAIL %-38s got %d, want %d  (%s)\n",
           name ? name : "(null)", got, expect, why);
    ++failures;
  }
}

int main(void) {
  /* Real entries, taken from the published manifests. All must pass. */
  want("CGI.20220904.cm",                 1, "ordinary set");
  want("ChromHMM.20220414.cm",            1, "ordinary set");
  want("CGI.20220904.cm.idx",             1, "seek index");
  want("cpg_nocontig.cr",                 1, "genome reference");
  want("MSA.ordering.tsv.gz",             1, "array probe ordering");
  want("SHA256SUMS",                      1, "the manifest itself");
  want("TFBS.20220921.cm",                1, "set with digits");
  want("XCILinkedWGBSSorted.20221121.cm", 1, "long mixed-case name");

  /* Traversal: the attack this exists to stop. */
  want("../evil",                         0, "parent traversal");
  want("../../../../.ssh/authorized_keys",0, "escape to home");
  want("..",                              0, "bare parent");
  want(".",                               0, "bare self");
  want("a/../../b",                       0, "traversal mid-name");
  want("sub/dir.cm",                      0, "subdirectory");
  want("/etc/passwd",                     0, "absolute path");
  want("\\windows\\path.cm",              0, "backslash separator");

  /* Hidden files: never legitimate here, and a dotfile in the store would be
   * invisible to the enumeration in kycg_store_find_cm. */
  want(".bashrc",                         0, "hidden file");
  want(".git",                            0, "hidden directory");

  /* Control characters, which would corrupt terminal output and could hide
   * the real name from anyone reading the fetch plan. */
  want("evil\nname.cm",                   0, "embedded newline");
  want("evil\rname.cm",                   0, "embedded carriage return");
  want("evil\033[2Jname.cm",              0, "embedded escape sequence");

  /* Degenerate input. */
  want("",                                0, "empty name");
  want(NULL,                              0, "null pointer");

  if (failures) {
    printf("  %d check(s) failed\n", failures);
    return 1;
  }
  printf("  all checks passed\n");
  return 0;
}
