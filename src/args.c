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
 * argv permutation, so options may follow operands on every platform.
 *
 * See args.h for why this is needed at all. This is deliberately a separate
 * module with no dependency on anything else in kycg: it is a pure function
 * over argv, and keeping it that way is what makes it testable without a
 * process.
 */

#include "args.h"

#include <stdlib.h>
#include <string.h>

/**
 * Does `c` take an argument, per `optstring`?
 *
 * The ':' that marks the previous option must not itself match, or "-:" would
 * be read as an option taking an argument; and a leading ':' or '+' in the
 * optstring is a mode flag rather than an option letter. kycg uses plain
 * optstrings today, but this is the kind of thing that silently mis-parses
 * one flag years later if it is assumed rather than checked.
 */
static int takes_arg(const char *optstring, char c) {
  if (c == ':' || c == '+' || c == '-') return 0;
  const char *p = strchr(optstring, c);
  return p && p[1] == ':';
}

void kycg_permute_args(int argc, char *argv[], const char *optstring) {
  if (argc < 3 || !optstring) return;   /* nothing could be out of order */

  char **opts = malloc((size_t)argc * sizeof(char *));
  char **ops  = malloc((size_t)argc * sizeof(char *));
  /* Out of memory here is not worth failing the command over: leaving argv
   * untouched simply restores the platform's own getopt behavior. */
  if (!opts || !ops) { free(opts); free(ops); return; }

  int n_opts = 0, n_ops = 0, i = 1;

  for (; i < argc; ++i) {
    const char *a = argv[i];

    /* "--" ends option parsing for getopt too; leave the tail alone. */
    if (a[0] == '-' && a[1] == '-' && a[2] == '\0') break;

    /* A lone "-" is an operand by convention, not an empty option. */
    if (a[0] != '-' || a[1] == '\0') { ops[n_ops++] = argv[i]; continue; }

    opts[n_opts++] = argv[i];

    /* Find the first argument-taking letter in the cluster. If anything
     * follows it in the same string, that text *is* the argument ("-dDIR");
     * otherwise the argument is the next element ("-d DIR"). Letters after an
     * argument-taking one are part of the argument, not options, so the scan
     * stops at the first match either way. */
    size_t len = strlen(a);
    int attached = 0, wants = 0;
    for (size_t k = 1; k < len; ++k) {
      if (takes_arg(optstring, a[k])) {
        wants = 1;
        attached = (k + 1 < len);
        break;
      }
    }
    if (wants && !attached && i + 1 < argc) opts[n_opts++] = argv[++i];
  }

  /* Everything from "--" onward keeps its position, which is why the
   * write-back below reaches exactly that far and no further. */
  int w = 1;
  for (int k = 0; k < n_opts; ++k) argv[w++] = opts[k];
  for (int k = 0; k < n_ops;  ++k) argv[w++] = ops[k];

  free(opts);
  free(ops);
}
