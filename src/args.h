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

#ifndef _KYCG_ARGS_H
#define _KYCG_ARGS_H

/**
 * Reorder argv so every option precedes every operand.
 *
 * WHY THIS EXISTS
 *   POSIX says getopt() stops at the first operand. glibc ignores that and
 *   permutes argv so options may be interleaved; the BSD getopt() shipped on
 *   macOS obeys it and stops. The same command line therefore means two
 *   different things on Linux and on macOS:
 *
 *     kycg fetch -f mm10:CGI -d /tmp/store
 *
 *   On Linux -d is an option. On macOS it is a *target name*, and so is
 *   "/tmp/store" -- kycg fetched into the default store, then reported two
 *   unknown targets, having already done half the work. Nothing warned that
 *   the flag had been reinterpreted, which is the part that makes this worth
 *   fixing rather than documenting.
 *
 *   Calling this before getopt() gives every platform glibc's behavior, which
 *   is what users expect and what every example in the README assumes.
 *
 * SEMANTICS
 *   Options (with their arguments) keep their relative order and move to the
 *   front; operands keep their relative order and follow. Both "-d DIR" and
 *   "-dDIR" are understood, as are clusters like "-fr" and "-frd DIR".
 *   A lone "-" is an operand, per convention. At "--" permutation stops and
 *   everything from there on is left exactly where it is, so a file genuinely
 *   named "-f" is still reachable.
 *
 *   argv is permuted in place and no string is copied, so the caller's
 *   pointers stay valid. argc must be the real count and argv[argc] is not
 *   touched.
 */
void kycg_permute_args(int argc, char *argv[], const char *optstring);

#endif /* _KYCG_ARGS_H */
