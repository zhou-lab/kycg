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

#ifndef _KYCG_H
#define _KYCG_H

#include "ui.h"

#define KYCG_VERSION "0.1.0"

/*
 * Shared styling for help text. Each expands to an empty string off a TTY (see
 * ui.c), so a redirected --help is byte-identical to plain text.
 */
#define KYCG_H_TITLE kycg_ui_bold()
#define KYCG_H_KEY   kycg_ui_cyan()
#define KYCG_H_NOTE  kycg_ui_dim()
#define KYCG_H_WARN  kycg_ui_yellow()
#define KYCG_H_OFF   kycg_ui_reset()

/**
 * Resolve a knowledgebase spec to file paths in the store.
 *
 * Accepts what `kycg fetch` accepts -- "mm10", "mm10:CGI", "MSA:CGI,ChromHMM"
 * -- so a set is named the same way whether it is being downloaded or tested
 * against. A spec that is an existing file is returned as-is, which keeps
 * plain paths working and leaves no ambiguity to resolve.
 *
 * Returns the number of paths and fills *out with a malloc'd vector of
 * malloc'd strings (free with kycg_free_specs). Zero means nothing matched;
 * the caller reports that, since only it knows what to say.
 */
size_t kycg_resolve_spec(const char *spec, const char *store, char ***out);
void kycg_free_specs(char **v, size_t n);

/* Subcommand entry points, dispatched from main() on argv[1]. */
int main_test(int argc, char *argv[]);
int main_info(int argc, char *argv[]);
int main_fetch(int argc, char *argv[]);
int main_list(int argc, char *argv[]);

#endif /* _KYCG_H */
