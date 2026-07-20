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

#ifndef _KYCG_UI_H
#define _KYCG_UI_H

#include <stdint.h>
#include <stddef.h>

/**
 * Terminal presentation and interaction.
 *
 * Everything here degrades to plain text off a TTY. See ui.c for why that is
 * a correctness property and not a nicety.
 */

/** Is the terminal capable of animation / color? (stderr is a TTY, TERM sane,
 *  NO_COLOR unset.) */
int kycg_ui_fancy(void);

/** Can we ask the user a question? (stdin AND stderr are TTYs.) */
int kycg_ui_interactive(void);

/** Does the locale advertise UTF-8? Selects glyphs vs ASCII fallbacks. */
int kycg_ui_unicode(void);

/* SGR sequences, or "" when color is off. */
const char *kycg_ui_dim(void);
const char *kycg_ui_bold(void);
const char *kycg_ui_green(void);
const char *kycg_ui_red(void);
const char *kycg_ui_yellow(void);
const char *kycg_ui_cyan(void);
const char *kycg_ui_reset(void);

/** Glyphs: check, cross, bullet. UTF-8 or ASCII depending on locale. */
const char *kycg_ui_check(void);
const char *kycg_ui_cross(void);
const char *kycg_ui_bullet(void);

/** Human-readable byte count into `buf`; returns buf. */
const char *kycg_ui_human(uint64_t bytes, char *buf, size_t n);

/* ------------------------------------------------------------- progress */

/**
 * A single file's transfer, rendered as one self-updating line while it runs
 * and one settled line when it finishes.
 */
typedef struct {
  char     label[128];
  uint64_t total;        /* 0 when the server does not say */
  uint64_t now;
  int      frame;        /* spinner position */
  int      active;       /* a line is currently on screen */
  double   started;      /* monotonic seconds, for the rate readout */
  double   last_draw;
} kycg_prog_t;

void kycg_prog_begin(kycg_prog_t *p, const char *label, uint64_t total);
/** Safe to call at any rate; redraw is throttled internally. */
void kycg_prog_update(kycg_prog_t *p, uint64_t now, uint64_t total);
/** Settle the line: check mark on success, cross on failure. */
void kycg_prog_done(kycg_prog_t *p, const char *detail, int ok);

/** A one-off settled line for work that did not transfer (e.g. cached). */
void kycg_ui_line(const char *glyph_color, const char *glyph,
                  const char *label, const char *detail);

/* -------------------------------------------------------------- prompts */

/**
 * Yes/no. Returns 1 for yes, 0 for no.
 *
 * Callers MUST check kycg_ui_interactive() first and choose a non-interactive
 * default themselves; this function assumes it may block on stdin.
 */
int kycg_ui_confirm(const char *question, int default_yes);

/**
 * Free-text answer with a default. Returns a malloc'd string (the default when
 * the user just hits return), or NULL on EOF.
 */
char *kycg_ui_ask(const char *question, const char *def);

/**
 * Single-choice menu. Returns the chosen index, or -1 if cancelled.
 * `notes` may be NULL; when given, notes[i] is shown dimmed after items[i].
 *
 * On a capable terminal this is an in-place cursor list navigated with the
 * arrow keys; otherwise it degrades to a numbered prompt.
 */
long kycg_ui_choose(const char *title, const char **items, const char **notes,
                    size_t n);

/**
 * Multi-select. Returns a malloc'd array of n flags, or NULL if cancelled.
 *
 * In-place on a capable terminal: arrows or j/k to move, space to toggle,
 * a/n for all/none, / to filter, enter to accept. Otherwise a numbered prompt
 * accepting "all", "none", and comma/range lists such as "1-5,8,12".
 */
int *kycg_ui_multiselect(const char *title, const char **items,
                         const char **notes, size_t n, int default_all);

/**
 * Scrollable in-place viewer for tabular output.
 *
 * `header` is a column header held fixed above the rows; it and `rows` are
 * tab-separated and rendered into aligned columns. Returns 0 when the viewer
 * ran, -1 when the terminal cannot support it and the caller should print
 * plainly instead. Never call this when stdout is redirected: piped output
 * must stay machine-readable.
 */
int kycg_ui_browse(const char *title, const char *header,
                   const char **rows, size_t n);

#endif /* _KYCG_UI_H */
