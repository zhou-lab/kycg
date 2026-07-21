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
 * Per-row emphasis.
 *
 * Passed alongside the rows rather than embedded in them: the widgets measure
 * and truncate row text to fit the terminal, and escape sequences inside that
 * text would be counted as visible cells and wreck the alignment.
 */
typedef enum {
  KYCG_ROW_PLAIN = 0,
  KYCG_ROW_HAVE,      /* present locally -- green */
  KYCG_ROW_MISSING,   /* not present -- dimmed */
} kycg_row_style_t;

/**
 * Scrollable in-place viewer for tabular output.
 *
 * `header` is a column header held fixed above the rows; it and `rows` are
 * tab-separated and rendered into aligned columns. `styles` may be NULL, or
 * one kycg_row_style_t per row. Returns 0 when the viewer ran, -1 when the
 * terminal cannot support it and the caller should print plainly instead.
 * Never call this when stdout is redirected: piped output must stay
 * machine-readable.
 */
int kycg_ui_browse(const char *title, const char *header,
                   const char **rows, const unsigned char *styles, size_t n);

/** Child rows of one expanded node, owned by the caller of the expand fn. */
typedef struct {
  char **rows;            /* preformatted lines; the tree indents, not aligns */
  char **keys;            /* optional; opaque, handed back on accept */
  unsigned char *styles;  /* optional, one per row */
  size_t n;
} kycg_ui_kids_t;

/**
 * Called once per checked row when the user accepts, before the tree returns.
 * `root` is the parent's row text, `key` the child's key from kycg_ui_kids_t.
 */
typedef void (*kycg_ui_accept_fn)(void *ctx, const char *root,
                                  const char *key);

/**
 * Fill `out` with the children of the node whose row is `row`. Called at most
 * once per node, the first time it is opened. Leaving out->n at 0 marks the
 * node as having nothing to show.
 */
typedef void (*kycg_ui_expand_fn)(void *ctx, const char *row,
                                  kycg_ui_kids_t *out);

/**
 * Two-level tree viewer: a table whose rows unfold in place.
 *
 * Right arrow (or l, or enter) opens the row under the cursor and splices its
 * children in beneath it; left (or h) closes it. Children are requested lazily
 * and kept, so opening a row twice costs one call.
 *
 * Children are rendered as given, only indented — the tree cannot align them,
 * because it sees them one parent at a time and column widths that shifted
 * with each expansion would be worse than none. Pre-format them to a common
 * width if they should line up.
 *
 * Same contract as kycg_ui_browse: returns -1 when the terminal cannot support
 * it, and must not be used when stdout is redirected.
 */
/**
 * When `accept` is non-NULL the tree is also a picker: rows carrying a key
 * and not already marked KYCG_ROW_HAVE get a checkbox, space toggles one,
 * space on a parent toggles all of its children, and `f` accepts.
 *
 * Rows styled KYCG_ROW_HAVE are shown as already present and cannot be
 * checked — there is nothing to ask for.
 *
 * Returns 1 if the user accepted a non-empty selection (every checked row
 * having been passed to `accept`), 0 if they quit, -1 if the terminal cannot
 * support the widget.
 */
int kycg_ui_tree(const char *title, const char *header,
                 const char **roots, const unsigned char *root_styles,
                 size_t n_roots, kycg_ui_expand_fn expand,
                 kycg_ui_accept_fn accept, void *ctx);

#endif /* _KYCG_UI_H */
