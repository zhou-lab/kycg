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
 * Terminal presentation: colors, spinners, progress bars, and prompts.
 *
 * WHY THIS IS A SEPARATE MODULE
 *   Every routine here has to answer the same question first -- is anyone
 *   actually watching? -- and getting that wrong is not a cosmetic bug.
 *   kycg is meant to run inside Nextflow pipelines and Docker builds, where
 *   stdin is closed and stderr is a log file. A prompt in that setting hangs
 *   the job forever with no indication of why, and an animated progress bar
 *   fills the log with thousands of carriage returns. Centralizing the
 *   capability checks means each caller cannot forget them independently.
 *
 * THE POLICY
 *   Interaction requires stdin AND stderr to both be TTYs. Animation requires
 *   stderr to be a TTY. Neither is ever assumed. Callers that want to prompt
 *   must check kycg_ui_interactive() and supply their own non-interactive
 *   behavior; nothing here silently blocks on a closed stdin.
 *
 *   This is the refinement of what DESIGN.md originally stated as "never
 *   prompt". The rule existed because a prompt would hang automation; gating
 *   prompts on an interactive terminal preserves that guarantee exactly, while
 *   letting a human at a keyboard get a usable interface.
 *
 * DEGRADATION
 *   Off a TTY: no escape sequences, no spinner, one plain line per event.
 *   Without a UTF-8 locale: ASCII glyphs ([ok], [xx], -) and an ASCII spinner.
 *   With NO_COLOR set or TERM=dumb: no color, everything else unchanged.
 *   The information content is identical in every mode -- only the ink differs.
 */

#include "ui.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>

/* ------------------------------------------------------- capability probes */

static int cached_fancy = -1, cached_inter = -1, cached_uni = -1;

int kycg_ui_fancy(void) {
  if (cached_fancy >= 0) return cached_fancy;

  const char *term = getenv("TERM");
  int ok = isatty(STDERR_FILENO)
        && !getenv("NO_COLOR")
        && !(term && strcmp(term, "dumb") == 0);
  cached_fancy = ok;
  return ok;
}

int kycg_ui_interactive(void) {
  if (cached_inter >= 0) return cached_inter;
  /* Both directions matter: a question we cannot display is as useless as one
   * whose answer we cannot read. */
  cached_inter = isatty(STDIN_FILENO) && isatty(STDERR_FILENO);
  return cached_inter;
}

int kycg_ui_unicode(void) {
  if (cached_uni >= 0) return cached_uni;
  const char *v = getenv("LC_ALL");
  if (!v || !*v) v = getenv("LC_CTYPE");
  if (!v || !*v) v = getenv("LANG");
  cached_uni = (v && (strcasestr(v, "utf-8") || strcasestr(v, "utf8"))) ? 1 : 0;
  return cached_uni;
}

/* --------------------------------------------------------------- palette */

#define COLOR(fn, seq) \
  const char *fn(void) { return kycg_ui_fancy() ? seq : ""; }

COLOR(kycg_ui_dim,    "\033[2m")
COLOR(kycg_ui_bold,   "\033[1m")
COLOR(kycg_ui_green,  "\033[32m")
COLOR(kycg_ui_red,    "\033[31m")
COLOR(kycg_ui_yellow, "\033[33m")
COLOR(kycg_ui_cyan,   "\033[36m")
COLOR(kycg_ui_reset,  "\033[0m")

const char *kycg_ui_check(void)  { return kycg_ui_unicode() ? "✓" : "ok"; }
const char *kycg_ui_cross(void)  { return kycg_ui_unicode() ? "✗" : "XX"; }
const char *kycg_ui_bullet(void) { return kycg_ui_unicode() ? "•" : "-"; }

/* Braille dots read as a smooth rotation; the ASCII fallback is the classic
 * four-frame spin. */
static const char *SPIN_U[] = {"⠋","⠙","⠹","⠸","⠼",
                               "⠴","⠦","⠧","⠇","⠏"};
static const char *SPIN_A[] = {"|","/","-","\\"};

static const char *spin_frame(int i) {
  if (kycg_ui_unicode()) return SPIN_U[((unsigned)i) % 10];
  return SPIN_A[((unsigned)i) % 4];
}

/* ---------------------------------------------------------------- helpers */

const char *kycg_ui_human(uint64_t bytes, char *buf, size_t n) {
  const char *unit[] = {"B", "KB", "MB", "GB", "TB"};
  double v = (double)bytes;
  int u = 0;
  while (v >= 1024.0 && u < 4) { v /= 1024.0; ++u; }
  snprintf(buf, n, "%.*f %s", (u == 0 || v >= 100) ? 0 : 1, v, unit[u]);
  return buf;
}

static double now_sec(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec / 1e6;
}

/* Terminal width, clamped to something sane. */
static int term_cols(void) {
  const char *c = getenv("COLUMNS");
  int w = c ? atoi(c) : 0;
  if (w <= 0) w = 80;
  if (w < 40) w = 40;
  if (w > 200) w = 200;
  return w;
}

/* Erase the current line, leaving the cursor at column 0. */
static void clear_line(void) {
  fputs("\r\033[2K", stderr);
}

/* ------------------------------------------------------------- progress */

#define LABEL_W 34
#define BAR_W   22
#define REDRAW_INTERVAL 0.08   /* seconds; ~12 fps is smooth and cheap */

void kycg_prog_begin(kycg_prog_t *p, const char *label, uint64_t total) {
  memset(p, 0, sizeof(*p));
  snprintf(p->label, sizeof(p->label), "%s", label);
  p->total = total;
  p->started = now_sec();
  p->last_draw = 0.0;
  p->active = kycg_ui_fancy();
}

void kycg_prog_update(kycg_prog_t *p, uint64_t now, uint64_t total) {
  if (!p->active) return;

  p->now = now;
  if (total) p->total = total;

  double t = now_sec();
  if (t - p->last_draw < REDRAW_INTERVAL) return;
  p->last_draw = t;
  ++p->frame;

  char hb_now[24], hb_tot[24], hb_rate[24];
  kycg_ui_human(p->now, hb_now, sizeof(hb_now));

  double elapsed = t - p->started;
  double rate = elapsed > 0.05 ? (double)p->now / elapsed : 0.0;
  kycg_ui_human((uint64_t)rate, hb_rate, sizeof(hb_rate));

  int cols = term_cols();
  clear_line();

  if (p->total) {
    double frac = (double)p->now / (double)p->total;
    if (frac > 1.0) frac = 1.0;
    int fill = (int)(frac * BAR_W + 0.5);

    kycg_ui_human(p->total, hb_tot, sizeof(hb_tot));

    /* The bar is drawn with block glyphs where the locale allows and '#'/'.'
     * where it does not; both are the same width in cells. */
    const char *full = kycg_ui_unicode() ? "█" : "#";
    const char *rest = kycg_ui_unicode() ? "░" : ".";

    fprintf(stderr, "%s%s%s %-*.*s ", kycg_ui_cyan(), spin_frame(p->frame),
            kycg_ui_reset(), LABEL_W, LABEL_W, p->label);

    fputs(kycg_ui_cyan(), stderr);
    for (int i = 0; i < BAR_W; ++i) fputs(i < fill ? full : rest, stderr);
    fputs(kycg_ui_reset(), stderr);

    fprintf(stderr, " %3d%%  %s%s/%s", (int)(frac * 100.0 + 0.5),
            kycg_ui_dim(), hb_now, hb_tot);
    if (cols > 96 && rate > 0.0) fprintf(stderr, "  %s/s", hb_rate);
    fputs(kycg_ui_reset(), stderr);
  } else {
    /* No Content-Length: show what has arrived, without a fake percentage. */
    fprintf(stderr, "%s%s%s %-*.*s %s%s%s",
            kycg_ui_cyan(), spin_frame(p->frame), kycg_ui_reset(),
            LABEL_W, LABEL_W, p->label,
            kycg_ui_dim(), hb_now, kycg_ui_reset());
  }

  fflush(stderr);
}

void kycg_prog_done(kycg_prog_t *p, const char *detail, int ok) {
  if (p->active) clear_line();
  kycg_ui_line(ok ? kycg_ui_green() : kycg_ui_red(),
               ok ? kycg_ui_check() : kycg_ui_cross(),
               p->label, detail);
  p->active = 0;
}

void kycg_ui_line(const char *glyph_color, const char *glyph,
                  const char *label, const char *detail) {
  fprintf(stderr, "%s%s%s %-*.*s %s%s%s\n",
          glyph_color ? glyph_color : "", glyph ? glyph : " ", kycg_ui_reset(),
          LABEL_W, LABEL_W, label ? label : "",
          kycg_ui_dim(), detail ? detail : "", kycg_ui_reset());
  fflush(stderr);
}

/* --------------------------------------------------------------- prompts */

/** Read a line from stdin, trimmed. Returns NULL on EOF. */
static char *read_line(void) {
  char buf[4096];
  if (!fgets(buf, sizeof(buf), stdin)) return NULL;

  size_t n = strlen(buf);
  while (n && (buf[n-1] == '\n' || buf[n-1] == '\r')) buf[--n] = '\0';

  char *p = buf;
  while (*p && isspace((unsigned char)*p)) ++p;
  size_t len = strlen(p);
  while (len && isspace((unsigned char)p[len-1])) p[--len] = '\0';

  return strdup(p);
}

int kycg_ui_confirm(const char *question, int default_yes) {
  for (;;) {
    fprintf(stderr, "%s%s%s %s [%s] ",
            kycg_ui_bold(), question, kycg_ui_reset(),
            kycg_ui_dim(), default_yes ? "Y/n" : "y/N");
    fputs(kycg_ui_reset(), stderr);
    fflush(stderr);

    char *ans = read_line();
    if (!ans) { fputc('\n', stderr); return 0; }   /* EOF declines */

    if (!*ans) { free(ans); return default_yes; }
    int y = (strcasecmp(ans, "y") == 0 || strcasecmp(ans, "yes") == 0);
    int n = (strcasecmp(ans, "n") == 0 || strcasecmp(ans, "no") == 0);
    free(ans);

    if (y) return 1;
    if (n) return 0;
    fprintf(stderr, "  %sPlease answer y or n.%s\n",
            kycg_ui_yellow(), kycg_ui_reset());
  }
}

char *kycg_ui_ask(const char *question, const char *def) {
  fprintf(stderr, "%s%s%s", kycg_ui_bold(), question, kycg_ui_reset());
  if (def && *def)
    fprintf(stderr, "\n  %s[%s]%s ", kycg_ui_dim(), def, kycg_ui_reset());
  else
    fputs(" ", stderr);
  fflush(stderr);

  char *ans = read_line();
  if (!ans) return NULL;
  if (!*ans && def) { free(ans); return strdup(def); }
  return ans;
}

/** Render a numbered list, one item per line. */
static void print_items(const char **items, const char **notes, size_t n) {
  int width = 1;
  for (size_t t = n; t >= 10; t /= 10) ++width;

  for (size_t i = 0; i < n; ++i) {
    fprintf(stderr, "  %s%*zu%s  %s",
            kycg_ui_cyan(), width, i + 1, kycg_ui_reset(), items[i]);
    if (notes && notes[i] && *notes[i])
      fprintf(stderr, "  %s%s%s", kycg_ui_dim(), notes[i], kycg_ui_reset());
    fputc('\n', stderr);
  }
}

long kycg_ui_choose(const char *title, const char **items, const char **notes,
                    size_t n) {
  if (!n) return -1;

  fprintf(stderr, "\n%s%s%s\n", kycg_ui_bold(), title, kycg_ui_reset());
  print_items(items, notes, n);

  for (;;) {
    fprintf(stderr, "  %sselect 1-%zu%s ", kycg_ui_dim(), n, kycg_ui_reset());
    fflush(stderr);

    char *ans = read_line();
    if (!ans) { fputc('\n', stderr); return -1; }

    char *end = NULL;
    long v = strtol(ans, &end, 10);
    int clean = (end && end != ans && *end == '\0');
    free(ans);

    if (clean && v >= 1 && v <= (long)n) return v - 1;
    fprintf(stderr, "  %sEnter a number between 1 and %zu.%s\n",
            kycg_ui_yellow(), n, kycg_ui_reset());
  }
}

/**
 * Parse "all" / "none" / "1-5,8,12" into the flag array.
 * Returns 0 on success, -1 if the spec is malformed.
 */
static int parse_selection(const char *spec, int *flags, size_t n) {
  if (strcasecmp(spec, "all") == 0 || strcmp(spec, "*") == 0) {
    for (size_t i = 0; i < n; ++i) flags[i] = 1;
    return 0;
  }
  if (strcasecmp(spec, "none") == 0) {
    for (size_t i = 0; i < n; ++i) flags[i] = 0;
    return 0;
  }

  for (size_t i = 0; i < n; ++i) flags[i] = 0;

  const char *p = spec;
  while (*p) {
    while (*p == ',' || isspace((unsigned char)*p)) ++p;
    if (!*p) break;

    char *end = NULL;
    long lo = strtol(p, &end, 10);
    if (end == p) return -1;
    long hi = lo;
    p = end;

    if (*p == '-') {
      ++p;
      hi = strtol(p, &end, 10);
      if (end == p) return -1;
      p = end;
    }
    if (lo > hi) { long t = lo; lo = hi; hi = t; }
    if (lo < 1 || hi > (long)n) return -1;

    for (long i = lo; i <= hi; ++i) flags[i-1] = 1;

    while (isspace((unsigned char)*p)) ++p;
    if (*p && *p != ',') return -1;
  }
  return 0;
}

int *kycg_ui_multiselect(const char *title, const char **items,
                         const char **notes, size_t n, int default_all) {
  if (!n) return NULL;

  int *flags = calloc(n, sizeof(int));
  if (!flags) return NULL;

  fprintf(stderr, "\n%s%s%s\n", kycg_ui_bold(), title, kycg_ui_reset());
  print_items(items, notes, n);

  for (;;) {
    fprintf(stderr,
            "  %sselect: 'all', 'none', or a list like 1-5,8  [%s]%s ",
            kycg_ui_dim(), default_all ? "all" : "none", kycg_ui_reset());
    fflush(stderr);

    char *ans = read_line();
    if (!ans) { fputc('\n', stderr); free(flags); return NULL; }

    const char *spec = *ans ? ans : (default_all ? "all" : "none");
    int rc = parse_selection(spec, flags, n);
    free(ans);

    if (rc == 0) {
      size_t chosen = 0;
      for (size_t i = 0; i < n; ++i) chosen += (size_t)(flags[i] != 0);
      if (chosen) return flags;
      fprintf(stderr, "  %sNothing selected.%s\n",
              kycg_ui_yellow(), kycg_ui_reset());
      continue;
    }
    fprintf(stderr, "  %sCould not read that. Use 'all', 'none', or "
                    "numbers between 1 and %zu.%s\n",
            kycg_ui_yellow(), n, kycg_ui_reset());
  }
}
