#!/usr/bin/env Rscript
#
# End-to-end acceptance check for `kycg test`.
#
# Reads a kycg TSV, recomputes every statistic with the same R calls
# knowYourCG uses (phyper / p.adjust), and reports the worst disagreement.
# Agreement on log10.p.value at ~1e-10 is the acceptance bar for Phase 1.
#
#   ./kycg test -m chromhmm.cm onecell.cg > res.tsv
#   Rscript tests/validate_vs_R.R res.tsv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: validate_vs_R.R <kycg-output.tsv>")

res <- read.table(args[1], header = TRUE, sep = "\t",
                  na.strings = "NA", stringsAsFactors = FALSE)

nU <- res$nU; nQ <- res$nQ; nD <- res$nD; nDQ <- res$overlap

nDmQ  <- nD - nDQ
nQmD  <- nQ - nDQ
nUmDQ <- nU - nQ - nD + nDQ

## --- p-value: knowYourCG's calculate_fisher_pvalue(), alternative="greater"
m <- nDQ + nQmD
n <- nUmDQ + nDmQ
k <- nDmQ + nDQ
r_log10p <- phyper(nDQ - 1, m, n, k, lower.tail = FALSE, log.p = TRUE) / log(10)

## --- effect sizes: calculate_odds_ratio() + calculate_effect_sizes()
odds <- (as.numeric(nDQ) * as.numeric(nUmDQ)) /
        (as.numeric(nQmD) * as.numeric(nDmQ))
odds[is.infinite(odds)] <- .Machine$double.xmax
odds[odds == 0] <- .Machine$double.xmin
odds[is.nan(odds)] <- NA_real_
r_estimate <- log2(odds)

r_jaccard <- as.numeric(nDQ) / as.numeric(nD + nQmD)
r_mcc <- (as.numeric(nDQ) * as.numeric(nUmDQ) -
          as.numeric(nQmD) * as.numeric(nDmQ)) /
         sqrt(as.numeric(nD) * (nU - nD) * nQ * (nU - nQ))
r_overlap <- as.numeric(nDQ) / as.numeric(pmin(nD, nQ))
r_npmi <- (log2(nD) + log2(nQ) - 2 * log2(nU)) / (log2(nDQ) - log2(nU)) - 1
r_dice <- (2 * as.numeric(nDQ)) / as.numeric(nD + nQ)

## --- FDR: p.adjust within group, as set_FDR(mtc_by_group=TRUE) does.
## Done on log10 p directly so the comparison is not itself destroyed by the
## underflow that motivated kycg's log-space BH.
bh_log10 <- function(lp) {
    o <- order(lp)
    m <- length(lp)
    ## p_(i) * m / i, in log10; then step up (running min from the largest
    ## p-value down), clamped at log10(1) = 0.
    v <- lp[o] + log10(m) - log10(seq_len(m))
    acc <- 0
    out <- numeric(m)
    for (i in m:1) { acc <- min(acc, v[i]); out[i] <- acc }
    q <- numeric(m)
    q[o] <- out
    q
}
## The stratum is (query_file, query, group): knowYourCG corrects one
## testEnrichment() frame at a time and that frame covers a single query, then
## p.adjust runs within group inside it.
stratum <- paste(res$query_file, res$query, res$group, sep = "\r")
r_fdr <- rep(NA_real_, nrow(res))
for (g in unique(stratum)) {
    i <- which(stratum == g)
    r_fdr[i] <- bh_log10(r_log10p[i])
}

## --- compare
cmp <- function(label, got, want) {
    ok <- (is.na(got) & is.na(want)) |
          (!is.na(got) & !is.na(want) &
           abs(got - want) <= 1e-9 + 1e-10 * abs(want))
    worst <- suppressWarnings(max(abs(got - want), na.rm = TRUE))
    if (!is.finite(worst)) worst <- 0
    cat(sprintf("  %-14s %s   worst abs diff %.3g\n",
                label, if (all(ok)) "OK  " else "FAIL", worst))
    all(ok)
}

cat(sprintf("Comparing %d rows against R\n", nrow(res)))
ok <- c(
    cmp("log10_p",     res$log10_p,      r_log10p),
    cmp("estimate",    res$estimate,     r_estimate),
    cmp("log10 fdr",   -res$neglog10_fdr, r_fdr),
    cmp("cf_jaccard",  res$cf_jaccard,   r_jaccard),
    cmp("cf_mcc",      res$cf_mcc,       r_mcc),
    cmp("cf_overlap",  res$cf_overlap,   r_overlap),
    cmp("cf_npmi",     res$cf_npmi,      r_npmi),
    cmp("cf_dice",     res$cf_dice,      r_dice)
)

if (all(ok)) {
    cat("PASS: kycg agrees with R on every column.\n")
    quit(status = 0)
} else {
    cat("FAIL: at least one column disagrees.\n")
    quit(status = 1)
}
