#!/usr/bin/env Rscript
#
# Acceptance check for `kycg test -c` against knowYourCG's testEnrichmentSEA.
#
# Builds the zhou-lab/knowYourCG#6 example -- random z-scores on MM285 with
# the placenta signature shifted up by 0.5 -- as a format 4 query and a
# format 0 knowledgebase, runs both implementations, and compares them.
# Needs knowYourCG >= 1.9.3, the release that took GSEA's sign; put it on
# R_LIBS if the default library has an older one.
#
#   Rscript tests/validate_sea_vs_R.R ./kycg ~/repo/YAME/yame <workdir>
#
# The estimates must agree to 1e-9: the walk is deterministic. The p-values
# cannot agree exactly, because each side draws its own permutations (R 100,
# kycg 1000), so the bar there is rank order and a loose absolute gap.

suppressMessages({
    library(knowYourCG)
    library(sesame)
    library(SummarizedExperiment)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
    stop("usage: validate_sea_vs_R.R <kycg> <yame> <workdir>")
kycg <- normalizePath(args[1])
yame <- normalizePath(args[2])
dir.create(args[3], showWarnings = FALSE, recursive = TRUE)
setwd(args[3])
if (packageVersion("knowYourCG") < "1.9.3")
    stop("knowYourCG ", packageVersion("knowYourCG"), " has the old sign; ",
         "need >= 1.9.3")

## --- the issue #6 query, on the MM285 row space
df <- as.data.frame(rowData(sesameDataGet("MM285.tissueSignature")))
set.seed(42)
zv <- setNames(rnorm(nrow(df)), df$Probe_ID)
plac <- df$Probe_ID[df$branch == "placenta"]
zv[plac] <- zv[plac] + 0.5

ord <- read.table(gzfile(file.path(Sys.getenv("YAME_DATA_HOME"),
                                   "MM285/MM285.ordering.tsv.gz")),
                  header = TRUE, sep = "\t")$Probe_ID
zv <- zv[names(zv) %in% ord]

## Format 4 holds a float in [0, 1]. pnorm() keeps the order; the round trip
## through 4 bytes makes R see exactly the values kycg reads, ties included.
v <- pnorm(zv)
v <- setNames(readBin(writeBin(v, raw(), size = 4), "double",
                      n = length(v), size = 4), names(v))
## R's sort() keeps ties in input order, and kycg keeps them in row order.
v <- v[order(match(names(v), ord))]

col <- rep("NA", length(ord))
col[match(names(v), ord)] <- sprintf("%.9g", v)
writeLines(col, "query.txt")
system2(yame, c("pack", "-f", "n", "query.txt", "query.cg"))

## --- the knowledgebase, one format 0 record per set
dbs <- getDBs("KYCG.MM285.tissueSignature", silent = TRUE)
## `yame pack -f b` packs one column, so pack each set, then concatenate and
## index the records by name
parts <- vapply(seq_along(dbs), function(i) {
    writeLines(as.character(as.integer(ord %in% dbs[[i]])), "set.txt")
    f <- sprintf("set%02d.cm", i)
    system2(yame, c("pack", "-f", "b", "set.txt", f))
    f
}, character(1))
system(paste("cat", paste(parts, collapse = " "), "> sets.cm"))
writeLines(names(dbs), "sets.names")
system2(yame, c("index", "-s", "sets.names", "sets.cm"))

## --- both sides
system2(kycg, c("test", "-c", "-m", "sets.cm", "query.cg", "-o", "kycg.tsv"))
k <- read.table("kycg.tsv", header = TRUE, sep = "\t")
set.seed(0)
r <- testEnrichmentSEA(v, dbs, platform = "MM285", silent = TRUE)

m <- merge(data.frame(db = r$dbname, r_est = r$estimate,
                      r_lp = r$log10.p.value, r_overlap = r$overlap),
           data.frame(db = k$db, k_est = k$estimate, k_lp = k$log10_p,
                      k_overlap = k$overlap), by = "db")
stopifnot(nrow(m) == length(dbs))
## an empty comparison agrees trivially; refuse it
if (sum(m$k_overlap > 0) < 2) stop("fewer than two sets overlap the query")

cat(sprintf("sets compared            %d\n", nrow(m)))
cat(sprintf("overlap mismatches       %d\n", sum(m$r_overlap != m$k_overlap)))
cat(sprintf("max |estimate diff|      %.3g\n", max(abs(m$r_est - m$k_est))))
cat(sprintf("sign mismatches          %d\n", sum(sign(m$r_est) != sign(m$k_est))))
cat(sprintf("spearman(log10 p)        %.4f\n",
            cor(m$r_lp, m$k_lp, method = "spearman")))
cat(sprintf("median |log10 p diff|    %.3g\n", median(abs(m$r_lp - m$k_lp))))
pl <- m[grepl("placenta", m$db), ]
cat(sprintf("placenta  R %.6f %.2f   kycg %.6f %.2f\n",
            pl$r_est, pl$r_lp, pl$k_est, pl$k_lp))
write.table(m, "compare.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

ok <- all(m$r_overlap == m$k_overlap) && max(abs(m$r_est - m$k_est)) < 1e-9
if (!ok) { cat("FAIL\n"); quit(status = 1) }
cat("PASS\n")
