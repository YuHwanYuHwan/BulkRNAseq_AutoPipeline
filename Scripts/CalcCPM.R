# CalcCPM.R <count_matrix.tsv> <out_logcpm.tsv>
# TMM-normalized CPM via edgeR, written as log2(CPM + 1). Reads the matrix from ReadCount.sh.
#
# The log is what the values get used for: PCA, clustering and heatmaps all need it, or a
# handful of highly expressed genes take the whole variance. The +1 keeps zeros at zero.
# edgeR's own cpm(log = TRUE, prior.count = 2) is the better choice where low counts matter,
# since it scales what it adds by library size instead of adding 1 to everything.
#
# Differential expression does not read this file. edgeR and DESeq2 model raw counts and are
# given count_matrix.tsv; limma-voom makes its own log-CPM.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)

suppressPackageStartupMessages(library(edgeR))

x <- read.delim(args[1], row.names = 1, check.names = FALSE)
stopifnot(ncol(x) > 0, all(sapply(x, is.numeric)))

y <- calcNormFactors(DGEList(counts = x))
write.table(log2(cpm(y, normalized.lib.sizes = TRUE) + 1), args[2],
            sep = "\t", quote = FALSE, col.names = NA)

cat(sprintf("[CPM ] %d genes x %d samples, log2(CPM+1) -> %s\n", nrow(x), ncol(x), args[2]))
