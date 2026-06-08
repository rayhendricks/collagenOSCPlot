#!/usr/bin/env Rscript
# DESeq2 median-of-ratios normalization of the GSE130811 raw count matrix.
suppressMessages({ library(DESeq2) })

infile  <- "GSE130811_expr.tab.gz"
# read.delim handles the R rownames offset: header has 1 fewer field than data rows,
# so the first data column (WBGene IDs) becomes rownames.
df <- read.delim(gzfile(infile), header = TRUE, row.names = 1, check.names = FALSE)
cat("Loaded matrix:", nrow(df), "genes x", ncol(df), "cols\n")
cat("First cols:", paste(head(colnames(df), 4), collapse = ", "), "...\n")

# Drop the gene-length 'width' column; everything else is an integer count column.
count_cols <- setdiff(colnames(df), "width")
counts <- as.matrix(df[, count_cols])
storage.mode(counts) <- "integer"
cat("Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

coldata <- data.frame(row.names = colnames(counts), sample = colnames(counts))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~1)
dds <- estimateSizeFactors(dds)

cat("\n=== Size factors ===\n")
print(round(sizeFactors(dds), 3))

norm <- counts(dds, normalized = TRUE)
out <- data.frame(wbgene = rownames(norm), round(norm, 4), check.names = FALSE)
write.table(out, file = "normalized_counts.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

sf <- data.frame(sample = names(sizeFactors(dds)), size_factor = sizeFactors(dds))
write.table(sf, file = "size_factors.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWrote normalized_counts.tsv and size_factors.tsv\n")
