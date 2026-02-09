In transcriptomic data analysis, **alternative splicing (AS) refers to the biological process that generates multiple transcript isoforms from a single gene, while exon usage refers to the quantifiable *result* of this process for a specific exon**. The primary difference lies in the scope of analysis and how the data are modeled: 

| **Feature** | **Exon Usage Analysis** |
| --- | --- |
| **Focus** | Identifies individual exons or genomic regions with statistically different expression levels relative to the parent gene's overall expression. |
| **Goal** | Pinpoints *where* within a gene the regulation is changing (the specific "counting bin"). |
| **Output** | A list of differentially used exonic regions, often with a p-value and fold change. |
| **Tools** | Exon-based methods (e.g., DEXSeq, edgeR, JunctionSeq). |

**Exon Usage (Differential Exon Usage - DEU)**

DEU analysis works by testing each predefined exonic region (or "counting bin") within a gene for differential abundance relative to the rest of the gene across different conditions.

1. **Count reads**: Reads mapping to each exon/counting bin are counted.
2. **Model with GLM**: Generalized linear models (GLMs) are used to model the read counts, accounting for biological variability and library size normalization.
3. **Test for interaction**: A statistical test (e.g., likelihood ratio test) is applied to detect an interaction effect between the condition and the exon counts, indicating that an exon's usage changes disproportionately to the gene's overall expression.

<aside>
👉🏻

Documentation of DEXSeq:
https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#9_Appendix

<aside>
👉🏻

Read this paper to understand exon usage 

</aside>

[https://pmc.ncbi.nlm.nih.gov/articles/PMC3460195/#:~:text=Abstract,Graveley 2010; Grabowski 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3460195/#:~:text=Abstract,Graveley%202010;%20Grabowski%202011)).

</aside>

### **STEP BY STEP TUTORIAL TO RUN FEATURECOUNTS and DEXSeq**

## **Install DEXSeq using BioConductor R**

In R/RStudio run the following:

```python
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DEXSeq")
```

Find python scripts’ location aftter installation:

**1- code lines you'll run on R**

```python
pythonScriptsDir = system.file("python_scripts", package="DEXSeq")
list.files(pythonScriptsDir)
```

**output you'll see:**
[1] "dexseq_count.py"              "dexseq_prepare_annotation.py"

**2- code lines you'll run on R**

```python
system.file( "python_scripts", package="DEXSeq", mustWork=TRUE )
```

output you’ll see somthing like this for your ssytem:

 "/home/sana/R/x86_64-pc-linux-gnu-library/4.5/DEXSeq/python_scripts"

<aside>
👉🏻

***COPY IT and replace below!!!***

</aside>

RUN the below in terminal 

```python
# <- replace with your path found above
DEXSEQ_PY="/usr/local/lib/R/site-library/DEXSeq/python_scripts/DEXSeq_prepare_annotation.py"  
GTF="/path/to/ref/annotation.gtf" #chnage path where your genom gtf file is present

python "$DEXSEQ_PY" \
  "$GTF" \
  "/path/to/ref/annotation_dexseq.gff" #chnage to where you want to keep this gff file
```

## Run featureCounts for **exon-bin** counts

This is the core command, now that Sambamba is done.

Decide strandedness:

- If library is **reverse-stranded** ➜ `s 2`
- If **forward-stranded** ➜ `s 1`
- If **unstranded** ➜ `s 0`

**CHange -s parameter according to the kit and strand**

```bash
featureCounts -T 12 -p --countReadPairs -B -C -s 2 \
  -t exonic_part -g exon_id --extraAttributes gene_id \
  -f -O \
  -a Homo_sapiens.GRCh38.115_dexseq.gff \
  -o exon_counts_dexseq_patients_s2.txt \
  *.markdup.bam

```

The output file is tab-delimited and looks like:

- Column 1: `Geneid` → with `g exon_id`, this will be the **exon bin ID** (best row identifier)
- Column 2–6: `Chr`, `Start`, `End`, `Strand`, `Length` → metadata (useful for sanity checks, not needed for the count matrix)
- Extra attribute column(s): e.g. `gene_id` (because of `-extraAttributes gene_id`) → mapping exon bin → gene
- Then one column per BAM file: sample count columns → **these are the counts matrix you use downstream**

***Once run this command for all 3 datasets (your smaples, and other 2 normal samples), combine them all to make one exon_counts.csv file and use the metadata file that you previously created and used for DESeq2!!!***

```r
suppressPackageStartupMessages({
  library(data.table)
  library(DEXSeq)
  library(AnnotationDbi)
  library(org.Hs.eg.db)  # human gene annotations
})

in_counts <- "exon_counts_dexseq_patients_s2.txt"

# -------------------------
# 1) Read featureCounts table
# -------------------------
# featureCounts output includes comment lines starting with '#'
fc <- fread(in_counts, comment.char = "#", data.table = FALSE, check.names = FALSE)

# Expected first columns: Geneid Chr Start End Strand Length gene_id <samples...>
stopifnot("Geneid" %in% colnames(fc))
stopifnot("gene_id" %in% colnames(fc))

# -------------------------
# 2) Filter ambiguous multi-gene bins (gene_id contains '+')
# -------------------------
fc$gene_id <- as.character(fc$gene_id)
keep <- !grepl("\\+", fc$gene_id) & !is.na(fc$gene_id) & fc$gene_id != ""
fc2 <- fc[keep, , drop = FALSE]

message("Rows total: ", nrow(fc), " | kept (single-gene bins): ", nrow(fc2))

# -------------------------
# 3) Clean Ensembl gene IDs (remove version if present like ENSG... .12)
# -------------------------
gene_ens <- sub("\\..*$", "", fc2$gene_id)

# -------------------------
# 4) Map Ensembl gene IDs -> HGNC symbols
# -------------------------
# org.Hs.eg.db keys for Ensembl are "ENSEMBL"
sym <- mapIds(
  org.Hs.eg.db,
  keys = unique(gene_ens),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

fc2$gene_ens <- gene_ens
fc2$gene_symbol <- unname(sym[gene_ens])

# If you want to keep unmapped symbols visible:
fc2$gene_symbol[is.na(fc2$gene_symbol) | fc2$gene_symbol == ""] <- "NA"

# -------------------------
# 5) Prepare count matrix and exon->gene annotation table
# -------------------------
# Identify sample columns: everything after 'gene_id' (and not our added cols)
# We'll locate the position of gene_id in the original featureCounts output.
gid_col <- which(colnames(fc2) == "gene_id")
sample_cols <- (gid_col + 1):ncol(fc2)

# If you appended columns, make sure sample_cols excludes them:
# We'll recompute sample_cols as columns that are numeric and not annotation fields.
ann_fields <- c("Geneid","Chr","Start","End","Strand","Length","gene_id","gene_ens","gene_symbol")
sample_cols <- which(!(colnames(fc2) %in% ann_fields))

# Count matrix (rows=exon bins, cols=samples)
count_mat <- as.matrix(fc2[, sample_cols, drop = FALSE])
mode(count_mat) <- "integer"
rownames(count_mat) <- fc2$Geneid

# Sample names as in header
sample_names <- colnames(fc2)[sample_cols]

# Exon annotation table (helpful for reporting)
exon_annot <- fc2[, c("Geneid","gene_ens","gene_symbol","Chr","Start","End","Strand"), drop = FALSE]
colnames(exon_annot)[1] <- "exon_id"

# Write out reusable files
write.table(
  cbind(exon_id = rownames(count_mat), count_mat),
  file = "dexseq_counts_matrix.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  exon_annot,
  file = "dexseq_exon_annotation.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -------------------------
# 6) DEXSeq dataset + differential exon usage
# -------------------------
# You must provide sample metadata with 'condition' (and optionally batch, etc.)
# Create samples.csv with columns: sample, condition
coldata <- read.csv("samples.csv", stringsAsFactors = TRUE)
stopifnot(all(sample_names %in% coldata$sample))

# reorder to match count matrix columns
rownames(coldata) <- coldata$sample
coldata <- coldata[sample_names, , drop = FALSE]

# groupID: gene ID per exon bin (required)
group_id <- fc2$gene_ens
names(group_id) <- fc2$Geneid

dxd <- DEXSeqDataSet(
  countData  = count_mat,
  sampleData = coldata,
  design     = ~ sample + exon + condition:exon,
  featureID  = rownames(count_mat),
  groupID    = unname(group_id[rownames(count_mat)])
)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

res <- DEXSeqResults(dxd)

# Attach gene symbols for easier reading
res_df <- as.data.frame(res)
res_df$gene_ens <- group_id[res_df$featureID]
res_df$gene_symbol <- unname(sym[res_df$gene_ens])

write.csv(res_df, "DEXSeq_results_with_symbols.csv", row.names = FALSE)

message("Done. Outputs: dexseq_counts_matrix.tsv, dexseq_exon_annotation.tsv, DEXSeq_results_with_symbols.csv")

# -------------------------
# 7) Visualizations
# -------------------------
dir.create("plots", showWarnings = FALSE)

# 7.1 Library sizes (raw)
lib_sizes <- colSums(counts(dxd))
pdf("plots/QC_library_sizes.pdf", width = 8, height = 4)
barplot(lib_sizes, las = 2, main = "Library sizes (sum of exon-bin counts)", ylab = "Reads (exon bins)")
dev.off()

# 7.2 Size factors
sf <- sizeFactors(dxd)
pdf("plots/QC_size_factors.pdf", width = 8, height = 4)
barplot(sf, las = 2, main = "DEXSeq size factors", ylab = "Size factor")
abline(h = 1, lty = 2)
dev.off()

# 7.3 Dispersion plot (built-in)
pdf("plots/QC_dispersion.pdf", width = 6, height = 5)
plotDispEsts(dxd, main = "Dispersion estimates")
dev.off()

# 7.4 P-value histogram
pdf("plots/Results_pvalue_hist.pdf", width = 6, height = 4)
hist(res$pvalue, breaks = 50, main = "P-value distribution (exon bins)", xlab = "p-value")
dev.off()

# 7.5 MA-style plot using exon fold-change (DEXSeq provides log2fold)
# We'll use res$log2fold and res$baseMean (if available) or mean normalized counts proxy
res_df <- as.data.frame(res)
if (!("baseMean" %in% colnames(res_df))) {
  # fallback: mean normalized counts per exon bin
  norm_cts <- counts(dxd, normalized = TRUE)
  res_df$baseMean <- rowMeans(norm_cts)
}

pdf("plots/Results_MA_exon_usage.pdf", width = 6, height = 5)
with(res_df, {
  plot(log10(baseMean + 1), log2fold,
       pch = 16, cex = 0.3,
       xlab = "log10(mean normalized count + 1)",
       ylab = "log2 exon usage fold change",
       main = "MA-style: exon usage changes")
  abline(h = 0, lty = 2)
})
dev.off()

# 7.6 Top genes: classic DEXSeq exon usage plots
# pick top N genes by minimal padj across their exon bins
res_df$featureID <- rownames(res_df)
# groupID is stored in the DEXSeqDataSet; we can pull it from dxd
gid <- DEXSeq::groupID(dxd)
names(gid) <- rownames(dxd)
res_df$gene_ens <- gid[res_df$featureID]

# drop NA padj
res_df2 <- res_df[!is.na(res_df$padj), ]
# gene-level ranking
top_genes <- head(names(sort(tapply(res_df2$padj, res_df2$gene_ens, min), decreasing = FALSE)), 12)

pdf("plots/TopGenes_DEXSeq_plots.pdf", width = 7, height = 5)
for (g in top_genes) {
  # DEXSeq plot wants gene identifier = groupID (here Ensembl gene id)
  plotDEXSeq(res, g, legend = TRUE, displayTranscripts = FALSE, main = paste0("DEXSeq: ", g))
}
dev.off()

# 7.7 Heatmap: exon log2 fold-changes for top genes (significant bins only)
sig <- res_df2[res_df2$padj < 0.05, ]
if (nrow(sig) > 0) {
  top_sig_genes <- head(names(sort(tapply(sig$padj, sig$gene_ens, min), decreasing = FALSE)), 25)
  sig_sub <- sig[sig$gene_ens %in% top_sig_genes, ]

  # Build a matrix: rows = exon bins, value = log2fold
  # For readability keep top K exon bins by padj
  sig_sub <- sig_sub[order(sig_sub$padj), ]
  sig_sub <- head(sig_sub, 200)

  mat <- matrix(sig_sub$log2fold, ncol = 1)
  rownames(mat) <- sig_sub$featureID
  colnames(mat) <- "log2FC_exonUsage"

  # attach gene symbol if you have sym mapping from earlier
  # (optional annotation in rownames)
  if ("gene_symbol" %in% colnames(res_df)) {
    # nothing special here; we can paste gene on rownames if you want
  }

  pdf("plots/Heatmap_top_sig_exon_bins_log2FC.pdf", width = 6, height = 10)
  pheatmap(mat, cluster_rows = TRUE, cluster_cols = FALSE,
           main = "Top significant exon bins (log2 fold-change)")
  dev.off()
} else {
  message("No significant exon bins at padj < 0.05; skipping heatmap.")
}

message("Plots saved in ./plots/")

```
