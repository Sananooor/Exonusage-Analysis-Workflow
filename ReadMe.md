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

***Once run this command for all 3 datasets, combine them all to make one exon_counts.csv file and use the metadata file that you previously created and used for DESeq2!!!***


And using the DEXSeq.r file perform the analysis.

