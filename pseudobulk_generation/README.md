# Tumor-Deconvolution-and-Immunogenicity-Manuscript
## Pseudobulk Random Generator
Generates a matrix of 1000 pseudobulk assemblies by randomly selecting the gene expression profiles of 10% of the total population of single cells available in the input data matrix (e.g, test_raw_counts.txt). 

### SYSTEM REQUIREMENTS:
- CentOS Linux 7 (Core)
- x86_64-redhat-linux-gnu (64-bit)
- R version 3.6.0 (2019-04-26)
- R base packages: stats, graphics, grDevices, utils, datasets, methods, base     

### USAGE:
```
Rscript pseudobulk_random_generator.r testout /path/to/test_raw_counts.txt /path/to/test_svm_predictions.tsv /path/to/PRG_output`
```

In this case, output will be written to the current directory ($PWD) with the "testout" prefix

### INPUTS:
- Raw read counts - tab-delimited matrix of raw single-cell sequence read counts (e.g., test_raw_counts.txt)
  - Rows: HUGO gene symbols (e.g., FAM87B)
  - Columns: Cell identifier (format: <sample identifier>_<cell barcode>.<rep>; e.g., e01_AAACGCTTCAACCCGG.1)
  - Values: Raw read counts

- Cell type predictions - single-column matrix of cell type predictions for each cell in the dataset as determined by a Support Vector Model (e.g. test_svm_predictions.tsv)
  - Rows: Cell identifier
  - Column: "label"
  - Values: Cell type identifier (e.g., bmem)

#### Cell Type Identifier Key for Test Data:
- Memory B-Cell (bmem)
- Naive B-Cell (bnv)
- M2 Macrophage (m2)
- Macrophage (mac)
- Monocyte (mon)
- Mast (mst)
- Natural Killer (nk)
- Plasma (plsm)
- CD4+ Memory T-Cell (t4mem)
- CD8+ T-Cell (t8)

### OUTPUTS:
NOTE: Actual output values will vary from run to run as each pseudobulk column is assembled from a random sample of all input cells.

- Pseudobulk matrix - a 1000-column tab-delimited matrix of raw read counts summed across the cells randomly selected for each pseudobulk (e.g. test.pseudobulk.txt)
  - Rows: HUGO gene symbols
  - Columns: Pseudobulk identifier (format: <cohort identifier>_pbulk_<pseudobulk set>; e.g., test_pbulk_123)
  - Values: Summed raw read counts

- Pseudobulk register - a 1000-column tab-delimited matrix of the cell-type-annotated cell identifiers randomly selected for each pseudobulk set (e.g. test.pseudobulk.labeled.txt)
  - Rows: Numbered 1:n, where n=#cells selected for each pseudobulk set (e.g., n=417)
  - Columns: Pseudobulk identifier
  - Values: Cell-type-annotated cell identifier (format: <cell identifier>_<cell type identifier>; e.g., e01_AAACGCTTCAACCCGG.1_bmem)

- Cell type population fraction matrix - a 1000-column tab-delimited matrix of pseudobulk population fractions represented by each of the available cell types (e.g. test.pseudobulk.fxn.txt)
  - Rows: Cell type identifier (one row per cell type; e.g., bmem)
  - Columns: Pseudobulk sample identifier
  - Values: Fraction of the pseudobulk population represented by each cell type

### EXAMPLE RUNTIME:
```
real     2m0.017s
user     1m52.316s
sys      0m6.426s
```

