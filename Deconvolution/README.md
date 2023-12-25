# Tumor Deconvolution and Immunogenicity Manuscript
## Deconvolution of bulk RNAseq expression profiles and generating integrated scores (iScores)
The ```Rscripts``` folder contains the code that was used to deconvolve the RNAseq gene expression matrix using multiple deconvolution methods and integrate these results as iScores. The ```Signatures``` folder contains gene signatures used for deconvolution and cell types labels used for integration of estimates across different tools into leukocyte iScores and cell types iScores.

### Deconvolution
```run_deconvolutions.R``` takes as input a gene experession file where genes are rows and samples are columns. The first column has [HUGO](https://www.genenames.org/) gene symbols. This script implements the following tools:

```
cibersort
xCell
decongsl
ssgsea
mcp_counter
epic
quantiseq
timer
```
[CIBERSORT](https://cibersort.stanford.edu/) is free for non-commerical use only. To run Cibersort, the cibersort R code needs to be requested from the authors of the tool and its full path added to the *tdi_deconvolve_methods.R*. Similarly, [CIBERSORTx](https://cibersortx.stanford.edu) can either be run directly via their website or a docker image can be requested from the authors.

Additional required R packages can be installed from the following links:

######
[xCell](http://xcell.ucsf.edu/)\
[Sparse group lasso (decongsl)](https://github.com/drisso/deconsgl)\
[ssGSEA](https://bioconductor.org/packages/release/bioc/html/GSVA.html)\
[MCPCounter](https://github.com/ebecht/MCPcounter)\
[EPIC](https://github.com/GfellerLab/EPIC)\
[Immunedeconv](https://github.com/omnideconv/immunedeconv/)

### Deconvolution Estimate Integration

```iScore_calc.R``` converts the deconvolution results into normalized scores for each cell type across all samples deconvolved and integrates them into iScores based on groupings defined in *Signatures/celltype_labels.txt*.


