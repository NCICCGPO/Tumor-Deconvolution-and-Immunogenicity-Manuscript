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

The outputs from this code are deconvolution estimates from individual tools and the corresponding output filenames contain the name of the tool  (e.g., *xCell_PBM.txt, Cibersort_PBM.txt*) and placed in the output directory defined by the user.


### Deconvolution Estimate Integration

```iScore_calc.R``` converts the deconvolution estimates into normalized scores for each cell type across all samples and integrates them into iScores based on groupings defined in *Signatures/celltype_labels.txt*. The input is a folder containing deconvolution result files for individual tools. Ideally, the input file names contain the name of the tool e.g., *Cibersort_PBM.txt*. In each file, samples are rows and celltypes are columns. 

```
Rscript iScore_calc.R COHORT_NAME /PATH/TO/DECONVOLUTION_OUTPUTS_FOLDER /PATH/TO/OUTPUTDIR /PATH/TO/Signatures
```

The outputs from this code are: a compiled file with normalized estimates from all tool, a file with leukocyt iScores, and a file with individual celltype iScores.




