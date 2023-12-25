# Tumor Deconvolution and Immunogenicity Manuscript
## Deconvolution of bulk RNAseq expression profiles and generating integrated scores (iScores)
The ```Rscripts``` folder contains the code that was used to deconvolve the RNAseq gene expression matrix using multiple deconvolution methods and integrate these results as iScores. The ```Signatures``` folder contains gene signatures used for deconvolution and cell types labels used for integration of estimates across different tools into leukocyte iScores and cell types iScores.

### Deconvolution
```run_deconvolutions.R``` implements the following deconvolution tools:

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
[CIBERSORT](https://cibersort.stanford.edu/) is free for non-commerical use only. To run Cibersort, the cibersort R code needs to be requested from the authors of the tool and its full path added to the *tdi_deconvolve_methods.R*.\
[CIBERSORTx](https://cibersortx.stanford.edu) can either be run directly via its website or a docker image can be requested from the authors.

Additionally, the required R packages can be installed from the following links:\
[xCell](http://xcell.ucsf.edu/)\
[Sparse group lasso (decongsl)](https://github.com/drisso/deconsgl)\
[ssGSEA](https://bioconductor.org/packages/release/bioc/html/GSVA.html)\
[MCPCounter](https://github.com/ebecht/MCPcounter)\
[EPIC](https://github.com/GfellerLab/EPIC)\
[Immunedeconv](https://github.com/omnideconv/immunedeconv/)\
[readr](https://cran.r-project.org/web/packages/readr/readme/README.html)

The deconvolution can be executed as follows:

```
Rscript run_deconvolutions.R COHORT_NAME \
/PATH/TO/MIXTURE_FILE \
/PATH/TO/OUTPUTDIR \
/PATH/TO/Signatures \
/PATH/TO/Rscripts \
CANCER_TYPE \
xcell ssgsea mcp sgl epic timer quantiseq
```
MIXTURE_FILE is a gene expression file to be deconvolved, where genes are rows and samples are columns. The first column has [HUGO](https://www.genenames.org/) gene symbols, and the following columns are samples. CANCER_TYPE is a cancer type as defined by the [TCGA](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) e.g., luad. Cancer type is a required parameter for Timer. The last argument takes a list of tools to be used for deconvolution. To include cibersort to this list, full path of the cibersort R code (once obtained from the source) needs to be updated in *tdi_deconvolve_methods.R* at ```source('Path/to/Cibersort/Code/CIBERSORT.R')``` and then the tool name (`cibersort`) can be added to the tool list above.

The outputs from this code are deconvolution estimates from individual tools and the corresponding output filenames contain the name of the tool  (e.g., *xCell_PBM.txt, Cibersort_PBM.txt*) and saved in the output directory defined by the user (OUTPUTDIR).


### Deconvolution Estimate Integration (iScores)

```iScore_calc.R``` converts the deconvolution estimates into normalized scores for each cell type across all samples and integrates them into iScores based on groupings defined in *Signatures/celltype_labels.txt*.

```
Rscript iScore_calc.R COHORT_NAME \
/PATH/TO/DECONVOLVED_DIR \
/PATH/TO/OUTPUTDIR \
/PATH/TO/Signatures
```
DECONVOLVED_DIR is the folder containing deconvolution result files for individual tools (output from the deconvolution step above). The input file names must contain the name of the tool with _PBM suffix as so *xCell_PBM.txt. In each file, samples are rows and celltypes are columns.

The outputs from this code are: a compiled file with normalized estimates from all tool, a file with leukocyte iScores, and a file with individual celltype iScores. These result files are saved in the output directory defined by the user (OUTPUTDIR).


### Troubleshooting
*Note*: quanTIseq may throws the following error which is due to the R version. Please see [here](https://github.com/omnideconv/immunedeconv/issues/146) for more details.
```
Error in `.rowNamesDF<-`(x, value = value) : 
  duplicate 'row.names' are not allowed
Calls: Quantiseq_run ... row.names<- -> row.names<-.data.frame -> .rowNamesDF<-
In addition: Warning message:
non-unique value when setting 'row.names': ‘entry withdrawn’ 
Execution halted
```


