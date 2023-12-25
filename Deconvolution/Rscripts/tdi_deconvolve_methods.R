# DECONVOLUTION R script v1 (last updated 12-16-2021)

#Load libraries and packages

#source('Path/to/Cobersort/Code/CIBERSORT.R') #add the path to cibersort R code obatined from https://cibersortx.stanford.edu

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(lsgl))
suppressPackageStartupMessages(library(xCell))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(MCPcounter))
suppressPackageStartupMessages(library(immunedeconv))
suppressPackageStartupMessages(library(EPIC))


## FUNCTIONS

get_expression_matrix<-function(data.file){
  exprMatrix <- read.table(data.file, header=TRUE,sep="\t", row.names=1, check.names = F)
  exprMatrix.noZero<- exprMatrix[rowSums(exprMatrix == 0) < 0.90*ncol(exprMatrix), ] 
  exprMatrix.noZeroVar<-exprMatrix.noZero[apply(exprMatrix.noZero, 1, var) > 0,]
  return(exprMatrix.noZeroVar)
}

Cibersort_run<-function(x, sig_file){
  results.abs <- CIBERSORT(sig_file, x, perm=100 , QN=FALSE, absolute=TRUE, abs_method='sig.score')
  return(results.abs)
}

xcell_run<-function(y){
  xcell.res.ipm<-xCellAnalysis(y)
  return(t(xcell.res.ipm))
}

ssgsea_run<-function(y, z){
  immune_sig_custom<-create_ssgsea_sig(z)
  ssgsea.res.rna<-gsva(as.matrix(y), gset.idx.list=immune_sig_custom, method="ssgsea",min.sz=1, max.sz=Inf)
  return(t(ssgsea.res.rna))
}

mcp_run<-function(y, m){
  mcp.sig<-read.table(m,sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
  mcp.res<-MCPcounter.estimate(y,featuresType="HUGO_symbols", genes=mcp.sig)
  return(t(mcp.res))
}

EPIC_run<-function(bulkSamplesMatrix){
  epic.out <- EPIC(bulk = bulkSamplesMatrix)
  return(epic.out)
}

Timer_run<-function(bulkSamplesMatrix, cancer){
  res.timer<-deconvolute(bulkSamplesMatrix, "timer",indications=c(rep(cancer,ncol(bulkSamplesMatrix))))
  dat.timer<-as.data.frame(res.timer); row.names(dat.timer)<-dat.timer[,1]; dat.timer<-t(dat.timer[,-1])
  return(dat.timer)
}

Quantiseq_run<-function(bulkSamplesMatrix){
  res.quantiseq<-deconvolute(bulkSamplesMatrix, "quantiseq", tumor = TRUE)
  datq<-as.data.frame(res.quantiseq); row.names(datq)<-datq[,1]; datq<-t(datq[,-1])
  return(datq)
}

sgl_run<-function(mat, sig_file){
  lm22 <- read.table(sig_file, header=TRUE, row.names=1, sep='\t')
  genes <- intersect(rownames(lm22), rownames(mat))
  tumor=log2(mat[genes,]+1)
  X <- log(as.matrix(lm22[genes,]))
  tumor<-tumor[order(match(rownames(X),rownames(tumor))),]
  sgl.res.p <- get_coefs(tumor, X)
  return(t(sgl.res.p))
}

get_coefs <- function(tumor,X, cores=6) {
  groups <- c(1, 1, 1, rep(2, 7), 3, 3, rep(4, 4), 5, 5, 6, 6, 7, 7) # length=22
  fit <- mclapply(as.data.frame(tumor), function(y) {
    fit_cv <- lsgl::cv(X, y, intercept = FALSE, grouping = groups,  alpha = 0.5, 
                       lambda = 0.1)
    fit3 <- lsgl::fit(X, y, intercept = FALSE, alpha = 0.5, lambda = 0.1)
    fit3$beta[[best_model(fit_cv)]]
  }, mc.cores = cores)
  coefs <- sapply(fit, as.vector)
  rownames(coefs) <- colnames(X)
  return(coefs)
}

read_tcga <- function(mfile, g) {
  data <- read.table(mfile, header=TRUE,sep="\t", as.is=TRUE) 
  duplicates <- names(which(table(data$Genes)>1))
  mat <- as.matrix(data[!(data$Genes %in% duplicates),-1])
  rownames(mat) <- data[!(data$Genes %in% duplicates),]$Genes
  genes <- intersect(g, rownames(mat))
  return(log2(mat[genes,]+1))
}

create_ssgsea_sig <-function(sig_file){
  gsea_imm<-read.table(sig_file,sep="\t",h=T)
  gsea_imm$Sig_name<-paste0(gsea_imm$Celltypes,"|" ,gsea_imm$Source)
  immsig<-aggregate(Genes ~ Sig_name, data=gsea_imm, FUN=function(x){ as.vector(x)});
  ssgsea_signature<-setNames(as.list(immsig$Genes),immsig$Sig_name)
  return(ssgsea_signature)
}




  
