# run_deconvolutions.R - deconvolves the bulk RNAseq FPKM expression profiles
# Copyright (C) 2022  Bhavneet Bhinder (bhb2003@med.cornell.edu)
#Cibersort requires the R code which needs to be requested form the authors of the tool https://cibersortx.stanford.edu

#Define Input file
args <- commandArgs(trailingOnly = TRUE)
cohort=args[1]  #Cohort identifier for use as a  prefix for output files (e.g. test) and labeling
mixture_file=args[2]  #Full file path to tab-delimited FPKM matrix. Rows are genes and columns are samples: first column with HUGO gene symbols is named "Genes"
outDir=args[3]  #Full path to desired output directory
sig_folder=args[4] # Full dir path for the folder containing the signatures
rscript_folder=args[5] # Full dir path to Rscripts which contains tdi_deconvolve_pkg.R
cancer=args[6] #TCGA cancer type, REQUIRED for TIMER
tools_to_run=args[-c(1:6)] #tools_to_run=c('cibersort','xcell','ssgsea','mcp','sgl', 'epic','quantiseq', 'timer')

source(paste0(rscript_folder,"/tdi_deconvolve_methods.R"))

######   Selecting tools to run
timer_available_cancers <- c(
  "acc", "blca", "brca", "cesc", "chol", "coad", "dlbc", "esca", "gbm", "hnsc",
  "kich", "kirc", "kirp", "lgg", "lihc", "luad", "lusc", "meso", "ov", "paad",  "pcpg",
  "prad", "read", "sarc", "skcm", "stad", "tgct", "thca", "thym", "ucec", "ucs", "uvm"
)

if(!cancer %in% timer_available_cancers){
  tools_to_run=tools_to_run[!tools_to_run %in%  c('timer')]
  cat("\nCancer type not found for TIMER analysis, allowed cancer types are: ",timer_available_cancers,"\n" )#;timer_available_cancers
}
cat("\nTools running: ",tools_to_run,"\n\n" )#;timer_available_cancers

######   Define signature file paths
LM22_sig_file<-paste0(sig_folder,"/LM22.txt")
mcp_sig_file=paste0(sig_folder,"/MCP_counter_gene_sig.txt")
gsea_sig_file<-paste0(sig_folder,"/gene_signatures.txt")

######   create output dirs
dir.create(outDir)
output_dir=paste0(outDir,"/deconvolutions")
dir.create(output_dir)

######   Read expression file
exprMatrix <- get_expression_matrix(mixture_file)

######   Run decovolution tools
#1: Cibersort
if ( 'cibersort' %in% tools_to_run) { 
  cat("\nTool :: Cibersort\n\n")
  results.abs<-Cibersort_run(mixture_file, LM22_sig_file)
  write.table(data.frame("SampleID"=row.names(results.abs),results.abs),paste0(output_dir,"/Cibersort_PBM.txt"), sep="\t",quote=F,row.names=F)
  
} 
#2: xCell
if ( 'xcell' %in% tools_to_run) { 
  cat("\nTool :: xCell\n\n")
  xcell.res<-xcell_run(exprMatrix)
  write.table(data.frame("SampleID"=row.names(xcell.res),xcell.res),paste0(output_dir,"/xCell_PBM.txt"),sep="\t",quote=F, row.names=F)
} 
#3: ssGSEA
if ( 'ssgsea' %in% tools_to_run) { 
  cat("\n\nTool :: ssGSEA\n\n")
  ssgsea.res.rnat<-ssgsea_run(exprMatrix, gsea_sig_file)
  write.table(data.frame("SampleID"=row.names(ssgsea.res.rnat),ssgsea.res.rnat),paste0(output_dir,"/ssGSEA_PBM.txt"),sep="\t",quote=F, row.names=F)
} 
#4: MCP counter
if ( 'mcp' %in% tools_to_run) { 
  cat("\nTool :: MCP counter\n\n")
  mcp.rest<-mcp_run(exprMatrix, mcp_sig_file) 
  write.table(data.frame("SampleID"=row.names(mcp.rest),mcp.rest),paste0(output_dir,"/MCP_PBM.txt"),sep="\t",quote=F, row.names=F)
} 
#5: sgl
if ( 'sgl' %in% tools_to_run) { 
  cat("\nTool :: sgl\n\n")
  sgl.res<-sgl_run(exprMatrix, LM22_sig_file)
  write.table(data.frame("SampleID"=row.names(sgl.res),sgl.res),paste0(output_dir,"/SGL_PBM.txt"), sep="\t",quote=F,row.names=F)
} 
#6: EPIC
if ( 'epic' %in% tools_to_run) { 
  cat("\nTool :: EPIC\n\n")
  epic.out<-EPIC_run(exprMatrix) 
  #write.table(data.frame("SampleID"=row.names(epic.out$cellFractions),epic.out$cellFractions),paste0(output_dir,"/EPIC_CellFractions.txt"), sep="\t",quote=F,row.names=F)
  write.table(data.frame("SampleID"=row.names(epic.out$mRNAProportions),epic.out$mRNAProportions),paste0(output_dir,"/EPIC-mRNA_PBM.txt"), sep="\t",quote=F,row.names=F)
  #write.table(data.frame("SampleID"=row.names(epic.out$fit.gof),epic.out$fit.gof),paste0(output_dir,"/EPIC_fit.gof.txt"), sep="\t",quote=F,row.names=F)
}
#7: Timer
if ( 'timer' %in% tools_to_run) { 
  cat("\nTool :: TIMER\n\n")
  dat.tim<-Timer_run (exprMatrix, cancer) 
  write.table(data.frame("SampleID"=row.names(dat.tim),dat.tim),paste0(output_dir,"/TIMER_PBM.txt"), sep="\t",quote=F,row.names=F)
}
#8: Qunatiseq
if ( 'quantiseq' %in% tools_to_run) { 
  cat("\nTool :: Quantiseq\n\n")
  datq<-Quantiseq_run(exprMatrix)
  write.table(data.frame("SampleID"=row.names(datq),datq),paste0(output_dir,"/quantiseq_PBM.txt"), sep="\t",quote=F,row.names=F)
}

