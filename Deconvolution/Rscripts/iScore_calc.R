# iScore_calc.R - Normalize and calculated iScores from deconvolved estimates
# Copyright (C) 2022  Bhavneet Bhinder (bhb2003@med.cornell.edu)

#Define Input file
args <- commandArgs(trailingOnly = TRUE)

cohort=args[1]  #Cohort identifier for use as a  prefix for output files (e.g. test) and labeling
decon_dir=args[2]  #Full file path to the dir containing deconvolution tool outputs. The files in this folder are names as [Tool]_PBM.txt i.e., CIBERSORTx_PBM.txt,Cibersort_PBM.txt,EPIC-mRNA_PBM.txt,MCP_PBM.txt,SGL_PBM.txt,TIMER_PBM.txt,quantiseq_PBM.txt,ssGSEA_PBM.txt,xCell_PBM.txt
outDir=args[3]  #Full path to desired output directory

sig_folder="Signatures" #Full dir path for the folder containing the signatures
cell_groupings_file=paste0(sig_folder,"/celltype_labels.txt")


#create output dirs
dir.create(outDir)


###### FUNCTIONS & Libraries ######
#calculate z-score for a given vector
zscore<- function(x){
  z<- (x - mean(x,na.rm = TRUE)) / sd(x,na.rm = TRUE)
  return(z)
}
#reads files and converts into column-wise z-scores i.e, z-scores for each cell type across all samples
read_decon_files <- function(dfile) {
  data <- read.table(dfile, header=TRUE,sep="\t", row.names=1) #Read File
  tool<-gsub(".*/|_PBM.*","",dfile) #Extract Tool Name
  if(tool == "Cibersort"){  data<-data[,1:22]}
  if(tool == "CIBERSORTx"){  data<-data[,1:22]}
  if(tool == "xCell"){  data<-data[,1:64]}
  data_mix<- data[, colSums(data != 0) > 0] #remove columns with zero values
  return(list(d1=data_mix, d2=tool))
}
#apply z-score
apply_zscore<- function(y, tool){
  data_zscore <- apply(y, 2, zscore) # Calc z-scores
  colnames(data_zscore) <- paste0(colnames(data_zscore),"_",tool)
  return(data_zscore)
}


###### Convert deconvolution scores to normalized scores ######
files <- list.files(decon_dir)

decon.z <- data.frame()
for (i in files){
  print(gsub("_PBM.txt","",i))
  dat.mat<-read_decon_files(paste0(decon_dir,"/",i))
  dat<-apply_zscore(dat.mat$d1, dat.mat$d2)
  if (nrow(decon.z) == 0) {
    decon.z<-dat
  } else {
    decon.z<-merge(decon.z, dat, by=0,all=T)
    row.names(decon.z)<-decon.z[,1]; decon.z<-decon.z[,-1]
  }
}

dim(decon.z); decon.z[1:6,1:6]
  
write.table(data.frame("SampleID"=row.names(decon.z),decon.z),paste0(outDir,"/All_tools_zscore_",cohort,".txt"), sep="\t",quote=F,row.names=F)


###################################################################################################################################
######### Aggregate into integrated scores or iScores ######
decoz<-read.table(paste0(outDir,"/All_tools_zscore_",cohort,".txt"),header=TRUE,sep="\t", row.names=1)
cell_gps<-read.table(cell_groupings_file,header=TRUE,sep="\t", row.names=1)

d<-merge(cell_gps[,1:2],t(decoz),by=0,all=F); row.names(d)<-d[,1]; d<-d[,-1]

#1. Leukocyte iScores
dleuk<-droplevels(d[d$CellGroup %in% c("Leukocyte"),])
dc<-aggregate(dleuk[, -c(1:2)], list(dleuk$CellGroup), mean,na.rm = TRUE)
row.names(dc)<- dc[,1]; dc<- t(dc[,-1]); colnames(dc)<-"Leukocyte_iScore"

write.table(data.frame("SampleID"=row.names(dc),dc),paste0(outDir,"/Leukocyte_iScore_",cohort,".txt"), sep="\t",quote=F,row.names=F)


#2. Celltype iScores
dcell<-aggregate(d[, -c(1:2)], list(d$CellType), mean,na.rm = TRUE)
colnames(dcell)[1]<-"CellTypes"
write.table(dcell,paste0(outDir,"/CellType_iScore_",cohort,".txt"), sep="\t",quote=F,row.names=F)

