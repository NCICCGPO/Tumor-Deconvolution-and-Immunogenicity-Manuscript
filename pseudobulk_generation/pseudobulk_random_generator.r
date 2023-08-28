# pseudobulk_random_generator.r - Constructs 1000 pseudobulk samples through 
# random selection of 10% of cells from the total input population.
# Copyright (C) 2021  Kami E. Chiotti (chiotti@ohsu.edu)
#
# Usage: Rscript pseudobulk_random_generator.r test $PWD/test_raw_counts.txt $PWD/test_svm_predictions.tsv $PWD
# In the test case, output will be written to the current directory ($PWD) with file names and identifiers prefixed "test".

args=commandArgs(TRUE)

cohort=args[1]  #Cohort identifier for use as a  prefix for output files (e.g. test) and labeling
rawFile=args[2]  #Full file path to tab-delimited matrix of HUGO gene symbols (rows) by cell identifiers (columns)
labFile=args[3]  #Full file path one-column matrix of unique cell identifiers (rows) by the SVM-predicted cell type (column; "label") abbreviation
outDir=args[4]  #Full path to desired output directory

dat=read.table(rawFile, header=T,row.names=1)
lab=read.table(labFile, header=T,row.names=1)

cell.types=names(table(lab$label))

out.pb=data.frame(row.names=rownames(dat))
out.labeled=data.frame(row.names=1:floor((ncol(dat)*0.1)))
out.fxn=data.frame(row.names=names(table(lab$label)))

for ( i in 1:1000){
  cn=paste(cohort,"_pbulk_",i,sep='')
  samp=dat[,sample(colnames(dat),size=ncol(dat)*0.1)]	#Select random 10% of all cells
  pb=rowSums(samp)
  out.pb=cbind(out.pb,pb)
  colnames(out.pb)[i]=cn

  labeled=paste(colnames(samp),lab[colnames(samp),"label"],sep='_')	# Annotate selected cells with SVM-predicted cell type identifier
  out.labeled=cbind(out.labeled,labeled)
  colnames(out.labeled)[i]=cn

  tbl=table(lab[colnames(samp),"label"])
  out.fxn=cbind(out.fxn,data.frame(round(tbl[names(tbl) %in% cell.types]/nrow(out.labeled),5))[,2])	# Calculate fraction of new pseudobulk population represented by each cell type present in the total population
  colnames(out.fxn)[i]=cn
}

write.table(out.pb,file.path(outDir, paste0(cohort,".pseudobulk.txt")), col.names=T, row.names=T,sep='\t', quote=F)
write.table(out.labeled,file.path(outDir, paste0(cohort,".pseudobulk.labeled.txt")), col.names=T, row.names=T,sep='\t', quote=F)
write.table(out.fxn,file.path(outDir, paste0(cohort,".pseudobulk.fxn.txt")), col.names=T, row.names=T,sep='\t', quote=F)
