#DEA1.0 - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2015)
#
#Differential Expression Analysis
#INPUT: # This script takes in input the file with the differential expressed genes for a specific phase-change
#OUTPUT: and extract two different files: 
#   - a file with genes upregulated in the change (Fold Change > 0)
#   - a file with genes downregulated (Fold Change < 0)
#

#R CMD BATCH --no-save --no-restore '--args <input_folder> <diff_exp_1> [.. <diff_exp_N>]' get_up_and_down.R
#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)
#The samples considered are 4 - please comment whenever you will not use all
length(args)
if (length(args)<3){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <input_folder> <out_folder> <diff_exp_1> [.. <diff_exp_N>]' get_up_and_down.R")
}

#The folder where the tables are
input_folder = args[1]
input_folder

#The output folder
out_folder = args[2]
out_folder

tables = args[3:length(args)]#Get the tables

#Output folder
out_folder = paste(input_folder,out_folder,sep="/")
dir.create(out_folder,showWarnings = TRUE, mode = "0777")

#read all the files in a list of data frames
ldf <- lapply(tables, function(x) read.table(file=paste(input_folder,x,sep="/"),sep="\t",head=T,quote='',comment.char='',stringsAsFactors=F,check.names=F))


get_up_down = function (data,name,i){
  #data = read.table(datafile,header=T,quote='',stringsAsFactors=F)
	#name
  #Extracts up and down-regulated
  upregulated = data[ data$logFC > 0, ]
  downregulated = data[ data$logFC < 0, ]

  #print a file with only the filtered ones
  #write.table(subset(upregulated, select=c('Row.names','logFC')),file=paste(out_folder,paste(gsub("\\.txt","",name),'upregulated.txt',sep='_'),sep="/") ,sep='\t',quote=F, row.names = FALSE, col.names = FALSE)
  #write.table(subset(downregulated, select=c('Row.names','logFC')),file=paste(out_folder,paste(gsub("\\.txt","",name),'downregulated.txt',sep='_'),sep="/") ,sep='\t',quote=F, row.names = FALSE,col.names = FALSE)
  
  write.table(upregulated,file=paste(out_folder,paste('upregulated',name,sep='_'),sep="/") ,sep='\t',quote=F, row.names = FALSE)
  write.table(downregulated,file=paste(out_folder,paste('downregulated',name,sep='_'),sep="/") ,sep='\t',quote=F, row.names = FALSE)
}

#Execute on all the tables
lapply(seq_along(ldf),function(x,name,i) get_up_down(x[[i]],name[[i]],i),x=ldf,name=tables)
