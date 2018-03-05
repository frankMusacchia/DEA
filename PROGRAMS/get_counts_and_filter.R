#DEA1.0 - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2015)
#
#INSTRUCTIONS FOR THE EXECUTION. PLEASE READ CAREFULLY TO USE PROPERLY
#EXAMPLE COMMAND:
#R CMD BATCH --no-save --no-restore '--args  <input_folder> <target_file> <organism_name> <cpm_min> <min_samples>' get_counts_and_filter.R

#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)

################HOW TO USE
##This R script has to basic purposes: 
#########The first purpose is (extract_counts_file). It creates a single
# file with counts from a set of file with counts as they have been created
# from samtools.

#The output file will have as much column as the number of files present
# in the folder
input_folder = args[1]

#if your file does not come from samtools but other software please keep 
#it tab separated and change this variable properly
column_with_count = args[2]

#The names of the file will be searched by taking the samples names from the file
target_file = args[3]#"sample_file" #it has 2columns: the first with the sample name and the second correspond to the experiment name

#The final file with counts will have the name of the organism you are using and giving in input
organism = args[4] 
counts_file = paste(paste(organism,"counts",sep="_"),"txt",sep=".") #"organism_counts.txt"


############The second function of this script (extract_filtered_percpm_counts_file), given the complete file with counts
#that it created before (counts_file) it will extract a list of names of transcripts by select only those that have:
# at least milfilt as cpm in at least minlibfilt random replicates
milfilt = as.numeric(args[5]) #1
minlibfilt = as.numeric(args[6]) #2

#it will give in output the file with the cpm filtered
cpm_filt = paste(paste(organism,"cpm","filtered",sep="_"),"txt",sep=".") #'organism_cpm_filtered.txt'


#the file with the counts filtered
filt_counts_file = paste(paste(organism,"counts","filtered",sep="_"),"txt",sep=".") #'organism_counts_filtered.txt'

#and finally the list of the names of the transcripts filtered
filt_tr_list_file = paste(paste(organism,"transcripts","filtered","list",sep="_"),"txt",sep=".") #"organism_filtered_transcripts_list.txt"




#######FUNCTIONS


#After the quantization done with an alignment program you can start this simple script
#to get the file with all the counts in a unique table
extract_counts_file = function (){

  samples = read.table(file=paste(input_folder,target_file,sep="/"),sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)

  #Get the filenames from the list depending from the fact that the extension is .count
  #counts_file_names = list.files(folder_path, pattern="*.counts", full.names=TRUE)

  #get the file names reading the sample names in the given file and concatenating ".counts"
  sorted_file_names = sapply(samples$name, function(x) paste(x,"counts",sep=".") )

  #read all the files in a list of data frames
  ldf <- lapply(sorted_file_names, function(x) read.table(file=paste(input_folder,x,sep="/"),sep="\t"))

  #extract the column with counts from each element of the list
  counts = sapply(ldf, function(x) x[column_with_count])

  #Bind the table of counts with a column of transcripts names
  counts_table = cbind(ldf[[1]][1],counts)

  colnames(counts_table) <- c('id',samples$name)

  #Write on a tab separated file
  write.table(counts_table,file=paste(input_folder,counts_file,sep="/"),sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
}


MyMerge <- function(x, y,column_with_count){
#	x[-c(2:(column_with_count-1),(column_with_count+1):ncol(x))]
#	y[-c(2:(column_with_count-1),(column_with_count+1):ncol(y))]

#x[,!c(1,column_with_count)]<- list(NULL)
#y[,!c(1,column_with_count)]<- list(NULL)
  x$V2<-NULL
  x$V4<-NULL
  y$V2<-NULL
  y$V4<-NULL    
  df <- merge(x, y, by="V1", all.x= F, all.y= F,sort=FALSE)
  return(df)
}


#After the quantization done with an alignment program you can start this simple script
#to get the file with all the counts in a unique table
merge_counts_file = function (){

  samples = read.table(file=paste(input_folder,target_file,sep="/"),sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)

  #get the file names reading the sample names in the given file and concatenating ".counts"
  sorted_file_names = sapply(samples$name, function(x) paste(x,"counts",sep=".") )

  #read all the files in a list of tables
  ldf <- lapply(sorted_file_names, function(x) read.table(file=paste(input_folder,x,sep="/"),sep="\t"))
  
  counts_table = Reduce(function(x,y) MyMerge(x,y,column_with_count),ldf)
  
  colnames(counts_table) <- c('id',samples$name)

  #Write on a tab separated file
  write.table(counts_table,file=paste(input_folder,counts_file,sep="/"),sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
}

 

#If you want to exctract a file of transcripts that has at least minlibfilt replicates
#with at least milfilt as cpm value, give to this function this two values and
#give it the counts_file and the target_file
extract_filtered_percpm_counts_file = function (counts_file){
  
  #Install edgeR
  source("http://bioconductor.org/biocLite.R")
  biocLite("edgeR")

  #Use the edgeR library
  library(edgeR)
  
  #get the table with the target samples and conditions
  samples = read.table(file=paste(input_folder,target_file,sep="/"),sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
  
  #Assing the rownames of the samples
  rownames(samples) = samples$name
  
  #Read the file of the counts for each replicate
  data = read.table(file=paste(input_folder,counts_file,sep="/"),head=T,quote='',check.names=F)
  
  #the rownames will be the ids of the transcripts
  rownames(data) = data$id
  
  #Order the samples data based on the names in the target files
  data = data[,samples$name]
  
  #Create an array with the for each replicate the sum of the counts of all the mapping reads
  tot = colSums(data)
  
  #
  lib.size=tot[samples$name]
  
  #Create the DGElist object to give to the cpm function of EdgeR
  dge = DGEList(data, group=samples$condition,lib.size=lib.size[samples$name])

  #Start EdgeR cpm counting
  cpm = cpm(dge,normalized.lib.sizes=F,log=FALSE,prior.count=0)
  
  #Select only those transcripts that has at least milfilt as cpm in at least minlibfilt random replicates
  ridx = rowSums(cpm >= milfilt) >= minlibfilt
  
  #Get the filtered dge file 
  dge = dge[ridx,]
  
  #Get the filtered cpm file
  cpm = cpm[ridx,]
  
  #Get dimensions of the filtered cpm file
  dim(cpm)
  
  #Write the file wit the filtered counts
  write.table(dge$counts,file=paste(input_folder,filt_counts_file,sep="/"),sep='\t',quote=F)
  
  
  #Write the file with the filtered cpm
  write.table(cpm,file=paste(input_folder,cpm_filt,sep="/"),sep='\t',quote=F,col.names=NA)#,row.names=FALSE,col.names=FALSE)
  
  #Write the file with the filtered cpm FOR DELIVERY
  #write.table(cpm,file=paste(input_folder,paste("DELIVERY",cpm_filt,sep="_"),sep="/"),sep='\t',quote=F)
  
  #Write only the column with the transcripts name
  write.table(row.names(cpm),file=paste(input_folder,filt_tr_list_file,sep="/"),quote=F, ,row.names=FALSE,col.names=FALSE)
}

#EXECUTION


#extract_counts_file()#OLD FUNCTION

merge_counts_file()
if (!(minlibfilt == 0 & milfilt == 0)){
	extract_filtered_percpm_counts_file(counts_file)
}

