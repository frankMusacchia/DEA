#DEA1.0 - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2015)
#
#creates Venn diagram with all the files given in input
#NEEDS: install.packages('VennDiagram') 
#R CMD BATCH --no-save --no-restore '--args <input_folder> <chars_to_use> <table1> <table2> [ .. <tableN>]' venn_generic.R
#(if you add parameters remember to change the (tables = args[3:length(args)]) 

#INPUT: the folder where the tables are 
#		the number of chars you want to keep from the file names as diagram titles
#		the list of tables (the first row must contain the transcripts names!)
#OUTPUT: generates a PDF file in the input folder with the Venn diagram showing all the common transcripts



#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)#Should only arguments after --args be returned?

#The samples considered are 4 - please comment whenever you will not use all
length(args)
if (length(args)<3){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <input_folder> <chars_to_use> <table1> <table2> [ .. <tableN>] ' venn_generic.R")
}

#The folder where the tables are
input_folder = args[1]
input_folder

#max chars to use for the names
maxChars = as.numeric(args[2])
maxChars
#output name
outName = "venn.pdf"


#creates a temporary directory and install there the package
#It wil be erased at the next system restart
tmp.install.packages <- function(pack, dependencies=TRUE) {
  path <- tempdir()
  ## Add 'path' to .libPaths, and be sure that it is not
  ## at the first position, otherwise any other package during
  ## this session would be installed into 'path'
  firstpath <- .libPaths()[1]
  .libPaths(c(firstpath, path))
  install.packages(pack, dependencies=dependencies,repos="http://cran.us.r-project.org", lib=path)
}


##Installs temporarily missing packages
if("VennDiagram" %in% rownames(installed.packages()) == FALSE) {
  tmp.install.packages("VennDiagram")
}


library(VennDiagram)

tables = args[3:length(args)]#Get the tables

#read all the files in a list of data frames
ldf <- lapply(tables, function(x) read.table(file=paste(input_folder,x,sep="/"),sep="\t",head=T,quote='',comment.char='',stringsAsFactors=F))

#extract the column with names from each element of the list
ldf.rowNames = sapply(ldf, function(x) x$Row.names)

#get the titles for each ball of the Venn
outFileNames = sapply(tables, function(x) substr(x,0,maxChars))

#create a venn diagram
venn<-venn.diagram(ldf.rowNames,filename=paste(input_folder,outName,sep="/"),category=outFileNames,
	cat.pos = 0,#keeps the names up on the balls
	height = 3000, width = 3000, resolution = 500, 
	main = outName,
	fill=rainbow(length(outFileNames)))#simply get a rainbow of colors
