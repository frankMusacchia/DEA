#DEA1.0 - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2015)
#
#Differential Expression Analysis
#INPUT: the sample_file with the grouping of replicates
#       the filtered table with transcripts counts 
#						(We need filtered file of counts since we make the analysis
#						 only on the significant transcripts)
#	
#R CMD BATCH --no-save --no-restore '--args <input_folder>  <counts_filtered> <organism_name> <target_file> <minFDR> <out_suffix> <analysis_type> <annocript_out>' DE_analysis_generic.R

#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)
#The samples considered are 4 - please comment whenever you will not use all
length(args)
if (length(args) < 5){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <input_folder> <counts_filtered> <organism_name> <target_file> <minFDR> <out_suffix>  <analysis_type> <annocript_out> ' DE_analysis_generic.R")
}
###############################################METTERE IL FOLD CHANGE SOGLIA COME PARAMETRO IN INPUT

#Folder where all is taken and output goes
input_folder = args[1]
input_folder

counts_filtered = args[2]#Raw counts filtered before (not CPM)!!
counts_filtered
organism_name = args[3]
organism_name
target_file = args[4]
target_file

# Minimum FDR Default: 0.1
minFDR = as.numeric(args[5])
minFDR

#Suffix of the output name (use "NULL" if not present)
out_suffix = args[6]

# Type of analysis: 
	#1 do differential expression.
	#2 only get the exactTest and extract table
	#3 get a stacked plot of expression 
analType = 	as.numeric(args[7])

#If the Annocript output is used
ann_out_used = 7
if (length(args) > ann_out_used){
	ann = args[8]#"chaetoceros_transcriptome_cpm1per2_uniref_2014_08_filt_ann_out.txt"
}


################################DOWNLOAD PACKAGES FUNCTION
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


########################

get_diff_exp = function (couple){

  exp_name = paste(paste(couple,collapse="_"),fdrString,out_suffix,sep="_")
  exp_name
  # Espressione differenziale tra le due fasi
  # Execute the statistical test (fisher-test) to see differences between different samples
  diff.ab = exactTest(dge,pair=couple)
  
  # Extract the most relevant results
  tt.ab = topTags(diff.ab,n=nrow(diff.ab))
  
  # Put the results into a table
  res.ab = tt.ab$table
  
  # Select the significant ones according to the cutoffs
  sig.ab = res.ab[res.ab$FDR<=minFDR & abs(res.ab$logFC )>=log2(2),]
  
  # Merge the statistical results with the normalized counts without sorting 
  sig.ab = merge(sig.ab, cpm, by=0, all.x=T, sort=F)
  
  # Write a table with all the info for the significantly differentially expressed genes
  write.table(sig.ab, file= paste(input_folder,exp_name,sep="/") , sep='\t', quote=F, row.names=F)
	
	if (length(args) > ann_out_used){
		join_with_annotation(exp_name)
	}
}

#Gets expression table of a given couple of samples by using the exact test
#with the file of counts
get_expression_table = function (couple){

  exp_name = paste(paste(couple,collapse="_"),out_suffix,sep="_")
  exp_name
  
  # Espressione differenziale tra le due fasi
  # Execute the statistical test (fisher-test) to see differences between different samples
  ####Annamaria suggests to use glm (generalized linear model).
  diff.ab = exactTest(dge,pair=couple)
  
  # Extract the most relevant results
  tt.ab = topTags(diff.ab,n=nrow(diff.ab))
  
  # Put the results into a table
  res.ab = tt.ab$table
  
  # Write a table with all the info for the significantly differentially expressed genes
  write.table(res.ab, file= paste(input_folder,exp_name,sep="/") , sep='\t', quote=F, row.names=T)
	
}


#Merge the significants table with the annocript annotation and print a new file
join_with_annotation = function (exp_name){

  cols = list("Row.names","logFC","FDR","HSPNameSP","HSPEvalueSP","DescriptionSP","EnzymeDescs","PwLev1","PwLev2","PwLev3","DescriptionUf","HSPEvalueUf","Taxonomy","BPDesc","MFDesc","CCDesc","CDDesc")
	#cols = list("Row.names","logFC","FDR","HSPNameSP","HSPScoreSP","DescriptionSP","EnzymeDescs","Pathways","DescriptionUf","HSPScoreUf","Taxonomy","BPDesc","MFDesc","CCDesc","CDDesc")#OLDANNOTATION
	
	tableA = read.table(file=paste(input_folder,exp_name,sep="/"),sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F,row.names=1)
	ann_out = read.delim(file=paste(input_folder,ann,sep="/"),sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F,row.names=1)	
	
	table_and_ann_out = merge(tableA,ann_out,by="row.names",all=F, sort=F)
	
	final_table = table_and_ann_out[,unlist(cols)]#get only the columns from cols
	
	colnames(final_table) = c("TranscriptNames",cols[2:length(cols)])

	write.table(final_table,file=paste(input_folder,paste(paste(unlist(strsplit(exp_name,"\\."))[1],"and_ann_out",sep="_"),"txt",sep="."),sep="/") ,sep='\t',quote=F,row.names=F)

}

#Print the hierarchical clustering of the counts
print_hclust = function(){
	# Clustering of the experiments
	exprs = dge$pseudo.counts

	ecTr = dist(t(exprs), method = "euclidean")
	hecTr = hclust(ecTr, method = "average")

	pdf(paste(input_folder,"hier_clustering_average.pdf",sep="/"),width=7,height=5)
	plot(hecTr, main = "Hierarchical clustering dendrogram for counts", xlab = "", sub = "Average linkage, Euclidean distance for counts")
	dev.off()

	hecTr = hclust(ecTr, method = "complete")
	pdf(paste(input_folder,"hier_clustering_complete.pdf",sep="/"),width=7,height=5)
	plot(hecTr, main = "Hierarchical clustering dendrogram for counts", xlab = "", sub = "complete linkage, Euclidean distance for counts")
	dev.off()
	#plot(hecTr, main = "Hierarchical clustering dendrogram for cpm", xlab = "", sub = "Average linkage, Euclidean distance for cpm")
}

stacked_plot = function (in_table) {

#NUOVO PLOT
#GGplot consists in a function that creates a plot, to which you can 
	#add different features
	###Installs temporarily missing packages
	if("reshape2" %in% rownames(installed.packages()) == FALSE) {
		#install.packages("RColorBrewer", repos="http://cran.us.r-project.org")
	  tmp.install.packages("reshape2")
	}
	#else{print("RColorBrewer already installed")}
 
 data <- read.table(file=paste(input_folder,in_table,sep="/"), header = F, stringsAsFactors = F,row.names=1,sep="\t")
 data <- data[nrow(data):1,] #Invert the order of rows
 colnames(data) <- c("up","nde","down")
 
 #The reshape library permits to reshape a table in a confortable way
 #to use whenever a function needs
 library(reshape2)#Call Reshape2 library
 data2 <- melt(t(data))
 
 library(ggplot2)
 pdf(file=paste(input_folder,paste(in_table,"pdf",sep="."),sep='/') ,width=15,height=10)
 #pdf_file = paste(input_folder,paste(in_table,"pdf",sep="."),sep="/")
 #pdf_file
 #jpeg(file=pdf_file, width = 680, height = 480, units = "px", pointsize = 12,
   #  bg = "white")
 colors_data <-c("red","grey","green")
 exp_plot <- ggplot(data=data2, aes(x=Var2, y=value, fill=Var1)) + 
	geom_bar(stat="identity") +
	xlab("\nGO classes") +
	ylab("Percentage\n") +
	coord_flip() +
	theme_bw() +
	guides(fill = guide_legend(title = "Legend")) +
	scale_fill_manual(values=colors_data)
 print(exp_plot)
 dev.off()
 #ggsave(filename=pdf_file,plot=exp_plot,)
}



###########MAIN

#Only if an analysis of expression is to be performed we run this part
if (analType == 1 || analType == 2){
	#FDR string to put in the file name
	fdrString = paste("minFDR",gsub("\\.","",toString(minFDR)),sep="")

	# Load edgeR library
	library("edgeR")


	# Load the sample file (name and condition has to be equal for the replicates)
	samples = read.table(file=paste(input_folder,target_file,sep="/") , sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F)

	# Load the raw counts file filtered (we use the raw since the normalization is done by edgeR function later)
	counts = read.table(file = paste(input_folder,counts_filtered,sep="/") , row.names=1, head=T, stringsAsFactors=F,check.names=F)

	# Put the counts file column in the same order as the sample file
	counts = counts[,samples$name]

	# Create the edgeR object
	dge = DGEList(counts, group=samples$condition)

	#Start EdgeR cpm counting
	cpm = cpm(dge,normalized.lib.sizes=F,log=FALSE,prior.count=0)
		
	# Normalization steps: they play with variance. In example the variance among the replicates and genes that
	# have similar expression value in the same sample (?)
	dge = calcNormFactors(dge)
	dge = estimateCommonDisp(dge)
	dge = estimateTagwiseDisp(dge)#

	# Clustering of the experiments
	#exprs = dge$pseudo.counts
	
	#Execute the differential expression analysis for each combination of replicate
	conditions = unique(samples$condition)
}



#1. Differential expression analysis of all the samples against all the samples
if (analType == 1){
	for (couple in combn(conditions,2,simplify=FALSE)){
		get_diff_exp(couple)
	}
	#Print hierarchical clustering plots of the counts
	print_hclust()
}

#2. Expression analysis only
if (analType == 2){
	for (couple in combn(conditions,2,simplify=FALSE)){
		get_expression_table(couple)
	}
}

#3. Expression plots of all the couple of samples
if (analType == 3){
		stacked_plot(counts_filtered)
}



