#DEA1.0 - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2015)
# Description: This script permits to execute a test of proportion giving
#							the enrichment of KO terms associated to a list of transcripts
#							given in input (significant_transcripts_file) compared with
#							the those in the complete Annocript output (path_to_annocript_filt_out)
### BE CAREFUL!!!!
# The Kegg Annotation system does not allow a free programmatical access. Hence you need
# to start the annotation using the web tool and later add the annotation to the Annocript table
##################
#-------------------------------------------------
# PARAMETERS TO SET UP BEFORE TO RUN THE ANALYSIS
#-------------------------------------------------
# input_folder: the folder where the files are and where the output should go
# annocript_filt_out: Annocript output
# significant_transcripts_file: table containing a column with a list of uniq transcripts
#	organism_name: name of the organism you are analyzing. It will be used for output name
# p-value: p-value to define the significance of the enrichments to show
#	min_transcr: minimum number of transcripts associated to a certain KO to make such term significant
#	min_respect_to: if the min_transcr should be referred to those in the Annocript output of in the list of significant

# Example run:
#
#R CMD BATCH --no-save --no-restore '--args  <out_folder> <annocript_filt_out> 
#<significant_transcripts_file> <organism_name> <p-value> <min_transcr> <min_respect_to[ann|sign]>' PW_analysis_1.R

#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)
#The samples consiedered are  - please comment whenever you will not use all
length(args)
if (length(args)<7){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <out_folder>
	<path_to_annocript_filt_out>  <significant_transcripts_file> <organism_name>
	 <p-value> <min_transcr> <min_respect_to[ann|sign]> ' PW_analysis_1.R")
}

#Folder where all is taken and output goes
out_folder = args[1]
out_folder

# Indicate name and path of the ANNOCRIPT TABLE containing the annotations of the whole transcriptome by ANNOCRIPT
anno =  args[2]#'cymodocea_filt_cpm1per2_uniref_2014_08_ann_out.txt'
anno

# Indicate name and path of the file containing a column with the id of a set of transcripts to test for GO enrichments
# They may be differentially expressed transcripts or stage specific transcripts or whatever
# The important thing is that they must derive from the same id present into the ANNOCRIPT TABLE
signifTable = args[3]#= 'NT_TT_significant.txt'
signifTable

#Organism name
organism = args[4] #
organism

#MOST CHANGED PARAMETERS

# Adjusted p-value cutoff to consider a class as significant (def: 0.1)
# 0.1 is that 10% are false positives
p.filt = as.numeric(args[5]) # 0.1
p.filt
# Indicate the minimum number of transcripts associated to a GO class to take it into consideration for the analysis
# It is relative to the number of transcript for each class into the "signifTable" variable (def: 5)
min = as.numeric(args[6]) #5 
min

#If it it "sign", the minimum number of transcripts is from the significant table  while if it is
#"ann" it is from the annocript output table
min_respect_to = args[7]
min_respect_to


#NOT OFTEN CHANGED
#To print the table with percentages
perc_tab = 0

# Indicate the maximum number of significant classes to display into the plot
topn.toplot = 15

# Indicate the column containing the id of the transcripts into the "signifTable" file
# Put 0 if it corresponds to rownames (def:1 with Annocript output)
col.signifTable = 18

# Indicate the fold the proportion of "diff" for a class should be higher than the universe to analyze it (def: 1)
mult = 1

# Indicate what kind of enrichment you are looking for in the "signifTable" table
# To look for enrichments use 'g', for impoverishment use 'l', for both use 't'
prop.alt = 'g'

#Parameters used
params = paste(organism,paste("min",min,sep=""),paste("pval",gsub("\\.","",toString(p.filt)),sep=""),sep="_")#'cymodocea_min5_pval01 NT_TT_significant'

#Experiment name
experName = strsplit(signifTable,"\\.")

#The title to assign to various files in ouput
title = paste(params,unlist(experName)[1]," ")

#-------------------------------------------------------
# END OF PARAMETERS TO SET UP BEFORE TO RUN THE ANALYSIS
#-------------------------------------------------------


#-------------------------------------------------------
# ANALYSIS
#-------------------------------------------------------

#########################
#Returns a list for the field ko_descs with the elements separated by the separator ]---[ else the complete string
map.pw.ko = function(table) {
  koid = strsplit(table$ko_descs,"]---[",fixed=T)
}

#Takes input a table with for each transcript its GO ids and returns a list with occurrences of each GO term
get.counts = function(counts) {
  counts = as.data.frame(table(unlist(counts)))#The table function get the frequencies of the terms (amount of '-' is in)
  colnames(counts) = c('ko_descs','count')
  counts = counts[counts$ko_descs != "-",]#Removes the amount for "-"
  counts
}


#Calculate the enrichment based on the count of transcripts with a specified GO term among the significants
calculate.enrichments = function(sel.ko,uni.ko,n.sel,n.uni,min_respect_to) {
  div = merge(sel.ko,uni.ko,by.x='ko_descs',by.y='ko_descs',all.y=T)
  div$count.x[is.na(div$count.x)] = 0
  
  div$pval = apply(div,1,function(x) prop.test(c(as.numeric(x[2]),as.numeric(x[3])),c(n.sel,n.uni),alternative=prop.alt)$p.value)
  div$padj = p.adjust(div$pval)
    
  #Add columns with totals n.sel and n.uni
  div$tot.sel <- n.sel
  div$tot.uni <- n.uni
  #Add columns with percentages 
  div$perc.sel <- (div$count.x / n.sel )*100
  div$perc.uni <- (div$count.y / n.uni )*100

	if (min_respect_to == "ann"){
		if(prop.alt == 'g') div = div[(div$count.x/n.sel >= mult*(div$count.y/n.uni)) & div$count.y>=min,]
		if(prop.alt == 'l') div = div[(div$count.x/n.sel <= mult*(div$count.y/n.uni)) & div$count.y>=min,]
		if(prop.alt == 't') div = div[div$count.y>=min,]
	}
	if (min_respect_to == "sel"){
		if(prop.alt == 'g') div = div[(div$count.x/n.sel >= mult*(div$count.y/n.uni)) & div$count.x>=min,]
		if(prop.alt == 'l') div = div[(div$count.x/n.sel <= mult*(div$count.y/n.uni)) & div$count.x>=min,]
		if(prop.alt == 't') div = div[div$count.x>=min,]
	}
  #Re-order the columns
  div <- div[,c(1,2,3,8,9,6,7,4,5)]
  #Change column names
  colnames(div)[c(2,3)] <- c("count.sel","count.uni") 
	
	if ( min_respect_to == "get_all"){
	  #Sort in increasing order
		div.sorted <- div[order(div$padj,decreasing=F),]
		#PRINTING A TABULAR FILE WITH ALL THE RESULTS
		write.table(div.sorted,file=paste(out_folder,paste(gsub(' ','_',title),'ko_PW_table.txt',sep='_'),sep="/") ,row.names=F,sep='\t',quote=F)
	}else{
		#Sort in increasing order
		div <- div[order(div$padj,decreasing=T),]
	}
		#return table
		div
}


	#Reading the tables
	uni.t = read.table(file=anno,sep='\t',header=T,quote='',comment.char='',stringsAsFactors=F,row.names=1)
	d = read.table(file=signifTable,sep='\t',header=T,quote='',comment.char='',stringsAsFactors=F)


	#Filtering using the Annocript Output
	sel = unique(d[,col.signifTable])
	sel[1]
	sel.t = uni.t[rownames(uni.t) %in% sel,]#Creating a filtered annocript output with only the elements present in sel


	n.uni = length(unique(rownames(uni.t)))#number of unique ids in the annocript output
	n.sel = length(sel)#number of rows in sel

	#############ko pathways
	map.uni.ko = map.pw.ko(uni.t)#the ko column from the unierse (including '-')
	uni.ko = get.counts(map.uni.ko)# a table with pathway -> occurrences for ko from uni
	map.sel.ko = map.pw.ko(sel.t)#the ko column from the 'sel' table
	sel.ko = get.counts(map.sel.ko)# a table with pathway -> occurrences for ko from sel

	res.ko = calculate.enrichments(sel.ko,uni.ko,n.sel,n.uni,min_respect_to)

	sig.ko = subset(res.ko,padj<=p.filt)#Filter based on the pvalue given in input

	#Plots will be print only if the enrichment were performed (min_respect_to = ann or sel)
	if ( min_respect_to != "get_all"){
		#PRINTING A TABULAR FILE WITH ALL THE RESULTS
		write.table(sig.ko,file=paste(out_folder,paste(gsub(' ','_',title),'ko_PW_enriched.txt',sep='_'),sep="/") ,row.names=F,sep='\t',quote=F)
		#PRINTING THE PLOTS
		pdf(file= paste(out_folder,paste(gsub(' ','_',title),'ko_enriched.pdf',sep='_'),sep="/") ,width=15,height=10)
		par(las=2,mar=c(5,25,5,5))

		sig.ko = tail(sig.ko,topn.toplot)
		if(nrow(sig.ko) >= 1)
		barplot(
			t(data.frame(selected=sig.ko$count.sel/n.sel*100,universe=sig.ko$count.uni/n.uni*100)),
			beside=T,horiz=T,names=unlist(lapply(sig.ko$ko_descs,substring,1,50)),xlab='percentage of transcripts',
			main=paste('Top',topn.toplot,'Kegg Pathways Enriched Among',title,sep=' '),
			cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Universe'))
		dev.off()
	}

