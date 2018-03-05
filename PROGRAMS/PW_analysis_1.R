#DEA1.0 - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2015)
#
#-------------------------------------------------
# PARAMETERS TO SET UP BEFORE TO RUN THE ANALYSIS
#-------------------------------------------------
#R CMD BATCH --no-save --no-restore '--args  <input_folder> <annocript_filt_out> <significant_transcripts_file> <organism_name> <p-value> <min_transcr> <min_respect_to[ann|sign]>' PW_analysis_1.R

#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)
#The samples consiedered are  - please comment whenever you will not use all
length(args)
if (length(args)<7){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <input_folder> 
	<path_to_annocript_filt_out>  <significant_transcripts_file> <organism_name> <p-value> <min_transcr> <min_respect_to[ann|sign]>' PW_analysis_1.R")
}

#Folder where all is taken and output goes
input_folder = args[1]
input_folder
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
col.signifTable = 1

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
#Returns a list for the field BPId with the elements separated by the separator ]---[ else the complete string
map.pw.l1 = function(table) {
  l1id = strsplit(table$PwLev1,"]---[",fixed=T)
}

#Returns a list for the field BPId with the elements separated by the separator ]---[ else the complete string
map.pw.l2 = function(table) {
  l2id = strsplit(table$PwLev2,"]---[",fixed=T)
}
#Returns a list for the field BPId with the elements separated by the separator ]---[ else the complete string
map.pw.l3 = function(table) {
  l3id = strsplit(table$PwLev3,"]---[",fixed=T)
}

#Takes input a table with for each transcript its GO ids and returns a list with occurrencies of each Pathway
get.counts = function(counts) {
  counts = as.data.frame(table(unlist(rapply(counts,function(x) unique(x)))))#The table function get the frequencies of the terms (amount of '-' is in)
  colnames(counts) = c('pwid','count')
  counts = counts[counts$pwid != "-",]#Removes the amount for "-"
  counts
}


#Calculate the enrichment based on the count of transcripts with a specified GO term in the output of Annocript
calculate.enrichments.ann = function(div.sel,div.uni,n.sel,n.uni) {
  div = merge(div.sel,div.uni,by.x='pwid',by.y='pwid',all.y=T)#Merge the two tables giving the same order
  div$count.x[is.na(div$count.x)] = 0#Remove NA values (they compare if something is missing)
  
  div$pval = apply(div,1,function(x) prop.test(c(as.numeric(x[2]),as.numeric(x[3])),c(n.sel,n.uni),alternative=prop.alt)$p.value)
  if(prop.alt == 'g') div = div[(div$count.x/n.sel >= mult*(div$count.y/n.uni)) & div$count.y>=min,]
  if(prop.alt == 'l') div = div[(div$count.x/n.sel <= mult*(div$count.y/n.uni)) & div$count.y>=min,]
  if(prop.alt == 't') div = div[div$count.y>=min,]
  div$padj = p.adjust(div$pval)
  div[order(div$padj,decreasing=T),]
}

#Calculate the enrichment based on the count of transcripts with a specified GO term among the significants
calculate.enrichments = function(div.sel,div.uni,n.sel,n.uni) {
  div = merge(div.sel,div.uni,by.x='pwid',by.y='pwid',all.y=T)
  div$count.x[is.na(div$count.x)] = 0
  
  div$pval = apply(div,1,function(x) prop.test(c(as.numeric(x[2]),as.numeric(x[3])),c(n.sel,n.uni),alternative=prop.alt)$p.value)
  if(prop.alt == 'g') div = div[(div$count.x/n.sel >= mult*(div$count.y/n.uni)) & div$count.x>=min,]
  if(prop.alt == 'l') div = div[(div$count.x/n.sel <= mult*(div$count.y/n.uni)) & div$count.x>=min,]
  if(prop.alt == 't') div = div[div$count.x>=min,]
  div$padj = p.adjust(div$pval)
  div[order(div$padj,decreasing=T),]
}
##########################

#Reading the tables
uni.t = read.table(file=paste(input_folder,anno,sep="/"),sep='\t',header=T,quote='',comment.char='',stringsAsFactors=F,row.names=1)
d = read.table(file=paste(input_folder,signifTable,sep="/"),sep='\t',header=T,quote='',comment.char='',stringsAsFactors=F)


#Filtering using the Annocript Output
sel = unique(d[,col.signifTable])
sel[1]
sel.t = uni.t[rownames(uni.t) %in% sel,]#Creating a filtered annocript output with only the elements present in sel


n.uni = length(unique(rownames(uni.t)))#number of unique ids in the annocript output
n.sel = length(sel)#number of rows in sel

#############Level 1 pathways
map.uni.l1 = map.pw.l1(uni.t)#the PwLev1 column from the 'ann' table (including '-')
uni.l1 = get.counts(map.uni.l1)# a table with pathway -> occurrences for level1 from uni
map.sel.l1 = map.pw.l1(sel.t)#the PwLev1 column from the 'sel' table
sel.l1 = get.counts(map.sel.l1)# a table with pathway -> occurrences for level1 from sel
if ( min_respect_to == "sign"){
	res.l1 = calculate.enrichments(sel.l1,uni.l1,n.sel,n.uni)
	min_respect_to
}
if ( min_respect_to == "ann"){
	res.l1 = calculate.enrichments.ann(sel.l1,uni.l1,n.sel,n.uni)
	min_respect_to
}
sig.l1 = subset(res.l1,padj<=p.filt)#Filter based on the pvalue given in input

#############Level 2 pathways
map.uni.l2 = map.pw.l2(uni.t)
uni.l2 = get.counts(map.uni.l2)
map.sel.l2 = map.pw.l2(sel.t)
sel.l2 = get.counts(map.sel.l2)
if ( min_respect_to == "sign"){
	res.l2 = calculate.enrichments(sel.l2,uni.l2,n.sel,n.uni)
	min_respect_to
}
if ( min_respect_to == "ann"){
	res.l2 = calculate.enrichments.ann(sel.l2,uni.l2,n.sel,n.uni)
	min_respect_to
}
sig.l2 = subset(res.l2,padj<=p.filt)

#############Level 3 pathways
map.uni.l3 = map.pw.l3(uni.t)
uni.l3 = get.counts(map.uni.l3)
map.sel.l3 = map.pw.l3(sel.t)
sel.l3 = get.counts(map.sel.l3)
if ( min_respect_to == "sign"){
	res.l3 = calculate.enrichments(sel.l3,uni.l3,n.sel,n.uni)
	min_respect_to
}
if ( min_respect_to == "ann"){
	res.l3 = calculate.enrichments.ann(sel.l3,uni.l3,n.sel,n.uni)
	min_respect_to
}
sig.l3 = subset(res.l3,padj<=p.filt)#

#PRINTING A TABULAR FILE WITH ALL THE RESULTS
restab = rbind(sig.l1,sig.l2,sig.l3)
write.table(restab,file=paste(input_folder,paste(gsub(' ','_',title),'PW_enriched.txt',sep='_'),sep="/") ,row.names=F,sep='\t',quote=F)


#PRINTING THE PLOTS
pdf(file= paste(input_folder,paste(gsub(' ','_',title),'PW_enriched.pdf',sep='_'),sep="/") ,width=15,height=10)
par(las=2,mar=c(5,25,5,5))

sig.l1 = tail(sig.l1,topn.toplot)
if(nrow(sig.l1) >= 1)
barplot(
  t(data.frame(selected=sig.l1$count.x/n.sel*100,universe=sig.l1$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.l1$pwid,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'Pathways Level 1 Enriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

sig.l2 = tail(sig.l2,topn.toplot)
if(nrow(sig.l2) >= 1)
barplot(
  t(data.frame(selected=sig.l2$count.x/n.sel*100,universe=sig.l2$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.l2$pwid,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'Pathways Level 2 Enriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

sig.l3 = tail(sig.l3,topn.toplot)
if(nrow(sig.l3) >= 1)
barplot(
  t(data.frame(selected=sig.l3$count.x/n.sel*100,universe=sig.l3$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.l3$pwid,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'Pathways Level 3 Enriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

dev.off()
