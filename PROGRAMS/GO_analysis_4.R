#DEA1.0 - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2015)
#
#-------------------------------------------------
# PARAMETERS TO SET UP BEFORE TO RUN THE ANALYSIS
#-------------------------------------------------
#R CMD BATCH --no-save --no-restore '--args <input_folder> <annocript_filt_out> <significant_transcripts_file> <organism_name> <p-value> <min_transcr> <min_respect_to[ann|sign]> <go mapping>' GO_analysis_4.R

#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)
#The samples consiedered are 4 - please comment whenever you will not use all
length(args)
if (length(args)<8){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <input_folder> 
	<path_to_annocript_filt_out>  <significant_transcripts_file> <organism_name> <p-value> <min_transcr> <min_respect_to[ann|sign]> <go mapping>' GO_analysis_4.R")
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
# 0.1 means 10% are false positives
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
# Indicate name and path of the file containing the mapping between the GO id and the definition
gomap = args[8]

#To print the table with percentages
perc_tab = 0

# Indicate the maximum number of significant classes to display into the plot
topn.toplot = 15

# Indicate the column containing the id of the transcripts into the "signifTable" file
# Put 0 if it corresponds to rownames (def:1 with Annocript output)
col.signifTable = 1

# Indicates the fold a class should be higher than the universe to analyze it (def: 1)
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

#Returns a list for the field BPId with the elements separated by the separator ]---[ else the complete string
map.go.bp = function(table) {
  bpid = strsplit(table$BPId,']---[',fixed=T)
}

#Returns a list for the field MFId with the elements separated by the separator ]---[ else the complete string
map.go.mf = function(table) {
  mfid = strsplit(table$MFId,']---[',fixed=T)
}

#Returns a list for the field CCId with the elements separated by the separator ]---[ else the complete string
map.go.cc = function(table) {
  ccid = strsplit(table$CCId,']---[',fixed=T)
}

#Takes input a table with for each transcript its GO ids and returns a list with occurrences of each GO term
get.counts = function(counts) {
  counts = as.data.frame(table(unlist(counts)))#The table function get the frequencies of the terms
  colnames(counts) = c('goid','count')
  counts = counts[grep('GO:',counts$goid),]
  counts
}

#Calculate the enrichment based on the count of transcripts with a specified GO term among the significants
calculate.enrichments = function(div.sel,div.uni,n.sel,n.uni,go) {
  div = merge(div.sel,div.uni,by.x='goid',by.y='goid',all.y=T)#Merge the two tables giving the same order
  div$count.x[is.na(div$count.x)] = 0#Remove NA values (they compare if something is missing)
  
  div$pval = apply(div,1,function(x) prop.test(c(as.numeric(x[2]),as.numeric(x[3])),c(n.sel,n.uni),alternative=prop.alt)$p.value)
  div = merge(div,go,by.x='goid',by.y='go_id')
  if(prop.alt == 'g') div = div[(div$count.x/n.sel >= mult*(div$count.y/n.uni)) & div$count.x>=min,]
  if(prop.alt == 'l') div = div[(div$count.x/n.sel <= mult*(div$count.y/n.uni)) & div$count.x>=min,]
  if(prop.alt == 't') div = div[div$count.x>=min,]
  div$padj = p.adjust(div$pval)
  div[order(div$padj,decreasing=T),]
}

#Calculate the enrichment based on the count of transcripts with a specified GO term in the output of Annocript
calculate.enrichments.ann = function(div.sel,div.uni,n.sel,n.uni,go) {
  div = merge(div.sel,div.uni,by.x='goid',by.y='goid',all.y=T)#Merge the two tables giving the same order
  div$count.x[is.na(div$count.x)] = 0#Remove NA values (they compare if something is missing)
  
  div$pval = apply(div,1,function(x) prop.test(c(as.numeric(x[2]),as.numeric(x[3])),c(n.sel,n.uni),alternative=prop.alt)$p.value)
  div = merge(div,go,by.x='goid',by.y='go_id')
  if(prop.alt == 'g') div = div[(div$count.x/n.sel >= mult*(div$count.y/n.uni)) & div$count.y>=min,]
  if(prop.alt == 'l') div = div[(div$count.x/n.sel <= mult*(div$count.y/n.uni)) & div$count.y>=min,]
  if(prop.alt == 't') div = div[div$count.y>=min,]
  div$padj = p.adjust(div$pval)
  div[order(div$padj,decreasing=T),]
}


#Reading the tables
go = read.table(file=gomap,sep='\t',header=T,quote='',comment.char='',stringsAsFactors=F)
uni.t = read.delim(file=paste(input_folder,anno,sep="/"),sep='\t',header=T,quote='',comment.char='',stringsAsFactors=F,row.names=1)
d = read.table(file=paste(input_folder,signifTable,sep="/"),sep='\t',header=T,quote='',comment.char='',stringsAsFactors=F)


#Filtering using the Annocript Output
sel = unique(d[,col.signifTable])
sel[1]#prints the first element #DEBUGCODE
sel.t = uni.t[rownames(uni.t) %in% sel,]#Creating a filtered annocript output with only the elements present in sel


n.uni = length(unique(rownames(uni.t)))#number of unique ids in the annocript output
n.sel = length(sel)#number of rows in sel

# Biological Process
map.uni.bp = map.go.bp(uni.t)#the GOBP column from the 'ann' table (including '-')
uni.bp = get.counts(map.uni.bp)# a table with GO:id -> occurrences for BP from uni
map.sel.bp = map.go.bp(sel.t)#the GOBP column from the 'sel' table
sel.bp = get.counts(map.sel.bp)# a table with GO:id -> occurrences for BP from sel
if ( min_respect_to == "sign"){
	res.bp = calculate.enrichments(sel.bp,uni.bp,n.sel,n.uni,go)
	min_respect_to
}
if ( min_respect_to == "ann"){
	res.bp = calculate.enrichments.ann(sel.bp,uni.bp,n.sel,n.uni,go)
	min_respect_to
}
sig.bp = subset(res.bp,padj<=p.filt)

#Creating a table with supplementary information
if(nrow(sig.bp) >= 1){
  bp.perc.tab = sig.bp
  bp.perc.tab$species = signifTable
  bp.perc.tab$perc.t = (bp.perc.tab$count.y/length(map.uni.bp))*100
  bp.perc.tab$perc.s = (bp.perc.tab$count.x/length(map.sel.bp))*100
  bp.perc.tab$pval=NULL
  #bp.perc.tab$padj=NULL
  bp.perc.tab$count.y=NULL
  bp.perc.tab$count.x=NULL
}

# Molecular Function
map.uni.mf = map.go.mf(uni.t)
uni.mf = get.counts(map.uni.mf)
map.sel.mf = map.go.mf(sel.t)
sel.mf = get.counts(map.sel.mf)
if ( min_respect_to == "sign"){
	res.mf = calculate.enrichments(sel.mf,uni.mf,n.sel,n.uni,go)
	min_respect_to
}
if ( min_respect_to == "ann"){
	res.mf = calculate.enrichments.ann(sel.mf,uni.mf,n.sel,n.uni,go)
	min_respect_to
}
sig.mf = subset(res.mf,padj<=p.filt)

#Creating a table with supplementary information
if(nrow(sig.mf) >= 1){
  mf.perc.tab = sig.mf
  mf.perc.tab$species = signifTable
  mf.perc.tab$perc.t = (mf.perc.tab$count.y/length(map.uni.mf))*100
  mf.perc.tab$perc.s = (mf.perc.tab$count.x/length(map.sel.mf))*100
  mf.perc.tab$pval=NULL
  #mf.perc.tab$padj=NULL
  mf.perc.tab$count.y=NULL
  mf.perc.tab$count.x=NULL
}

# Cellular Component
map.uni.cc = map.go.cc(uni.t)
uni.cc = get.counts(map.uni.cc)
map.sel.cc = map.go.cc(sel.t)
sel.cc = get.counts(map.sel.cc)
if ( min_respect_to == "sign"){
	res.cc = calculate.enrichments(sel.cc,uni.cc,n.sel,n.uni,go)
	min_respect_to
}
if ( min_respect_to == "ann"){
	res.cc = calculate.enrichments.ann(sel.cc,uni.cc,n.sel,n.uni,go)
	min_respect_to	
}
sig.cc = subset(res.cc,padj<=p.filt)

#Creating a table with supplementary information
if(nrow(sig.cc) >= 1){
  cc.perc.tab = sig.cc
  cc.perc.tab$species = signifTable
  cc.perc.tab$perc.t = (cc.perc.tab$count.y/length(map.uni.cc))*100
  cc.perc.tab$perc.s = (cc.perc.tab$count.x/length(map.sel.cc))*100
  cc.perc.tab$pval=NULL
  #cc.perc.tab$padj=NULL
  cc.perc.tab$count.y=NULL
  cc.perc.tab$count.x=NULL
}

#Finalizing the creation of a table with all the enrichments obtained for all bp, mf, cc
restab = rbind(sig.bp,sig.mf,sig.cc)
write.table(restab,file=paste(input_folder,paste(gsub(' ','_',title),'GO_enriched.txt',sep='_'),sep="/") ,row.names=F,sep='\t',quote=F)

perctab = NULL

if (nrow(sig.bp) >= 1)  
  perctab = bp.perc.tab
if (nrow(sig.mf) >= 1)
  perctab = rbind(perctab,mf.perc.tab)   
if (nrow(sig.cc) >= 1)
  perctab = rbind(perctab,cc.perc.tab)

if (exists("perctab") & perc_tab == 1)  
  write.table(perctab,file=paste(input_folder,paste(gsub(' ','_',title),'perctab.txt',sep='_'),sep="/"),row.names=F,sep='\t',quote=F)



pdf(file=paste(input_folder,paste(gsub(' ','_',title),'GO_enriched.pdf',sep='_'),sep="/") ,width=15,height=10)
par(las=2,mar=c(5,25,5,5))

sig.bp = tail(sig.bp,topn.toplot)
if(nrow(sig.bp) >= 1)
barplot(
  t(data.frame(selected=sig.bp$count.x/n.sel*100,universe=sig.bp$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.bp$definition,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'GO Biological Process Enriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

sig.mf = tail(sig.mf,topn.toplot)
if(nrow(sig.mf) >= 1)
barplot(
  t(data.frame(selected=sig.mf$count.x/n.sel*100,universe=sig.mf$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.mf$definition,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'GO Molecular Function Enriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

sig.cc = tail(sig.cc,topn.toplot)
if(nrow(sig.cc) >= 1)
barplot(
  t(data.frame(selected=sig.cc$count.x/n.sel*100,universe=sig.cc$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.cc$definition,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'GO Cellular Component Enriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

dev.off()
