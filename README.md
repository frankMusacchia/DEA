# DEA1.0 - Differential  Expression analysis using Annocript
Author: Francesco Musacchia (2018)

**Mission**
This software permits to execute differential expression analysis using many samples from RNA-seq experiment.
It needs in input the files as they come from the count of reads which after mapping with transcripts. Forexample you could use Bowtie, STAR, etc to align your raw reads against a transcriptome and then SAMTOOLS to execute the count of reads.

If you have N samples, then you should get N files with counts. Output file from samtools is one with a column with the transcript name and count of reads aligning on it (3rd column). You can change the index of this column in CONFIGURATION/config_program.txt.

Starting from the counts files, DEA0.1 permits to:
- (A) merge the counts columns in a unique file
- (B) filter your assembled transcriptome on the base of the expression
- (C) perform the differential expression (DE) analysis of the samples using the file of counts
- (D) extract up and down regulated transcripts from each DE analysis
- (E) get GO and pathways enrichments TABLE and PLOTS of the DE transcripts using the output from Annocript
- (F) extract a Venn diagram showing intersections among the DE analyses
- (G) get the expression plots of differentially expressed GO terms


 To use it you need:

1. A target file with an header which is a list of samples with an experimental condition code associated with. Please use the same field names as in the following example:
Example:
		  name	condition
		SampleA1	A
		SampleA2	A
		SampleB1	B
		SampleB2	B

2. A folder where you put all the counts files. The counts files should 
	have be named: sample_name.counts. sample_name should be the same as in the target file. 
 i.e.: SampleA1.counts,SampleA2.counts,SampleB1.counts,SampleB2.counts

3. A configuration file with all the parameters to use for the script (see the one included in DEA as an example).

**OPTIONALLY**

4. The filtered output from Annocript. Its usage is mandatory for E and G. Optionally you can use it
when running C, in this case a file with the differentially expressed transcripts and few columns of
annotation will be output.



# Start with DEA


**Working folder configuration**

To use DEA please create a folder to work with and run the install.pl script in the program folder.

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~$ mkdir dea_works 
              frank@compaq2:~$ cd DEA0.1 
              frank@compaq2:~$ perl install.pl 
              ---------------- TERMINAL ------------------

This will ask to you the full path to the working folder you just created: 

              ---------------- TERMINAL ------------------ 
              This script will prepare DEA to work in your directory. 
                
              Are you sure you want to enjoy this software?(y or n) 
                y 
              Write an existing complete path where you want to play with DEA: /home/francesco/dea_works/ 
                
              Done! Now you can start DEA from this folder! 
              ---------------- TERMINAL ------------------ 

Notice: The install.pl script simply creates a file (folders.txt) where are written the folders that DEA will use. 
If the folders.txt file is not there DEA cannot work. 


**Run DEA**

To run DEA you can directly call it with PERL console command. Three parameters are mandatory: --input_folder, --config_file and --target_file (read the first part of this README to understand their
usage).
To start more consecutive or single analyses use the user_config file. A copy of it (from CONFIGURATION folder) is created by the installation program in the working folder.

*Get global counts and filter transcriptome*

To start with DEA you need the full fasta transcriptome assembled and the files with counts named as specified before. 
If you want to get a filtered transcriptome to annotate you must use the first three parameters and
the parameter --full_transcriptome= the path to the transcriptome created with the assembling program. 

		SELECT : filter_transcriptome = YES

The transcriptome is filtered by first computing the Count Per Million (CPM) value a
nd then fetching only those transcripts which have CPM greater than a given threshold (min_cpm) for at least a given number of replicates (min_num_repl). Please change the default parameters in the config_user.txt file.

This step expect you to run Annocript after getting the filtered transcriptome, thus it stops. But you can run DEA also without the Annocript annotation.

Parameters:
The followings are to be reach from a transcript to be defined as expressed
min_cpm: minimum count per million value
min_num_repl: minimum number of replicates

*Differential gene expression analysis with EdgeR*

DEA uses EdgeR to perform the differential gene expression analysis. Thus it applies the statistical t-test to get differences in expression. Genes are considered differentially expressed whenever the False Discovery rate and the Fold change cross or are under defined thresholds.


If you want to perform a differential gene expression analysis you need only the file of filtered counts as it comes from the *filter_transcriptome* step. You do not need to give it in input since it will be taken from the folder you choose.

		SELECT: de_analysis = YES

Parameters: 
The followings are the constrains for a transcript to be defined as differentially expressed
- minFDR, the minimum false discovery rate admitted
- logFCThr, the minimum Fold Change of difference between the two comparisons 

*Enrichment analysis*

If you want to get the enrichment of GO terms and/or pathways you must give to DEA also the filtered Annocript output

		SELECT: de_go_enrichments = YES
			de_pw_enrichments = YES

Parameters:
- p_val: Adjusted p-value cutoff to consider a class as significant (def: 0.1 means 10% are false positives)
- min_transcripts: indicates the minimum number of transcripts associated to a GO class to take it into consideration for the analysis. It is relative to the number of transcript for each class into the "signifTable" variable (def: 5)
- min_respect_to: if it it "sign", the minimum number of transcripts is from the significant table  while if it is "ann" it is from the annocript output table
- topn_toplot: number of enrichment classes to plot
- col_signifTable: column where transcripts names should be find. Files are created by DEA, you can leave this as default (1)
- mult: indicates the fold a class should be higher than the universe to analyze it (def: 1)
- prop_alt: indicates what kind of enrichment you are looking for in the "signifTable" table
To look for enrichments use 'g', for impoverishment use 'l', for both use 't'


*Venn diagrams*

Venn diagrams are built using venn.diagram R function. The R script uses the files from the differential expression analysis bringing only first column of each of these files. Only 5 files can be analyzed per time as in the subject function. DEA searches all the files which finish with the variable *de_transcripts* in the config_program.txt file (def: significant.txt). All the files found are used in the plot. To plot less samples please run the task alone and reduce the amount of files in the working folder.

		SELECT: de_venn_diagram = YES

Parameters: 
- venn_chars_to_use: the number of chars to be used as title of each set in the plot. change accordingly to graphic input
- max_files_venn: maximum number of files to analyze


*Expression plots*

Starting from a differential expression analysis between two samples, DEA can make a plot which shows for each GO term the amount of transcripts which are up or downregulated.

		SELECT expression_plots = YES

It is mandatory to specify:
- the Annocript output in input to the script (use --annocript_out)
- the UniProt database used in the Annocript output (use: ann_out_db_used[uniref/trembl])
- if GO terms should be associated to 'proteins' or 'domains'

Parameters:
ann_out_db_used: UniProt database used in Annocript
go_terms_association: where should be take GO terms (proteins or domains)
topN_exp_plots: number of GO terms to be shown
maxLengthDescs_exp_plots: maximum length of the descriptions


To start dea you can use the following command
perl  dea.pl --input_folder <input folder> --config_file <file_with_parameters> --target_file  <file with file names and sorting> [--annocript_out <annocript_output> --full_transcriptome <trinity output transcriptome>]

You may want to run it in background and get a log file:

nohup perl  dea.pl --input_folder <input folder> --config_file <file_with_parameters> --target_file  <file with file names and sorting> [--annocript_out <annocript_output> --full_transcriptome <trinity output transcriptome>] > dea_exec_info.log &
