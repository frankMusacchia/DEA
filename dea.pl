#!/usr/bin/perl
# Package DEA - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2016)
#
#This software starts from the counts files from samtools and creates a complete table of
#all the counts from which it starts the differential expression analysis

#To use it you need:

#1. A target file with an header which is a list of samples with an experimental condition code associated with.
#Example:
#		  name	condition
#			SampleA1	A
#			SampleA2	A
#			SampleB1	B
#			SampleB2	B
#
#2. A folder where you put all the counts files. The counts files should 
#	have be named: sample_name.counts. As the samples in the target file 
# were named. Example: SampleA1.counts,SampleA2.counts,SampleB1.counts,SampleB2.counts
#
#3. A configuration file with all the parameters to use for the script. See the sample file to understand.
#
#Example of usage:
# perl  dea.pl --input_folder <input folder> --config_file <file_with_parameters> --target_file  <file with file names and sorting> 
#[--annocript_out <annocript_output> --full_transcriptome <trinity output transcriptome>]\n ";

#Use the configuration file to select what you want to perform:
#
#   If you want to get a filtered transcriptome to annotate you must use the first three parameters and
#			the parameter --full_transcriptome= the path to the transcriptome created with the assembling program. 
#			SELECT : filter_transcriptome = YES
#
#	  If you want to perform a differential gene expression analysis you can use only the counts files
#			SELECT: de_analysis = YES
#
#		If you want to get the enrichment of GO terms and/or pathways you must give to DEA also the filtered Annocript output
#			SELECT: de_go_enrichments = YES
#							de_pw_enrichments = YES
#

use strict;
use warnings;
use Getopt::Long;#To use the parameters in input
use Pod::Usage;#Used to write an usage
use Data::Dumper;#To print the hashes
use File::Copy;#To manage files
use FindBin;
use lib $FindBin::Bin;

#use File::Basename qw(dirname);
#use Cwd  qw(abs_path);
#use lib dirname(dirname abs_path $0) . '/lib';
use USEFUL::utilities qw(try_exec_command extract_name checkVariable extract_fasta_from_fasta checkConfigVariables);
use LIBS::get_expression_plots qw(execute_expression_plots);

my $programName = "dea.pl";

#GLOBAL VARIABLES COMMON TO MANY PROGRAMS
my $configHash;#This hash substain all the life of the program
my $configFile = '';#The path to the config user. This must be give in input
my $configProgram = "CONFIGURATION/config_program.txt";#The path to the configuration file for the program
my $variablesFile = "CONFIGURATION/variables.txt";
my $workingFolder = '';#Folder where to work
my $programFolder = '';#Folder of the program
my $foldersFile = "folders.txt";
my $programVersion = '0.1.2';
my $readme = "README.me";
my $run = 1;#Determines if the program should execute tasks or not (used to print the help)
my $logFolder = "LOG";

#GLOBAL VARIABLES SPECIFIC FOR PROGRAM 
my $inputFolder = '';#Folder with raw reads
my $paramFile = '';#File with parameters
my $targetFile = '';#File with information about the samples
my $annocriptOut = '';#Annocript output name
my $fullTranscriptome = "";#The full transcriptome


##DECLARATION OF GLOBAL INDEXES OF THE FIELDS IN THE ANNOCRIPT OUTPUT
my $transcrId ; my $transcrLength ;

#my $count;
my $Name1; my $Length1;
my $Score1; my $HitLength1;
my $QCoverage1; my $HCoverage1; 
my $Descr1;
my $Strand1;
my $EnzId; my $EnzDescr; 
my $PwLev1;my $PwLev2;my $PwLev3;

my $Name2 ; my $Length2; 
my $Score2; my $HitLength2;
my $QCoverage2 ; my $HCoverage2;
my $Strand2;
my $Descr2;
my $closerOS;

my $BPId; my $BPDesc;
my $MFId; my $MFDesc;
my $CCId; my $CCDesc;

my $domBPId; my $domBPDesc;
my $domMFId; my $domMFDesc;
my $domCCId; my $domCCDesc;

my $CDName; my $CDStartEnd;
my $CDScore ; my $CDDescr;

my $riboName; my $riboScore;
my $riboDesc;
my $longORFLength; my $longORFS; my $longORFF;

my $probNC; my $NC4Annocript; my $seq;

#Defining an hash to keep all the indexes of the columns in the Annocript
# output as they are found
my $headHash;

#LOG nice printings
my $niceClosure = "\n=================================================\n";
my $configClosure = "\n*************************************************\n";


=head2 settingsForUniref

 Title : settingsForUniref
 Usage : settingsForUniref(  );

 Function: this function sets the indexes of the array  containing Annocript output for the case of using the
            uniref

 Returns : nothing

=cut
sub setFieldNames {
  my $outFileName = shift;
  my $dbInUse = shift;
  
  #Open the output file and read the header
  open(OUT_ANN,"<$outFileName") or die "Cannot open the file $outFileName: ?\n";
  my $row = <OUT_ANN>;
  chomp($row);
  #print "Just read the first line of $outFileName: ".$row;
  #Put the headers name in an hash and associate them the field number
  my @fields = split(/\t/, $row);
  #print_array(\@fields);#DEBUGCODE
  #print "Number of fields :".scalar(@fields)."\n";#DEBUGCODE
  my $numField = 0;
  #Here we fill the hash with the field's names and assign a number
  foreach my $field (@fields){
	  $headHash->{$field} = $numField;
	  #print "Assigning $numField to the field '$field': ".$headHash->{$field}."\n";
	  $numField++;
	  
  }	
  close(OUT_ANN);
    
  ##INDEXES OF THE FIELDS IN THE ANNOCRIPT OUTPUT
  ## IF THE FIELDS CHANGE THEN YOU MUST CHANGE THE NAMES ACCORDINGLY
  ## SEE ALSO HEADERS IN GFF3_MANAGER.pm
  $transcrId = 'TranscriptName'; $transcrLength = 'TransLength';

  $Name1 = 'HSPNameSP'; $Length1 = 'HSPLengthSP'; 
  $Score1 = 'HSPEvalueSP'; $HitLength1 = 'HITLengthSP';
  $QCoverage1 = 'QCoverageSP';  $HCoverage1 = 'HCoverageSP'; 
  $Strand1 = 'StrandSP';
  $Descr1 = 'DescriptionSP';   
  $EnzId = 'EnzymeIds';  $EnzDescr = 'EnzymeDescs'; 
  $PwLev1 = 'PwLev1'; $PwLev2 = 'PwLev2'; $PwLev3 = 'PwLev3';
  
  if ($dbInUse eq "uniprotkb"){ 
		$Name2 = 'HSPNameTR'; $Length2 = 'HSPLengthTR' ; 
		$Score2 = 'HSPEvalueTR'; $HitLength2 = 'HITLengthTR';
		$QCoverage2 = 'QCoverageTR';  $HCoverage2 = 'HCoverageTR'; 
		$Strand2 = 'StrandTR';
		$Descr2 = 'DescriptionTR'; 
		$closerOS = 'OSName';
  }
  if ($dbInUse eq "uniref"){ 
		$Name2 = 'HSPNameUf'; $Length2 = 'HSPLengthUf' ; 
		$Score2 = 'HSPEvalueUf'; $HitLength2 = 'HITLengthUf';
		$QCoverage2 = 'QCoverageUf';  $HCoverage2 = 'HCoverageUf'; 
		$Strand2 = 'StrandUf';
		$Descr2 = 'DescriptionUf'; 
		$closerOS = 'Taxonomy';
	}   
 
  $CDName = 'CDName'; 
  $CDStartEnd = 'CDStartEnd';
  $CDScore = 'CDEvalue';  $CDDescr = 'CDDesc';

  $BPId = 'BPId';  $BPDesc = 'BPDesc';
  $MFId = 'MFId';  $MFDesc = 'MFDesc';
  $CCId = 'CCId';  $CCDesc = 'CCDesc';

  $domBPId = 'domBPId';  $domBPDesc = 'domBPDesc';
  $domMFId = 'domMFId';  $domMFDesc = 'domMFDesc';
  $domCCId = 'domCCId';  $domCCDesc = 'domCCDesc';

  $riboName = 'OtherNCName'; $riboScore = 'OtherNCEvalue'; $riboDesc = 'OtherNCDesc'; 
  $longORFLength = 'LongOrfLength';
  $longORFS = 'LongOrfStrand';  $longORFF = 'LongOrfFrame';

  $probNC = 'ProbToBeNonCoding';  $NC4Annocript = 'lncRNA4Annocript'; $seq = 'Sequence';

}


=head2 checkConfigVariables

 Title   : checkConfigVariables
 Usage   : checkConfigVariables( - configUser -> file with the user configuration
                              - configAnnocript -> file with the basic parameters of Annocript
                              - variablesFile -> the path to a file with all variables written
          );

 Function: this subroutine reads the config files and check if all variables are there and are well written.
            The variables.txt file is needed fot this operation.
 
 Returns : nothing

=cut
sub checkSingleConfigVariables {
  
  my $hashUserCheck;
          
  if (! open(VARF,"<".$programFolder."/".$variablesFile)){ die "ERROR: Failure opening ".$programFolder."/".$variablesFile.". Your $programName version is corrupted - $!";}
  if (! open(CUSER,"<$configFile")){ die "ERROR: Cannot find '$configFile' ";}
  
  #Stores the variables in the config user file inside the hash
  my $start = 0;
  while (my $line = <CUSER>){ 
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
     #print $line."\n";#DEBUGCODE
     $hashUserCheck->{$1} = "OK";
    }
  }	

  close(CUSER);
    
  my @userConfVars = ();  
  
  #Variables that are in the configuration file must be also in the variables.txt file
  my $errors=0;
  my $lines=0;
  
  my $line = <VARF>;
  
  $line =~ s/\n//;#Remove \n in the end of the line
    
  #get the variables in the line
  my @variables = split (/;/,$line);
    
  $lines++;
  #For each of the variables in the variables file
  foreach my $var (@variables){
    #print "Variable: $var - - value: ".$hashUserCheck->{$var}."\n";#DEBUGCODE
    #put the variable inside an array for ann config
    push (@userConfVars, $var);
    if( !(defined($hashUserCheck->{$var})) ){
      die "ERROR: in $configFile variable $var is missing. Please check the file. Closing...\n "; 
      $errors=1;
    }#else{ annoPrint ("From the hash: ".$hashCheck->{$var}."\n";)}
  }
   
  #print_array(\@allVars);	
  #print Dumper\$hashCheck;
  #Now check if all the elements in the hash are also in the array
  foreach my $key (keys %$hashUserCheck){
     # print "Search $key in array...\n";#DEBUGCODE
     if (!(grep {/$key/} @userConfVars )){
       die "ERROR: Variable $key is in the user config files and not in $variablesFile file. This is completely wrong. $programName will not work...\n ";
     }
  }
 
  close(VARF);
}

=head2 configFile2Hash

 Title   : configFile2Hash
 Usage   : configFile2Hash( - configFilePath = path of the config file
                               );

 Function:  gets the hash table with all the path and names in input from the config file in input
 Returns : nothing

=cut
sub configFile2Hash{  
  my $configFilePath=shift;
  
  my $start = 0;
  #Here we open config file and read all its line to find elements belonging to each of the executers
  open (configFile,$configFilePath) or annoDie("ERROR: The file $configFilePath doesn't exists. Annocript will exit..\n");
  while (my $line = <configFile>){ 
    if ($line =~ /#########/){
      $start = 1;
    }
    if( ($line =~ /(\S+)\s*=\s*(\S+)/) and ($start == 1) and !($line =~ /#/) ){
      $configHash->{$1} = $2; 
      #annoPrint ("$1 = $2\n") ;#DEBUGCODE     
    }
  }	 
  close(configFile);	
	#annoPrint (Dumper\$configHash); #DEBUGCODE	
}


  
=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.
  my $VERSION = 0;# Shows version number and exit.
  
  
  my $howToUse = "$programName$programVersion the software for the gene expression analysis using Annocript.\n".
  "Please use with: perl  $programName --input_folder <input folder> --config_file <file_with_parameters> --target_file ".
  " <file with file names and sorting> [--annocript_out <annocript_output> ]\n ";
  #  Parse options
  GetOptions(
           "help|h" =>        \$HELP,
           "version|v" =>      \$VERSION,
           "input_folder|i=s" =>  \$inputFolder,#It's mandatory 
           "config_file|c=s" => \$configFile,#It's mandatory 
           "target_file|t=s" => \$targetFile,   #It's mandatory
           "annocript_output|a=s" =>  \$annocriptOut,   #It's mandatory for some     
           "full_transcriptome|f=s" => \$fullTranscriptome #Full transcriptome to use fo filtering
          );
  
  #Let's create an array with mandatory parameters depending from the execution
  #for each program executed, it will ask for certain parameters. Finally, they will be checked
  my @mandatoryPar =  ();
  push(@mandatoryPar,'target_file','input_folder','config_file'); #These parameters are always needed
          
  #Print a little help  
  if ( $HELP ){
		#print $howToUse;
		$run = 0;
		#pod2usage(1);
    #exit;
  }
  
  #Print version
  if ( $VERSION ){
    print "$programName Version $programVersion \n";
    exit;
  }

	#If the program should execute some task do these checks
  if ($run){
		print "Input folder: $inputFolder\n";
		#Check if the folder with reads exists
		if ( $inputFolder ne '' ){
				if (!(-d $inputFolder) ){
					print "Please give an existing input folder!\n $howToUse";
					exit; #Without folder... nothing. It should exist!
				}else{
						if ($inputFolder =~ /\/$/){
							chop($inputFolder);
						}
					}
		}elsif( grep {/\binput_folder\b/} @mandatoryPar){
			print "I need a folder with input files. Please give one existent!\n $howToUse";
				exit; #Without reads... nothing...
		}
	 
		 #Check if the target file exists and is needed
		if ( $targetFile ne '' ){
				if (!(-e $targetFile) ){
					print "I need a target file. Please give one existent!\n $howToUse";
					exit; #Without target_file... nothing. It should exist!
				}else{
			 #Move the target file in the input folder and get the name only
			 copy($targetFile,$inputFolder);
			 $targetFile = extract_name($targetFile,0);#Get only the name 
			}
		}elsif ( grep {/\btarget_file\b/} @mandatoryPar ) {
			print "I need a target file. Please give one existent!\n $howToUse";
				exit; #Without target_file... nothing...
		}
		
		#Check if the configuration file exists and is needed
		if ( $configFile ne '' ){
				if (!(-e $configFile) ){
					print "I need a configuration file. Please give one existent!\n $howToUse";
					exit; #Without configFile... nothing. It should exist!
				}
		}elsif ( grep {/\bconfig_file\b/} @mandatoryPar ) {
			print "I need a configuration file. Please give one existent!\n $howToUse";
				exit; #Without target_file... nothing...
		}
	}
 }

=head2 check_configuration

 Title   : check_configuration
 Usage   : check_configuration(    );

 Function: checks integrity of variables in the configuration file and dependencies
 Returns : nothing

=cut
sub check_configuration{ 
	
	
  #Programs to execute
  #Check the filter transcriptome program
  #needs the files with counts and the full transcriptome as from the assembly program, target file (input)
  print $configClosure;
  checkVariable($configHash->{'filter_transcriptome'},'filter_transcriptome', "counts file will be used to filter the transcriptome.\n");
  if ($configHash->{'filter_transcriptome'} eq 'YES'){
		if ( $fullTranscriptome eq ''){
					 print "When you want to filter the full transcriptome, full_transcriptome parameter must be used. Exiting..\n";
					 exit;
				}elsif (!(-e $fullTranscriptome) or (-z $fullTranscriptome) ){
					print "$fullTranscriptome does not exist. Please execute the assembly program and get the transcriptome. Exiting..\n";
					 exit;
			}else{
				#Move the annocript output in the input folder and get the name only
				copy($fullTranscriptome,$inputFolder);
				$fullTranscriptome = extract_name($fullTranscriptome,0);
			}  
			#Check the presence of the counts files using the target file
			check_counts();
  }
  
	#Check the de analysis program used
	#needs: counts filtered (produced by filter transcriptome),target file (input), annocript output (to get with Annocript)
	checkVariable($configHash->{'de_analysis'},'de_analysis', "differential expression analysis wil be performed.\n");
	if ($configHash->{'de_analysis'} eq 'YES'){
			if ( $annocriptOut eq ''){
					print "annocript_out parameter has not been used. DEA will not attach annotation to the de analysis..\n";
					#exit;
			}elsif (!(-e $annocriptOut) or (-z $annocriptOut) ){
					print "$annocriptOut does not exist. Please execute Annocript and get the filtered annotation table. Exiting..\n";
					exit;
			}else{
				#Move the annocript output in the input folder and get the name only
				copy($annocriptOut,$inputFolder);
				$annocriptOut = extract_name($annocriptOut,0);
			}
			my $counts = $inputFolder."/".$configHash->{'organism'}."_counts_filtered.txt";
			if (! (-e $counts) and ($configHash->{'filter_transcriptome'} eq 'NO') ){
				print "There should be a file with counts inside the input folder: $inputFolder ".
							"named $counts. It is missing! Please chech the organism name in the configuration ".
							"file or run filter_transcriptome before.Exiting...\n";
				exit;
			}
	}
		
  
  #Programs which use the Annocript output will not be xhecked/executed if the user starts the transcriptome filtering
  if ($configHash->{'filter_transcriptome'} eq 'NO'){
		
		
		#Check the expression plots variable used
		#needs: counts filtered (produced by filter transcriptome),target file (input), annocript output (to get with Annocript)
		checkVariable($configHash->{'expression_plots_GO'},'expression_plots', "Construction of expression plots for GO terms will be performed.\n");
		checkVariable($configHash->{'expression_plots_PW'},'expression_plots', "Construction of expression plots for pathways will be performed.\n");
		if ($configHash->{'expression_plots_GO'} eq 'YES' or $configHash->{'expression_plots_PW'} eq 'YES'){
			
			if ( !($configHash->{'ann_out_db_used'} eq 'uniref') and !($configHash->{'ann_out_db_used'} eq 'trembl') ){
				die "ERROR! Variable ann_out_db_used in config_user file can be only 'uniref' or 'trembl'";
			}
			if ( !($configHash->{'go_terms_association'} eq 'proteins') and !($configHash->{'go_terms_association'} eq 'domains') ){
				die "ERROR! Variable go_terms_association in config_user file can be only 'proteins' or 'domains'";
			}
			if ( $annocriptOut eq ''){
						print "annocript_out parameter has not been used. DEA needs the output from Annocript to perform this taks...\n";
						exit;
				}elsif (!(-e $annocriptOut) or (-z $annocriptOut) ){
						print "$annocriptOut does not exist. Please execute Annocript and get the filtered annotation table. Exiting..\n";
						exit;
				}else{
					if( ( ($configHash->{'ann_out_db_used'} eq 'uniref') or ($configHash->{'ann_out_db_used'} eq 'trembl') ) ){
					#Move the annocript output in the input folder and get the name only
					copy($annocriptOut,$inputFolder);
					$annocriptOut = extract_name($annocriptOut,0);
					}else{
							print "Please use the parameter ann_out_db_used in the configuration file to specify what kind of database was used in the output".
							" of Annocript [uniref|trembl]. Exiting...\n";
							exit;
					}
				}
				my $counts = $inputFolder."/".$configHash->{'organism'}."_".$configHash->{'counts_filt_suffix'};
				if (! (-e $counts) and ($configHash->{'filter_transcriptome'} eq 'NO') ){
					print "There should be a file with counts inside the input folder: $inputFolder ".
								"named $counts. It is missing! Please chech the organism name in the configuration ".
								"file or run filter_transcriptome before. Exiting...\n";
					exit;
				}
		}		
		
		  #Check the de analysis program used
  #needs: counts filtered (produced by filter transcriptome),target file (input), annocript output (to get with Annocript)
  checkVariable($configHash->{'de_go_enrichments'},'de_go_enrichments', "GO enrichment analysis of differentially expressed will be performed.\n");
  if ($configHash->{'de_go_enrichments'} eq 'YES'){
	   if ( $annocriptOut eq ''){
          print "annocript_out parameter must be used. Exiting..\n";
          exit;
      }elsif (!(-e $annocriptOut) or (-z $annocriptOut) ){
		 print "$annocriptOut does not exist. Please execute Annocript and get the filtered annotation table. Exiting..\n";
          exit;
	  }else{
		#Move the annocript output in the input folder and get the name only
		copy($annocriptOut,$inputFolder);
		$annocriptOut = extract_name($annocriptOut,0);
	  }
	  
	 #If de analysis is not executed, here it should control if some file with enrichments is present
	 #print "\nTo build the plots all the DE analysis files in the $inputFolder directory will be used. To get a different diagram please ".
	#"remove from this folder the files you are not interested in.\n";
	chdir($inputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
	my @de_files = glob("*".$configHash->{'de_transcripts'});
	chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
	
	if ( ($configHash->{'de_analysis'} eq 'NO') and (scalar (@de_files) == 0)){
		print "No differential expression results in the input folder: $inputFolder. Please run de_analysis\n";
		exit;
	}
	 $configHash->{'gomap'}  = $programFolder."/".$configHash->{'gomap'};
  }
  
  checkVariable($configHash->{'de_pw_enrichments'},'de_pw_enrichments', "Pathways enrichment analysis of differentially expressed will be performed.\n");
  if ($configHash->{'de_pw_enrichments'} eq 'YES'){
	  if ( $annocriptOut eq ''){
          print "annocript_out parameter must be used. Exiting..\n";
          exit;
      }elsif (!(-e $annocriptOut) or (-z $annocriptOut) ){
		 print "$annocriptOut does not exist. Please execute Annocript and get the filtered annotation table. Exiting..\n";
          exit;
	  }else{
		#Move the annocript output in the input folder and get the name only
		copy($annocriptOut,$inputFolder);
		$annocriptOut = extract_name($annocriptOut,0);
	  }
		#If de analysis is not executed, here it should control if some file with enrichments is present
		print "\nTo build the plots all the DE analysis files in the $inputFolder directory will be used. To get a different diagram please ".
		"remove from this folder the files you are not interested in.\n";
		chdir($inputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
		my @de_files = glob("*".$configHash->{'de_transcripts'});
		chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
		
		if ( ($configHash->{'de_analysis'} eq 'NO') and (scalar (@de_files) == 0)){
			print "No differential expression results in the input folder: $inputFolder. Please run de_analysis\n";
			exit;
		}
  }
		
	}
		
	checkVariable($configHash->{'de_venn_diagram'},'de_venn_diagram', "Venn diagrams for differentially expressed genes will be generated.\n");
  if ($configHash->{'de_venn_diagram'} eq 'YES'){
		print "\nTo build the Venn diagram all the DE analysis files in the $inputFolder directory will be used. To get a different diagram please ".
		"remove from this folder the files you are not interested in.\n";
		chdir($inputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
		my @de_files = glob("*".$configHash->{'de_transcripts'});
		chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
		
		if ( ($configHash->{'de_analysis'} eq 'NO') and (scalar (@de_files) == 0)){
			print "No differential expression results in the input folder: $inputFolder. Please run de_analysis\n";
			exit;
		}elsif (scalar (@de_files) > $configHash->{'max_files_venn'}){
			print "There is a limit on the files to use for the Venn diagram. It is ".$configHash->{'max_files_venn'}.". The diagram will be ".
			"generated with the first $configHash->{'max_files_venn'} files.\n";
			}
  }
  
  checkVariable($configHash->{'get_up_down'},'get_up_down', "Up and down regulated will be separated from differentially expressed genes.\n");
  if ($configHash->{'get_up_down'} eq 'YES'){
		 my $updownfolder = $configHash->{'up_down_folder'};
		 
		if ( ($inputFolder  =~ /$updownfolder/)){
			print "Up and down regulated will not be generated since you are in: $inputFolder.\n";
			$configHash->{'get_up_down'} = 'NO';
		}else{
			print "All the DE analysis files in the $inputFolder directory will be used. To get a different analysis please ".
			"remove from this folder the files you are not interested in.\n";
			chdir($inputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
			my @de_files = glob("*".$configHash->{'de_transcripts'});
			chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
			
			if ( ($configHash->{'de_analysis'} eq 'NO') and (scalar (@de_files) == 0)){
				print "No differential expression results in the input folder: $inputFolder. Please run de_analysis\n";
				exit;
			}
		}
  }
  print $configClosure;
}

=head2 check_counts

 Title   : check_counts
 Usage   : check_counts(    );

 Function:  checks the counts files given in the input folder comparing their names with the samples names in the target file
 Returns : 1 is the check is good, -1 otherwise

=cut
sub check_counts{ 
  my $nRows = 1;
  my $retVal = 0;
  
  print "$targetFile\n";
  open (TARG, "<$inputFolder/$targetFile") or die "ERROR: cannot open $targetFile";
  print "Checking counts files..";
  #Opens the dir with all the databases
  opendir DIR, $inputFolder or die  "ERROR : cannot open dir $inputFolder";
  my @files= readdir DIR;#All the files are in this list
  #print "Found ".scalar(@files)."files\n";	
	while (my $row = <TARG>){
		#print $row."\n";
		if ($nRows > 1){
			my @pieces = split("\t",$row);
			my $sample = $pieces[0].".counts";
			if (!(grep {/\b$sample\b/} @files)){
				 die "$sample not present in $inputFolder but needed. Check the target file and the input folder.\n";
			}
		}
		$nRows++;	
	}
  close(TARG);
  print "OK!\n";
}


 
=head2 execute_analyses

 Title   : execute_analyses
 Usage   : execute_analyses(    );

 Function:  execute all the analyses indicated in the configuration file
 Returns : nothing

=cut
sub execute_analyses{ 
	my $rbatch = "R CMD BATCH  --no-save --no-restore '--args ";
	my $command = '' ;
	my $programsPath = "PROGRAMS/";
	
	#Set the LOG folder
	$logFolder =  $inputFolder."/".$logFolder;
	#Creates the log folder if it does not exists
	unless(-d $logFolder){
		print $logFolder." doesn't exist. Creating folder...\n";
		mkdir $logFolder or die "ERROR: can't create folder $logFolder\n ";
	}
	
	#STEP (A) merge the counts columns in a unique file AND 
	#STEP (B) filter your assembled transcriptome on the base of the expression
	if ($configHash->{'filter_transcriptome'} eq 'YES'){
		print $niceClosure;
		$configHash->{'counts'} = $configHash->{'organism'}."_".$configHash->{'counts'};
		$configHash->{'transcripts_filtered_list'} = $configHash->{'organism'}."_".$configHash->{'transcripts_filtered_list'};
		$configHash->{'filtered_transcriptome'} = $configHash->{'organism'}."_mincpm".$configHash->{'min_cpm'}."per".$configHash->{'min_num_repl'}.".fasta";
		print "Getting a file with selected transcripts. Functions which use the Annocript output will not be run...\n";
		$command = $rbatch.
		" $inputFolder ".$configHash->{'column_with_count'}." $targetFile ".$configHash->{'organism'}." ".$configHash->{'min_cpm'}." ".$configHash->{'min_num_repl'}."' ".
		$programFolder."/".$programsPath.$configHash->{'get_counts_and_filter_sw'};		
	   	print "Executing command: $command\n";
	   	try_exec_command($command) or die "Unable to execute command: $command\n";	
	    extract_fasta_from_fasta($fullTranscriptome,$inputFolder."/".$configHash->{'transcripts_filtered_list'},$inputFolder."/".$configHash->{'filtered_transcriptome'});
		print "$programName produced a fasta file with all the transcripts that have a good number of reads that mapped. Specifically\n".
		" all the transcripts with cpm > ".$configHash->{'min_cpm'}." for at least ".$configHash->{'min_num_repl'}." samples, were taken.";
		print "Please use now Annocript to get the annotation and come back to $programName to get the differential expression analysis of samples.\n";
		exit;#At the end of this process nothing more can be performed
	}
	
	
	
	#STEP (C) perform the differential expression (DE) analysis of the samples using the file of counts
	if ($configHash->{'de_analysis'} eq 'YES'){
		print $niceClosure;
		print ">>Getting the differential expression analysis\n ";
		$configHash->{'counts_filt_file'} = $configHash->{'organism'}."_".$configHash->{'counts_filt_suffix'};
	  my $parameters = '';
	  
	  if ( $annocriptOut ne ''){
			$parameters = " $inputFolder ".$configHash->{'counts_filt_file'}." ".$configHash->{'organism'}." $targetFile  ".$configHash->{'minFDR'}."  ".$configHash->{'de_transcripts'}." 1 ".$annocriptOut."' ";
		}else{
			$parameters = " $inputFolder ".$configHash->{'counts_filt_file'}." ".$configHash->{'organism'}." $targetFile  ".$configHash->{'minFDR'}."  ".$configHash->{'de_transcripts'}." 1' ";
		}	
		
		my $log = $logFolder."/de_analysis.Rout";
	  $command = $rbatch.$parameters.$programFolder."/".$programsPath.$configHash->{'de_analysis_sw'}." $log";		
	  print "Executing command: $command\n";
	  try_exec_command($command) or die
			"ERROR [$!]: an error occurred while executing $command. Please check the R log file: $log\n";	
	}
	
	if ($configHash->{'filter_transcriptome'} eq 'NO'){
		
		#STEP(G) get the expression plots of most differentially expressed GO terms
		#Needs the file of counts and the output from Annocript	
		if ($configHash->{'expression_plots_GO'} eq 'YES'){
			print $niceClosure;
			getExpressionPlots('go',$rbatch,$programsPath);
		}
		#STEP(G) get the expression plots of most differentially expressed GO terms
		#Needs the file of counts and the output from Annocript	
		if ($configHash->{'expression_plots_PW'} eq 'YES'){
			print $niceClosure;
			getExpressionPlots('pathways',$rbatch,$programsPath);
		}
		#At the end of this process nothing more can be performed
		
		if ($configHash->{'expression_plots_GO'} eq 'YES' or $configHash->{'expression_plots_PW'} eq 'YES'){
			
			exit;
			print $niceClosure;
		}
		
		#STEP (E) get GO and pathways enrichments TABLE and PLOTS of the DE transcripts using the output from Annocript
	if ($configHash->{'de_go_enrichments'} eq 'YES'){
		print $niceClosure;
		print "\n\n>>Getting the GO terms enriched from your differential expression analysis\n ";
		#Get all the files of the differential expression analysis and execute the GO and Pathways enrichment
		# Go into the directory containing the  files
		chdir($inputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
		my @de_files = glob("*".$configHash->{'de_transcripts'});
		chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
		
		if (scalar(@de_files) > 0){
			foreach my $de_file (@de_files){
				print "Enrichment for $de_file\n";
				
				my $log = $logFolder."/de_go_enrichments.Rout";#Log file
				$command = $rbatch.
				" $inputFolder $annocriptOut $de_file ".$configHash->{'organism'}." ".$configHash->{'p_val'}." ".$configHash->{'min_transcripts'}.
				" ".$configHash->{'min_respect_to'}." ".$configHash->{'gomap'}."' ".
				$programFolder."/".$programsPath.$configHash->{'go_enrichments_sw'};		
				print "Executing command: $command\n";
				try_exec_command($command) or die "ERROR [$!]: an error occurred while executing $command. Please check the R log file: $log\n";
			}
		}else{print "No significant transcripts to analyze here!\n";}
		print $niceClosure;
	}
	# STEP (E) get GO and pathways enrichments TABLE and PLOTS of the DE transcripts using the output from Annocript
	if ($configHash->{'de_pw_enrichments'} eq 'YES'){
		print $niceClosure;
		print "\n\n>>Getting the pathways enriched from your differential expression analysis\n ";
		#Get all the files of the differential expression analysis and execute the GO and Pathways enrichment
		# Go into the directory containing the files
		chdir($inputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
		my @de_files = glob("*".$configHash->{'de_transcripts'});
		chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
		if (scalar(@de_files) > 0){
			foreach my $de_file (@de_files){
				print "Enrichments for $de_file\n";
				
				my $log = $logFolder."/de_pw_enrichments.Rout";#Log file
				$command = $rbatch.
				" $inputFolder $annocriptOut $de_file ".$configHash->{'organism'}." ".$configHash->{'p_val'}." ".$configHash->{'min_transcripts'}.
				" ".$configHash->{'min_respect_to'}."' ".
				$programFolder."/".$programsPath.$configHash->{'pw_enrichments_sw'};		
				print "Executing command: $command\n";
				try_exec_command($command) or die "ERROR [$!]: an error occurred while executing $command. Please check the R log file: $log\n";
			}
		}else{print "No significant transcripts to analyze here!\n";}
		print $niceClosure;
	}
	}
		
	
	
	#STEP (D) extract up and down regulated transcripts from each DE analysis
	if ($configHash->{'get_up_down'} eq 'YES'){
		print $niceClosure;
		print "\n\n>>Getting the up and down regulated from your differential expression analysis\n ";
		#Get all the files of the differential expression analysis and execute the GO and Pathways enrichment
		# Go into the directory containing the files
		chdir($inputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
		my @de_files = glob("*".$configHash->{'de_transcripts'});
		chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
		if (scalar(@de_files) > 0){
			my $all_de_files = '';
			foreach my $de_file (@de_files){
					$all_de_files .= $de_file." ";#build a string with all the files to use
			}
			print "Generating up and down\n";
			my $log = $logFolder."/get_up_down.Rout";#Log file
				$command = $rbatch.
				" $inputFolder ".$configHash->{'up_down_folder'}." $all_de_files' ".
				$programFolder."/".$programsPath.$configHash->{'get_up_down_sw'};		
				print "Executing command: $command\n";
				try_exec_command($command) or die "ERROR [$!]: an error occurred while executing $command. Please check the R log file: $log\n";
		}else{print "No significant transcripts to analyze here!\n";}
		print $niceClosure;
	}
	
	#(F) extract a Venn diagram showing intersections among the DE analyses
	if ($configHash->{'de_venn_diagram'} eq 'YES'){
		print $niceClosure;
		print "\n\n>>Generating Venn diagrams for differentially expressed genes\n ";
		#Get all the files of the differential expression analysis and generate a Venn diagram
		# Go into the directory containing the  files
		chdir($inputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
		my @de_files = glob("*".$configHash->{'de_transcripts'});
		chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
		
		if (scalar(@de_files) > 0){
			my $all_de_files = '';
			my $de_files_num = 1;
			foreach my $de_file (@de_files){
				if ($de_files_num <= $configHash->{'max_files_venn'}){
					$all_de_files .= $de_file." ";#build a string with all the files to use
				}
				$de_files_num++;
			}
			print "Venn diagram\n";
			my $log = $logFolder."/de_venn_diagram.Rout";#Log file
				$command = $rbatch.
				" $inputFolder ".$configHash->{'venn_chars_to_use'}." $all_de_files' ".
				$programFolder."/".$programsPath.$configHash->{'venn_diagrams_sw'};		
				print "Executing command: $command\n";
				try_exec_command($command) or die "ERROR [$!]: an error occurred while executing $command. Please check the R log file: $log\n";
		}else{print "No significant transcripts to analyze here!\n";}
		print $niceClosure;
	}
	
	if ($configHash->{'extr_ph_expressed'} eq 'YES'){
	  		
	}
	
	if ($configHash->{'extr_ph_specific'} eq 'YES'){
	  		
	}
	
	if ($configHash->{'ph_enrichments'} eq 'YES'){
	  		
	}
	
}


sub getExpressionPlots {
		my $dataType = shift;
		my $rbatch = shift;
		my $programsPath = shift;
				
		#Creating a supplementary output folder for these plots
		my $outputFolder = $inputFolder."/".$configHash->{'expression_plots_folder'};
		unless (-e $outputFolder) { mkdir $outputFolder or die "Cannot create folder $outputFolder check permissions...\n";}
		#Getting the expression tables
		print ">>Getting the expression tables...\n ";
		$configHash->{'counts_filt_file'} = $configHash->{'organism'}."_".$configHash->{'counts_filt_suffix'};
	  my $parameters = '';
	  
	  #We copy these two files (they could be removed later, because they are in the session folder, too
	  copy($inputFolder."/".$configHash->{'counts_filt_file'},$outputFolder."/".$configHash->{'counts_filt_file'});
	  copy($inputFolder."/".$targetFile,$outputFolder."/".$targetFile);
	  
	  #Add the output of Annocript to the list of parameters
	  if ( $annocriptOut ne ''){
			$parameters = " $outputFolder ".$configHash->{'counts_filt_file'}." ".$configHash->{'organism'}." $targetFile  ".$configHash->{'minFDR'}." ".$configHash->{'transcripts_exp'}." 2 ".$annocriptOut."' ";
		}else{
			die "Unable to execute the analysis without the output from Annocript\n";		
		}	
		#Creating a command and executing the R script
	  my $command = $rbatch.$parameters.$programFolder."/".$programsPath.$configHash->{'de_analysis_sw'};		
	  print "Executing command: $command\n";
			try_exec_command($command) or die "Unable to execute command: $command\n";
	  
	  print ">>Getting the tables with the abundances...\n ";	
	  #Set the names of the columns from the Annocript output
	  setFieldNames($inputFolder."/".$annocriptOut,$configHash->{'ann_out_db_used'});
	  
	  #Decide which indexes to use in the table depending by the GO terms wanted
	  my $indexes;
	  #Suffixes to save and search files
		my $suffixes;
		
	  if ($dataType  eq 'go'){
			$suffixes->{'index1'} = $configHash->{'bp_suffix'};
			$suffixes->{'index2'} = $configHash->{'mf_suffix'};
			$suffixes->{'index3'} = $configHash->{'cc_suffix'};
			if ( $configHash->{'go_terms_association'}  eq 'proteins'){
				$indexes->{'index1'} = $headHash->{$BPId};
				$indexes->{'index2'} = $headHash->{$MFId};
				$indexes->{'index3'} = $headHash->{$CCId};
			}elsif ( $configHash->{'go_terms_association'}  eq 'domains'){
					$indexes->{'index1'} = $headHash->{$domBPId};
					$indexes->{'index2'} = $headHash->{$domBPId};
					$indexes->{'index3'} = $headHash->{$domCCId};
				}
		}elsif ($dataType  eq 'pathways'){
				$indexes->{'index1'} = $headHash->{$PwLev1};
				$indexes->{'index2'} = $headHash->{$PwLev2};
				$indexes->{'index3'} = $headHash->{$PwLev3};
				$suffixes->{'index1'} = $configHash->{'pwl1_suffix'};
				$suffixes->{'index2'} = $configHash->{'pwl2_suffix'};
				$suffixes->{'index3'} = $configHash->{'pwl3_suffix'};
			}
		$indexes->{'tr_id_index'} = $headHash->{$transcrId};	
		
	  #Open the output folder and search the tables coming from the t-test
	  #It will search all those files ending with $configHash->{'transcripts_exp'}
	  chdir($outputFolder) or die "\nCannot change in directory $outputFolder: $!\n\n";
		my @de_files = glob("*".$configHash->{'transcripts_exp'});
		chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
		if (scalar(@de_files) > 0){
			foreach my $de_file (@de_files){
				print ">>>>Executing the expression plot creation for $de_file\n";
				print "Indexes: ".$indexes->{'index1'}." ".$indexes->{'index2'}." ".$indexes->{'index3'}."\n" ;
				execute_expression_plots($inputFolder,$annocriptOut,$de_file,$configHash,$outputFolder,$indexes,$suffixes,$dataType);
			}
		}else{print "ERROR: no expression analysis files found in $outputFolder...\n";}	
		
		#The previous function created new tables with the abundances of expression
		#relatively to bp, mf and cc by using the same name as the expression tables and adding
		#a suffix. This suffix is searched here to plot all these tables
	  print "Printing plots for biological processes...\n";
	  chdir($outputFolder) or die "\nCannot change in directory $inputFolder: $!\n\n";
		my @files1 = glob("*".$suffixes->{'index1'});
		my @files2 = glob("*".$suffixes->{'index2'});
		my @files3 = glob("*".$suffixes->{'index3'});
		chdir($workingFolder) or die "\nCannot change in directory $workingFolder: $!\n\n";
		#all the files name will be contained in an array
		my @all_cl_files = ();
		push(@all_cl_files,@files1);
		push(@all_cl_files,@files2);
		push(@all_cl_files,@files3);
		if (scalar(@all_cl_files) > 0){
			#Print a barplot for each file
			foreach my $all_cl_file (@all_cl_files){
				$parameters = " $outputFolder $all_cl_file ".$configHash->{'organism'}." $targetFile  ".$configHash->{'minFDR'}." NULL 3 ".$annocriptOut."' ";
				$command = $rbatch.$parameters.$programFolder."/".$programsPath.$configHash->{'de_analysis_sw'};		
				print "Executing command: $command\n";
				try_exec_command($command) or die "Unable to execute command: $command\n";
			}
		}else{print "ERROR: No file of expression abundances found $outputFolder...\n";}	
	
}
=head2 initialize_folders

 Title   : initialize_folders
 Usage   : initialize_folders(  );

 Function:  this subroutine initializes the basic folders. If the folders.txt file is not present 
            means that the program has not been installed in the folder used.
 
 Returns : nothing

=cut
sub initialize_folders{

  #If the file with folders references is installed the program can start
  open (FOLD, $foldersFile) or die " ERROR: $programName has not been installed here. Please execute install.pl in the program folder.\n";
  my $line = <FOLD>;
  
  #Extracts folders paths
  my @folders = split(" ", $line);
  
  #The working folder is the first path
  $workingFolder = $folders[0];
  
  #Cancel the final slash..if it is there it is removed
  if($workingFolder =~ /\/$/){
    chop($workingFolder);
  }

  $programFolder = $folders[1];

  close(FOLD);
  

}

#Reads a file using the more command
sub read_and_print {
		my $file = shift;
		my $command = "more $file";
		try_exec_command($command) or die "Unable to execute command: $command\n";
}

###MAIN
#Get the parameters in input and does a first check
parse_command_line_args();
#get the variables indicating what are the fundamental folders
initialize_folders();

#To read an help I used this if
if (!$run){
	read_and_print($programFolder."/".$readme);
	exit;#Exit the program
}

#Get and check variables from the configuration files
checkConfigVariables($configFile,$programFolder."/".$configProgram,$programFolder."/".$variablesFile);
configFile2Hash($configFile);
configFile2Hash($programFolder."/".$configProgram);

#print Dumper\$configHash; #DEBUGCODE

#Copy the configuration file inside the working folder
copy($configFile,$inputFolder."/".extract_name($configFile,0));

#Check the variables to see the configuration the user chosen
check_configuration();

#Execute the chosen programs
execute_analyses();
