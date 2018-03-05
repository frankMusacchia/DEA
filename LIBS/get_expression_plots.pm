#DEA1.0 - Differential  Expression analysis using Annocript
# Author: Francesco Musacchia (2015)
#
#This script uses a file with a differential expression analysis and the
#output from Annocript to build a regulation plot which shows most differentially
#expressed GO terms. They will be shown sorted by the highest number of 
#differentially expressed transcripts. In brackets you can read the total
# number of transcripts (DE and NDE) associated with that GO term. 
# NOTICE: With NDE are denoted the transcripts which are evaluated as not
#    			differentially expressed even if we cannot say that they are totally
# 				instead we use a threshold to determine if the fold change and
#					the false discovery rate are enough to be reputed as evidently DE.

#HOW IT WORKS
#1) Fills an hash with all the differential expression of the transcripts 
		#using the input file
#2) For each row of the Annocript output, takes the transcript expression
		#and increases a counter in the hash for GO classesBP, MF and CC for 
		#the associated upregulated, downregulated or not differentially
		# expressed genes

package LIBS::get_expression_plots;
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( execute_expression_plots );


use strict;
use warnings;
use Data::Dumper;
use USEFUL::utilities qw(deallocate_hash);


#Hashes to be used to store the informations belonging to 
#the GO terms
my $hash1;
my $hash2;
my $hash3;
#Where should I find the information? Which columns?
my $trIdIndex;
my $idIndex1;
my $idIndex2;
my $idIndex3;

#Indexes in the table of differentially expressed
my $trNameIndex = 0;
my $fcIndex = 1;
my $fdrIndex = 4;
#What is used to separate the GO terms in the Annocript output?
my $separator;

#I will use the same parameters always used to decide if 
#a transcript is differentially expressed or not
my $minFDR;
my $logFCThr;

#Top to show in the table
my $topN;

#Max length of the GO descriptions
my $maxLengthDescs ;

#Hash with all the differential expression info
my $de_hash;

sub execute_expression_plots{
	my $inputFolder = shift;#Folder where there are some input files
	#Annocript filtered output table
	my $annOut = shift;
	#Complete table with the t-test executed. This script will decide which genes
	#are differentially expressed
	my $de_table = shift;#Table with differentially expressed genes
	my $configHash = shift;
	my $outputFolder = shift;
	my $indexes = shift;#hash of indexes with which to access the Annocript output
	my $suffixes = shift;#suffixes to use for the filenames
	my $dataType = shift; #pathways or go terms
	
	$annOut = $inputFolder."/".$annOut;
	#Other globals
	$separator = $configHash->{'separator'};
	$topN = $configHash->{'topN_exp_plots'};
	$maxLengthDescs = $configHash->{'maxLengthDescs_exp_plots'};
	$minFDR = $configHash->{'minFDR'};
	$logFCThr = $configHash->{'logFCThr'};
	
	#MAIN EXECUTION
	#Uses an hash to store all the information concerning the expression
	print "Storing expression values for each transcript from file $de_table inside an hash...\n";
	deInfo2Hash($outputFolder."/".$de_table);
	#Function which fills three hashes for the GO classes with amount of
	#transcripts belonging to each term
	
	#Define the indexes of the Pathways ids in the Annocript output
	$idIndex1 = $indexes->{'index1'};
	$idIndex2 = $indexes->{'index2'}; 
	$idIndex3 = $indexes->{'index3'};
	if ($dataType  eq 'pathways'){
		$trIdIndex = $indexes->{'tr_id_index'};
		getPathwaysExpressionAmounts($annOut);
	}else{
		$trIdIndex = $indexes->{'tr_id_index'};
		getExpressionAmounts($annOut);
	}

	#Prints sorted tables to be used in a final stacked plot
	print "Printing tables sorted..\n";
	my $classTable = $outputFolder."/".$de_table."_".$suffixes->{'index1'};
	printTableSorted($classTable,"count",$hash1);
	$classTable = $outputFolder."/".$de_table."_".$suffixes->{'index2'};
	printTableSorted($classTable,"count",$hash2);
	$classTable = $outputFolder."/".$de_table."_".$suffixes->{'index3'};
	printTableSorted($classTable,"count",$hash3);
	
	deallocate_hash(\$hash1);#Deallocate hash
	deallocate_hash(\$hash2);#Deallocate hash
	deallocate_hash(\$hash3);#Deallocate hash
}

#Filling the hash with all the differential expression of the transcripts 
sub deInfo2Hash {
	my $de_table = shift;
	
	print "Filling the hash with all the differential expression of the transcripts \n";
	open (DE,"<$de_table");
	#my $de_hash;
	while (my $row_de = <DE>){
			next if $row_de =~ /logFC/;#Jumps the header
			chomp ($row_de);
			my @parts = split("\t",$row_de);	
			$de_hash->{$parts[$trNameIndex]}->{'fc'}= $parts[$fcIndex];#Get the fold change
			$de_hash->{$parts[$trNameIndex]}->{'fdr'}= $parts[$fdrIndex];#Get the FDR
	}
	#print Dumper $de_hash;
close(DE);
}

#Function to compute the log in base 2
sub log2 {
my $n = shift;
return log($n)/log(2);
}

sub getExpressionAmounts {
	my $annOut = shift;
	
	print "Counting the differential expressed for each GO term in each row ".
				" of the Annocript output $annOut\n";
  #print Dumper $de_hash;
	#For each row of the Annocript output, takes the transcript expression and increases
	#a counter in the hash for GO classesBP, MF and CC for the associated upregulated, downregulated
	#or not differentially expressed genes
	open(ANN, "<$annOut") or die "Cannot open $annOut";
	while (my $row = <ANN>){
		chomp($row);
		next if $row =~ /TranscriptName/;
		my @parts = split("\t",$row);
		my $deType = '';#sring for up,down or not changing regulation
			
		#Define the change or not in expression of this transcript
		#To be defined as upregulated or downregulated the FDR should be 
		#very low and the Fold Change will decide if it is up-down
		#otherwise it is NDE
		my $transcript = $parts[0];
		
		#To be up or downregulated the transcript should have
		#a given minimum False Discovery Rate and
		#depending by the fold change sign, if it is greater or lower than a given
		#threshold it is up or down regulated
		if ($de_hash->{$transcript}->{'fdr'} <= $minFDR){
			if ( $de_hash->{$transcript}->{'fc'} >= log2($logFCThr) ){
				$deType = 'up';
			}elsif ($de_hash->{$transcript}->{'fc'} <= -(log2($logFCThr))){
				$deType = 'down';
			 }else{
				 $deType = 'nde';
				} 
		}else{
			#print "Analyzing: $transcript\n";
			#print "FDR: ".$de_hash->{$transcript}->{'fdr'}."\tFC: ".$de_hash->{$transcript}->{'fc'}."\n";
			$deType = 'nde';
		}
		
		
		# Use these transcripts with the hash of the differentially expressed 
		# and get the fold change. Depending by the fold change, sum 1 in the hash 
		# of the GO term assigned to upregulated, down regulated or NDE (if it is 0) 
		# (here we should see if this value is very little and retain zero.  
		#print "GOs: ".$parts[$idIndex1]."\n";
		if ( $parts[$idIndex1] ne '-'){
			my @goTermsBP = split (/\Q$separator\E/,$parts[$idIndex1]);
			my @goDescs = split (/\Q$separator\E/,$parts[$idIndex1+1]);
			my $GOpos = 0;
			foreach my $goTerm (@goTermsBP){
				#This variable is needed to count the number of transcripts which
				#have this GO term but it will not keep in count the 
				#transcripts which are not differentially expressed (NDE)
				#If them all are NDE, this count will be zero
				if ( not defined ($hash1->{$goTerm}->{'count'}) ){
					if ($deType ne 'nde'){ 
						$hash1->{$goTerm}->{'count'} = 1;
					}else{
						$hash1->{$goTerm}->{'count'} = 0;
					}			
				}else{
					#To count only the total of differentially expressed genes
					$hash1->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
				}
				#reducing the description to be used in a plot in R
				if (length($goDescs[$GOpos]) > $maxLengthDescs+10){
					my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs-10;
					$hash1->{$goTerm}->{'desc'} = $reducedDesc." ($goTerm)";
				}else{
					my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs;
					$hash1->{$goTerm}->{'desc'} = $reducedDesc;
				}   
				#Sum this transcript to the count for this GO term
				if (not defined $hash1->{$goTerm}->{$deType}){ 
					$hash1->{$goTerm}->{$deType} = 1;
				}else{
					$hash1->{$goTerm}->{$deType}++;
				}
				$GOpos++;
			}
		}
		if ( $parts[$idIndex2] ne '-'){
			my @goTermsMF = split (/\Q$separator\E/,$parts[$idIndex2]);
			my @goDescs = split (/\Q$separator\E/,$parts[$idIndex2+1]);
			my $GOpos = 0;
			foreach my $goTerm (@goTermsMF){
				#Counting only differentially expressed terms
				if ( not defined $hash2->{$goTerm}->{'count'} ){
					if ($deType ne 'nde'){ 
						$hash2->{$goTerm}->{'count'} = 1;
					}else{
						$hash2->{$goTerm}->{'count'} = 0;
					}
				}else{
					$hash2->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
				}
				#reducing the description to be used in a plot in R
				if (length($goDescs[$GOpos]) > $maxLengthDescs+10){
					my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs-10;
					$hash2->{$goTerm}->{'desc'} = $reducedDesc." ($goTerm)";
				}else{
					my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs;
					$hash2->{$goTerm}->{'desc'} = $reducedDesc;
				} 
				#Sum this transcript to the count for this GO term
				if (not defined $hash2->{$goTerm}->{$deType}){ 
					$hash2->{$goTerm}->{$deType} = 1;
				}else{
					$hash2->{$goTerm}->{$deType}++;
				}
				$GOpos++;
			}
		}
		if ( $parts[$idIndex3] ne '-'){
			#Splits the GO terms and descriptions using the separator
			my @goTermsCC = split (/\Q$separator\E/,$parts[$idIndex3]);
			my @goDescs = split (/\Q$separator\E/,$parts[$idIndex3+1]);
			
			#This is used to determine which is the position to fetch the description		
			my $GOpos = 0;
			#For each GO term sum this transcript to their count of transcripts
			#associated
			foreach my $goTerm (@goTermsCC){
				if ( not defined $hash3->{$goTerm}->{'count'}){
					if ($deType ne 'nde'){ 
						$hash3->{$goTerm}->{'count'} = 1;
					}else{
						$hash3->{$goTerm}->{'count'} = 0;
					}
				}else{
					$hash3->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
				}
				#reducing the description to be used in a plot in R
				if (length($goDescs[$GOpos]) > $maxLengthDescs+10){
					my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs-10;
					$hash3->{$goTerm}->{'desc'} = $reducedDesc." ($goTerm)";
				}else{
					my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs;
					$hash3->{$goTerm}->{'desc'} = $reducedDesc;
			 } 
				#Sum this transcript to the count for this GO term
				if (not defined $hash3->{$goTerm}->{$deType}){ 
					$hash3->{$goTerm}->{$deType} = 1;
				}else{
					$hash3->{$goTerm}->{$deType}++;
				}
				$GOpos++;
			}
		}
	} 
	close(ANN);
}

sub getPathwaysExpressionAmounts {
	my $annOut = shift;
	
	print "Counting the differential expressed for each Pathways in each row ".
				" of the Annocript output $annOut\n";
  #print Dumper $de_hash;
	#For each row of the Annocript output, takes the transcript expression and increases
	#a counter in the hash for Pathways Levels 1,2 and 3 for the associated upregulated, downregulated
	#or not differentially expressed genes
	open(ANN, "<$annOut") or die "Cannot open $annOut";
	while (my $row = <ANN>){
		chomp($row);
		next if $row =~ /TranscriptName/;
		my @parts = split("\t",$row);
		my $deType = '';#sring for up,down or not changing regulation
			
		#Define the change or not in expression of this transcript
		#To be defined as upregulated or downregulated the FDR should be 
		#very low and the Fold Change will decide if it is up-down
		#otherwise it is NDE
		my $transcript = $parts[0];
		
		#To be up or downregulated the transcript should have
		#a given minimum False Discovery Ratea and
		#depending by the fold change sign, if it is greater or lower of a given
		#threshold it is up or down regulated
		if ($de_hash->{$transcript}->{'fdr'} <= $minFDR){
			if ( $de_hash->{$transcript}->{'fc'} >= log2($logFCThr) ){
				$deType = 'up';
			}elsif ($de_hash->{$transcript}->{'fc'} <= -(log2($logFCThr))){
				$deType = 'down';
			 }else{
				 $deType = 'nde';
				} 
		}else{
			#print "Analyzing: $transcript\n";
			#print "FDR: ".$de_hash->{$transcript}->{'fdr'}."\tFC: ".$de_hash->{$transcript}->{'fc'}."\n";
			$deType = 'nde';
		}
		
		
		# Use these transcripts with the hash of the differentially expressed 
		# and get the fold change. Depending by the fold change, sum 1 in the hash 
		# of the GO term assigned to upregulated, down regulated or NDE (if it is 0) 
		# (here we should see if this value is very little and retain zero.  
		#print "GOs: ".$parts[$idIndex1]."\n";
		if ( $parts[$idIndex1] ne '-'){
			my @goTermsBP = split (/\Q$separator\E/,$parts[$idIndex1]);
			
			foreach my $goTerm (@goTermsBP){
				#This variable is needed to count the number of transcripts which
				#have this Pathway but it will not keep in count the 
				#transcripts which are not differentially expressed (NDE)
				#If them all are NDE, this count will be zero
				if ( not defined ($hash1->{$goTerm}->{'count'}) ){
					if ($deType ne 'nde'){ 
						$hash1->{$goTerm}->{'count'} = 1;
					}else{
						$hash1->{$goTerm}->{'count'} = 0;
					}			
				}else{
					#To count only the total of differentially expressed genes
					$hash1->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
				}
				#reducing the description to be used in a plot in R
				if (length($goTerm) > $maxLengthDescs+10){
					my $reducedDesc = substr $goTerm, 0, $maxLengthDescs-10;
					$hash1->{$goTerm}->{'desc'} = $reducedDesc;
				}else{
					my $reducedDesc = substr $goTerm, 0, $maxLengthDescs;
					$hash1->{$goTerm}->{'desc'} = $reducedDesc;
				}   
				#Sum this transcript to the count for this GO term
				if (not defined $hash1->{$goTerm}->{$deType}){ 
					$hash1->{$goTerm}->{$deType} = 1;
				}else{
					$hash1->{$goTerm}->{$deType}++;
				}
			}
		}
		if ( $parts[$idIndex2] ne '-'){
			my @goTermsMF = split (/\Q$separator\E/,$parts[$idIndex2]);

			foreach my $goTerm (@goTermsMF){
				#Counting only differentially expressed terms
				if ( not defined $hash2->{$goTerm}->{'count'} ){
					if ($deType ne 'nde'){ 
						$hash2->{$goTerm}->{'count'} = 1;
					}else{
						$hash2->{$goTerm}->{'count'} = 0;
					}
				}else{
					$hash2->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
				}
				#reducing the description to be used in a plot in R
				if (length($goTerm) > $maxLengthDescs+10){
					my $reducedDesc = substr $goTerm, 0, $maxLengthDescs-10;
					$hash2->{$goTerm}->{'desc'} = $reducedDesc;
				}else{
					my $reducedDesc = substr $goTerm, 0, $maxLengthDescs;
					$hash2->{$goTerm}->{'desc'} = $reducedDesc;
				} 
				#Sum this transcript to the count for this GO term
				if (not defined $hash2->{$goTerm}->{$deType}){ 
					$hash2->{$goTerm}->{$deType} = 1;
				}else{
					$hash2->{$goTerm}->{$deType}++;
				}
			}
		}
		if ( $parts[$idIndex3] ne '-'){
			#Splits the GO terms and descriptions using the separator
			my @goTermsCC = split (/\Q$separator\E/,$parts[$idIndex3]);

			#For each GO term sum this transcript to their count of transcripts
			#associated
			foreach my $goTerm (@goTermsCC){
				if ( not defined $hash3->{$goTerm}->{'count'}){
					if ($deType ne 'nde'){ 
						$hash3->{$goTerm}->{'count'} = 1;
					}else{
						$hash3->{$goTerm}->{'count'} = 0;
					}
				}else{
					$hash3->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
				}
				#reducing the description to be used in a plot in R
				if (length($goTerm) > $maxLengthDescs+10){
					my $reducedDesc = substr $goTerm, 0, $maxLengthDescs-10;
					$hash3->{$goTerm}->{'desc'} = $reducedDesc;
				}else{
					my $reducedDesc = substr $goTerm, 0, $maxLengthDescs;
					$hash3->{$goTerm}->{'desc'} = $reducedDesc;
			 } 
				#Sum this transcript to the count for this GO term
				if (not defined $hash3->{$goTerm}->{$deType}){ 
					$hash3->{$goTerm}->{$deType} = 1;
				}else{
					$hash3->{$goTerm}->{$deType}++;
				}
			}
		}
	} 
	close(ANN);
}

#Prints a table with for each row a GO term and three columns represent
#the number of upregulated, downregulated, not differentially expressed
#transcripts associated to such term.
sub printTableSorted {
	my $expTable = shift;
	my $field = shift;
	my $hash = shift;

	open (EXP_T,">$expTable") or die "Unable to open $expTable";
	my $numKeys = 0;
	
	 # The hash is sorted putting before the GO terms which have more
	 # differentially expressed transcripts: the 'count' is given as the 
	 # sum of both the up and down regulated.
	 # GO terms build a final table with 
	 # each column being a GO term and three rows represent the 
	 # number of UpRegulated, DownRegulated and NDE genes
	 foreach my $key ( #
		 sort { $hash->{$b}->{$field} <=> $hash->{$a}->{$field} } #
		 keys %{$hash}
	 )
	 {
			#These variables could be not initialized if the particular differential expression
			#is not found. Thus the following is initialization
			if (not defined $hash->{$key}->{'up'}){ 
				$hash->{$key}->{'up'} = 0;
			}
			if (not defined $hash->{$key}->{'down'}){ 
				$hash->{$key}->{'down'} = 0;
			}
			if (not defined $hash->{$key}->{'nde'}){ 
				$hash->{$key}->{'nde'} = 0;
			}
			
			my $totCountDE = $hash->{$key}->{'count'};#Count of only the differentially expressed
			my $overallCount = $totCountDE + $hash->{$key}->{'nde'};#Overall count with also NDE
			#print "GO: $key Overall count: $overallCount - Total count diff expr: $totCountDE \n";
			#print $hash->{$key}->{'up'}." ".$hash->{$key}->{'down'}." ".$hash->{$key}->{'nde'}."\n";
			my $upPerc = 0;
			my $downPerc = 0;
			my $ndePerc = 0;
			
			#To avoid the division by zero I used this commented check
			#But this should never be zero
			#if ( $overallCount != 0){#DEBUGCODE
			 $upPerc = ($hash->{$key}->{'up'}/$overallCount)*100;
			 $downPerc = ($hash->{$key}->{'down'}/$overallCount)*100;
			 $ndePerc = ($hash->{$key}->{'nde'}/$overallCount)*100;
			#}
			
			my $desc = $hash->{$key}->{'desc'}." [$overallCount]";
			#$desc =~ s/ /_/g;
			
			print EXP_T join("\t",$desc,$upPerc,$ndePerc,$downPerc);
			#print EXP_T join("\t",$key,$upPerc,$ndePerc,$downPerc);
			print EXP_T "\n";
			$numKeys ++;
			last unless $numKeys < $topN;
	 } 
	close (EXP_T); 
}

1;
