#!/usr/bin/perl
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

#=======================================================================================================================================================
#					PE_illumina_assembly_pipe.pl
#
#   1:  If you want to assemble more than one genome it is easier to use the PE_illumina_assembly_prepipe.pl script to create automated assembly of them all.
#   2:  Start this script with: "perl PE_illumina_assembly_pipe.pl fastq1 fastq2 trim/notrim 40 (wanted input coverage X) 1641481 (genome size Mbp) pilon/nopilon ariba/noariba kraken/nokraken

# All the parameters have to be given every time and in the same order as above, otherwise the script will return error message.


# IMPORTANT! The script extracts a short filename for the assembly from the fastq-file. It starts from beginning and stops at first underscore "_" so if your files all start with 1_ for instance, then all the assemblies will be named "1" and overwrite each other to the end of the pipe. So chose a unique filename before the first underscore. This is done automatically if using the Miseq which puts the sample ID first on FASTQs. But the Hiseq result files all start the same so you have to rename them.

# IMPORTANT! If you put this script on new computer, update the script (under arguments) with how many CPU-threads and RAM-size available. Also update path to Kraken-database, Pilon and trimmomatic if you want to run them.

# IMPORTANT! If no assembly is needed and for instance AMR-genes and/or species determination is the only goal, set the wanted coverage to "0". No SPAdes assembly will then be performed.



# This pipe does the following:

# 1: Extract filename from FASTQ-file.
# 2: Start a log-file where one can see what the pipe did and also find metrics for assembly.
# 3: Test if FASTQs are zipped, unzip if they are.
# 4: Run Trimmomatic 0.36 and then zip original files to save space (OPTIONAL).
# 5: Download CARD-reference AMR database and run ARIBA on it (OPTIONAL). 
# 6: Download ResFinder-reference AMR database and run ARIBA on it (OPTIONAL).
# 7: Run Kraken on reads. Results are extracted, sorted and printed in log-file.(OPTIONAL).
# 8: Gzip original reads-files
# 9: Count bases in FASTQs and extract as many bases needed to achieve wanted input coverage (If wanted_coverage > 0). Gzip original files to save space (if trim = no).
# 6: Start SPAdes with "-careful" parameter, change in script if you don't want careful assembly. (OPTIONAL).
# 7: Move extracted reads into assembly-folder that is named "[filename-short]_[wanted-coverage]_spades".
# 8: Pilon correction, starts with: Bowtie2 alignment, SAM-> BAM-> Pilon corrected assembly (OPTIONAL).
# 9: Analyse original SPAdes-assembly and print metrics in log (If wanted_coverage > 0).
# 10: Analyse Pilon-corrected assembly and print metrics in log (OPTIONAL).
# 11: Complete ARIBA results and move to assembly-folder
# 12:  Cleanup, removes BAMs, corrected reads from SPAdes, trimmomatic-outputs etc. Change manually in bottom of script if you want to keep some of these files.


# "trim" uses Trimmomatic. Pilon corrects assemblies. Ariba searches the reads for AMR-genes and mutations using CARD and ResFinder. Kraken classifies reads   
# according to taxon which can be used to quality check the reads, for species-identification and for metagenomics. If using ARIBA, also use Kraken to make sure AMR-genes are from correct genome.



# Latest modification to script: 190130  (mod for Robert that puts all assemblies in one folder) there is newer version on 72-threds
# //Joakim Skarin


#=================================================================================================================================================
#		Arguments and tests - update here when using script on new computer
#=================================================================================================================================================

$fastqfile1 = $ARGV[0];
$fastqfile2 = $ARGV[1];
$trim = $ARGV[2];
$wanted_coverage = $ARGV[3];
$genome_size = $ARGV[4];
$run_pilon = $ARGV[5];
$run_ariba = $ARGV[6];
$run_kraken = $ARGV[7];


$skip_assembly = 0;

if ($wanted_coverage eq 0) {

	$skip_assembly = 1;

}


if ($run_kraken eq "kraken") {

	$kraken_DB = "/home/asgeir/software/minikraken_8GB_20200312";	# UPDATE HERE
	
		if (-e $kraken_DB) {}

		else { die("Fatal error: The Kraken database was not found. Please update path in perl script"); }

}

$threads_available = "42";		# UPDATE HERE. Set this to number of simultaneaus threads CPU can handle, usually twice the number of CPU-cores
$RAM_available = "224";			# UPDATE HERE. Set this to max amount of GB RAM-memory you want to use, usually a couple of GB below total number of GB RAM


if ($trim eq "trim") {

	$illum_adaptors_path = "/home/asgeir/miniconda3/envs/campy/share/trimmomatic/adapters/NexteraPE-PE.fa";	# UPDATE HERE

	if (-e $illum_adaptors_path) {}
	else { die("Fatal error: trimmomatic illumina adaptor fasta file was not found. Please update path in perl script"); }
	
}

$my_dir = getcwd;

if (-e "$my_dir/all_assemblies") { }

else  { 
	mkdir ("$my_dir/all_assemblies") ;
}


#==========================================================================================
#			sub timestamp
#==========================================================================================

sub getLoggingTime {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}

#==========================================================================================
#			EXTRACT FILENAME FROM READS-FILE
#==========================================================================================


($filename1_short) = $fastqfile1 =~ /^(.*?)_/;


$filename2_short = $filename1_short;



#===========================================================================================
#			START A LOG
#===========================================================================================



open(HANDLE_LOG,">"."$filename1_short\_log_and_metrics.log") || die("FATAL ERROR: could not open file:  log-file ");

$timestamp = getLoggingTime();

print HANDLE_LOG "Assembly-pipe started at at $timestamp in folder $my_dir\n\n";

print HANDLE_LOG "\n==============================================================================================================\n\n";
print HANDLE_LOG "Arguments:      Fastq file 1:\t\t $fastqfile1\n";
print HANDLE_LOG "	        Fastq file 2:\t\t $fastqfile2\n";
print HANDLE_LOG "	        Trim using trimmomatic:\t $trim\n";
print HANDLE_LOG "	        Wanted coverage:\t $wanted_coverage\n";
print HANDLE_LOG "		Genome Size:\t\t $genome_size\n";
print HANDLE_LOG "		Run Pilon:\t\t $run_pilon\n";
print HANDLE_LOG "		Run ARIBA:\t\t $run_ariba\n";
print HANDLE_LOG "		Run Kraken:\t\t $run_kraken\n\n\n\n";

print HANDLE_LOG "Max threads:\t $threads_available\n\n";
print HANDLE_LOG "Max RAM used:\t $RAM_available GB\n\n";


if ($skip_assembly eq 1) {

	print HANDLE_LOG "Wanted coverage set to 0, skipping assembly (and Pilon if chosen)\n";
	print HANDLE_LOG "The reads will still be extracted to enable Kraken and ARIBA\n\n";
}


if ($run_kraken eq "kraken") {

	print HANDLE_LOG "Kraken database path: $kraken_DB\n\n";		

}

print HANDLE_LOG "\n==============================================================================================================\n\n";

#============================================================================================
#		TEST IF READS ARE GZIPPED
#============================================================================================


$fastqfile1_zipped = 0;
$fastqfile2_zipped = 0;


if($fastqfile1 =~ /gz/ ) { 
	$fastqfile1_zipped = 1;
	system ("gunzip $fastqfile1"); 
	print HANDLE_LOG  "Fastq file 1 was unzipped using gzip\n\n";
 	($fastqfile1) = $fastqfile1 =~ /^(.*?).gz/;	# Must change name to without .gz

}

else { 
	$fastqfile1_zipped = 0;
	 
	print HANDLE_LOG  "Fastq file 1 was already unzipped\n\n";
}



if($fastqfile2 =~ /gz/ ) { 
	$fastqfile2_zipped = 1;
	system ("gunzip $fastqfile2"); 
	print HANDLE_LOG  "Fastq file 2 was unzipped using gzip\n\n";
 	($fastqfile2) = $fastqfile2 =~ /^(.*?).gz/;	# Must change name to without .gz

}

else { 
	$fastqfile2_zipped = 0;
 
	print HANDLE_LOG  "Fastq file 2 was already unzipped\n\n";

}

#=======================================================================================================
#			Trimmomatic
#=======================================================================================================

if ($trim eq "trim") {

	print HANDLE_LOG  "Trimming reads...";


	system ("trimmomatic PE -threads $threads_available -phred33 $fastqfile1 $fastqfile2 $fastqfile1.trimmed $fastqfile1.unpaired_trimmed $fastqfile2.trimmed $fastqfile2.unpaired_trimmed ILLUMINACLIP:$illum_adaptors_path:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");


	print HANDLE_LOG  "done\n";

}

#ARIBA and KRAKEN moved up here to circumvent running on trimmed and/or subset extracted reads since AMR-findings differed when using only read-subset
# ARIBA results are moved to assembly folder later and results are written in log at the same time

#========================================================================================================
#				ARIBA with CARD-database
#========================================================================================================


if($run_ariba eq "ariba")  {

	print HANDLE_LOG "\n==============================================================================================================\n";

	chdir ("$my_dir");	

	$timestamp = getLoggingTime();

	print HANDLE_LOG "\n\nARIBA (CARD) started at $timestamp ";

	$test_CARD = "CARD_reference_dataset_downloaded";

	if(-e $test_CARD) {

		print HANDLE_LOG "\n\nCARD reference dataset present, will not download new dataset";
				
         }
	
	else { 

		print HANDLE_LOG "\n\nNo CARD reference dataset present";
  		
		print HANDLE_LOG "\n\nDownloading CARD reference dataset";

		system ("ariba getref card out.card");

		print HANDLE_LOG "\n\nPreparing dataset for analysis";

		system ("ariba prepareref -f out.card.fa -m out.card.tsv cardref.out");

		print HANDLE_LOG "\n\nCreating file 'CARD_reference_dataset_downloaded'";		

				open(HANDLE_CARD_REF,">"."$test_CARD") || die("FATAL ERROR: could not create file:  CARD_reference_dataset_downloaded "); # Create empty file when database has been downloaded so that download only happens for the first genome


	}

	print HANDLE_LOG "\n\nRunning ARIBA (CARD database)...";
	
	system ("ariba run cardref.out $fastqfile1 $fastqfile2 $filename1_short\_ARIBA\_CARD\_results/");  



#========================================================================================================
#			ARIBA with ResFinder database
#========================================================================================================

print HANDLE_LOG "\n==============================================================================================================\n";


chdir ("$my_dir");	

	$timestamp = getLoggingTime();

	print HANDLE_LOG "\n\nARIBA (ResFinder) started at $timestamp ";

	$test_ResFinder = "ResFinder_reference_dataset_downloaded";

	if(-e $test_ResFinder) {

		print HANDLE_LOG "\n\nResFinder reference dataset present, will not download new dataset";
				
         }
	
	else { 

		print HANDLE_LOG "\n\nNo ResFinder reference dataset present";
  		
		print HANDLE_LOG "\n\nDownloading ResFinder reference dataset";

		system ("ariba getref resfinder out.resfinder");

		print HANDLE_LOG "\n\nPreparing dataset for analysis";

		system ("ariba prepareref -f out.resfinder.fa -m out.resfinder.tsv resfinderref.out");

		print HANDLE_LOG "\n\nCreating file 'ResFinder_reference_dataset_downloaded'";		

		open(HANDLE_RESFINDER_REF,">"."$test_ResFinder") || die("FATAL ERROR: could not create file:  ResFinder_reference_dataset_downloaded "); # Create empty file when database has been downloaded so that download only happens for the first genome


	}

	print HANDLE_LOG "\n\nRunning ARIBA (ResFinder database)...";
	
	system ("ariba run resfinderref.out $fastqfile1 $fastqfile2 $filename1_short\_ARIBA\_ResFinder\_results/");

}



#=======================================================================================================
#			  Kraken
#=======================================================================================================


if ($run_kraken eq "kraken") {

	#chdir "$path_assembly";

	$timestamp = getLoggingTime();

		
	if (-e $kraken_DB )  {
	
		print HANDLE_LOG "\n\nKraken started at $timestamp with database $kraken_DB and $fastqfile1 and $fastqfile2 as input\n\n";

		system ("kraken2 --db $kraken_DB --threads $threads_available --output $filename1_short.kraken.out --report $filename1_short.kraken.report --paired $fastqfile1 $fastqfile2");

		$timestamp = getLoggingTime();

		open(KRAKEN_REPORT,"$filename1_short.kraken.report") || die("FATAL ERROR: could not open file: $filename1_short.kraken.report");
	
		@kraken1 = ();		#Extract data from columns 1, 3 and 6 from .report-file
		@kraken3 = ();
		@kraken6 = ();
		@index = ();


		while(<KRAKEN_REPORT>)  {

			@columns = split (/\t/);
			push @kraken1, $columns[0];
			push @kraken3, $columns[2];
			push @kraken6, $columns[5];

			@columns = ();


		}		

		close KRAKEN_REPORT;

		#sort array1 containing "Percentage of reads covered by the clade". Then apply sorting on array 3 and 6 to

		@index = sort { $kraken1[$b] <=> $kraken1[$a] } 0 .. $#kraken1;


		@kraken1 = @kraken1[@index];
		@kraken3 = @kraken3[@index];
		@kraken6 = @kraken6[@index];


		s{^\s+|\s+$}{}g foreach @kraken1;
		s{^\s+|\s+$}{}g foreach @kraken3;


		$max_kraken_results_showing = 40;   #Raise this number if you want to show more lines of result from Kraken


		print HANDLE_LOG "\n==============================================================================================================\n";
		print HANDLE_LOG "\t Results of Kraken analysis of reads, sorted by % of reads covered by clade, showing max $max_kraken_results_showing rows\n";
		print HANDLE_LOG "==============================================================================================================\n\n";

		print HANDLE_LOG "% of reads covered by clade         Reads assigned directly to this taxon         Name\n";
		print HANDLE_LOG "--------------------------------------------------------------------------------------------------------------\n";

		$i = 0;

		foreach (@kraken1) {

			print HANDLE_LOG "$kraken1[$i]                               $kraken3[$i]                                            $kraken6[$i]";
			$i++;

			if ($i eq $max_kraken_results_showing) { last; }		

		}

		print HANDLE_LOG "--------------------------------------------------------------------------------------------------------------\n\n";


	}

	else { print "\n\nCould not find krakenDB, check path\nSkipping kraken-analysis" }


}
	

	#======  Gzip original reads-files if Trim was used, otherwise they are zipped after read-extraction =======

if ($trim eq "trim") {
	if ($fastqfile1_zipped eq 1) {

		system ("gzip $fastqfile1");
		print HANDLE_LOG "\n\n$fastqfile1 was zipped using gzip\n";
	}

	if ($fastqfile2_zipped eq 1) {
	
		system ("gzip $fastqfile2");
		print HANDLE_LOG "\n$fastqfile2 was zipped using gzip\n";
	}



	#================================= CHANGE NAMES OF READS-FILES ========================================

	$fastqfile1 = "$fastqfile1.trimmed";	#Now, when fastqfile1/2 are zipped, point at ".trimmed-file to continue with assembly using trimmed reads"
	$fastqfile2 = "$fastqfile2.trimmed";
}




#====================================================================================================
#			NUMBER OF READS NEEDED TO MATCH WANTED COVERAGE
#====================================================================================================


if ($skip_assembly eq 0) {


	$bases_needed = ($wanted_coverage * $genome_size) / 2; #divide by 2 due to PE



	$bases_needed = int($bases_needed);  					# rounding down

	print HANDLE_LOG  "\nTo achieve $wanted_coverage X, $bases_needed bases are needed from each fastq-file\n"; 

	open(FASTQ1_FILE,$fastqfile1) || die("FATAL ERROR: could not open file: $fastqfile1 ");


	$total_bases = 0;
	$read_counter = 0;
	$radtyp = 0;

 	 while ( <FASTQ1_FILE> )
	 {

 	   $row = $_;

 	   $radtyp++;

 	   if($radtyp == 2)     { 
         
  	     chomp $row;
	       $total_bases += length($row);

 	   }

  	  if($radtyp >= 4){
      
	
		$read_counter++;
	
		$radtyp = 0;
	
		if ($total_bases >= $bases_needed) {
	
			last;
 		}

 	   }

 	 }

	close FASTQ1_FILE;  #Error if not closed here and opened later


	if ($total_bases > $bases_needed) {

		print HANDLE_LOG "\n\nNeeded bases: $bases_needed (per direction) which amounts to $read_counter reads from fastq_1 which is $total_bases bases\n\n";

	}

	if($total_bases < $bases_needed) {

		print HANDLE_LOG "\n\nThere are not enough bases to achieve $wanted_coverage X coverage.\n";

		$wanted_coverage = ($total_bases / $genome_size) * 2; #a new coverage is calculated

		$wanted_coverage = int($wanted_coverage);
  
		print HANDLE_LOG "\nUsing an estimated coverage of $wanted_coverage X instead which amounts to $read_counter reads from fastq_1 which is $total_bases bases\n";
	}

$reads_needed = $read_counter;



#==============================================================================================
#			EXTRACT READS FROM FASTQ
#==============================================================================================


	open(FASTQ1_FILE,$fastqfile1) || die("FATAL ERROR: could not open file: $fastqfile1 ");
	open(FASTQ2_FILE,$fastqfile2) || die("FATAL ERROR: could not open file: $fastqfile2 ");


	$counter = 0;
	while (<FASTQ1_FILE>) { $counter += 1; }		#count number of reads in fastq1
	$reads_in_fastq1 = $counter / 4;
	print HANDLE_LOG "\nCounting reads in fastq-file 1 ... $reads_in_fastq1\n";

	$counter = 0;
	while (<FASTQ2_FILE>) { $counter += 1; }		#count number of reads in fastq2
	$reads_in_fastq2 = $counter / 4;
	print HANDLE_LOG "\n\nCounting reads in fastq-file 2 ... $reads_in_fastq2\n";


	if ($reads_in_fastq1 ne $reads_in_fastq2) { 
		print HANDLE_LOG "\nFATAL ERROR: The number of reads in file 1 did not match file 2\n";
		die "\nThe number of reads in $fastqfile1 did not match $fastqfile2\n";
	}


	$read1_output = "$filename1_short\_R1\_$wanted_coverage\X.fastq";	#create filename for extracted reads

	$read2_output = "$filename2_short\_R2\_$wanted_coverage\X.fastq";



	$lines_to_extract = $reads_needed * 4;

	system ("head -n $lines_to_extract $fastqfile1 > $read1_output");

	print HANDLE_LOG "\n$reads_needed reads ($lines_to_extract lines) were extracted from $fastqfile1 and put in $read1_output\n";

	system ("head -n $lines_to_extract $fastqfile2 > $read2_output");

	print HANDLE_LOG "\n$reads_needed reads ($lines_to_extract lines) were extracted from $fastqfile2 and put in $read2_output\n";

	close FASTQ1_FILE;
	close FASTQ2_FILE;

} #End if ($skip_assembly eq 0)


if ($skip_assembly eq 1) {


	$path_assembly = "$filename1_short\_$wanted_coverage\X_no-spades";

	mkdir ($path_assembly);


}



#========================================================================================================
#			zip original files if trim was not used
#========================================================================================================

if ($trim eq "notrim") {

	if ($fastqfile1_zipped eq 1) {

		system ("gzip $fastqfile1");
		print HANDLE_LOG "\n\n$fastqfile1 was zipped using gzip\n";
	}

	if ($fastqfile2_zipped eq 1) {
	
		system ("gzip $fastqfile2");
		print HANDLE_LOG "\n$fastqfile2 was zipped using gzip\n";
	}
}

print HANDLE_LOG "\n==============================================================================================================\n";

#=======================================================================================================
#			  start Spades
#=======================================================================================================



if ($skip_assembly eq 0) {

	$timestamp = getLoggingTime();

#december_2021_updated for SPAdes 3.15.3

	system ("spades.py --careful -o $filename1_short\_$wanted_coverage\X_spades --pe1-1 $read1_output --pe1-2 $read2_output -t $threads_available -m $RAM_available");

	print HANDLE_LOG "\nSPAdes 3.15.3 was started at $timestamp with parameters: --careful -o $filename1_short\_$wanted_coverage\X_spades --pe1-1 $read1_output --pe1-2 $read2_output -t $threads_available -m $RAM_available\n";


	$path_assembly = "$filename1_short\_$wanted_coverage\X_spades";

	rename ("$path_assembly/contigs.fasta", "$path_assembly/$filename1_short.fasta") or die "\nrename of fasta failed: $!"; #rename must happen here to fit both pilon and nopilon

	$timestamp = getLoggingTime();

	print HANDLE_LOG "\nSPAdes finished at $timestamp\n";

	print HANDLE_LOG "\nThe fasta file in the spades output was renamed from 'contigs.fasta' to '$filename1_short.fasta'\n";


	


#========================================================================================================
#		Move read-files to assembly-folder
#========================================================================================================


move($read1_output, "$path_assembly/$read1_output") or die "The move operation failed: $!";      

move($read2_output, "$path_assembly/$read2_output") or die "The move operation failed: $!";


} #End if ($skip_assembly eq 0) Run SPAdes


#========================================================================================================
#		 Pilon
#========================================================================================================

if ($skip_assembly eq 0) {

if ($run_pilon eq "pilon") {

	chdir "$path_assembly";

	$timestamp = getLoggingTime();


	system ("bowtie2-build -f --threads $threads_available --quiet $filename1_short.fasta $filename1_short");

	print HANDLE_LOG "\n\nBowtie2 started at $timestamp with parameter --very-sensitive-local";

	system ("bowtie2 -x $filename1_short -1 $read1_output -2 $read2_output -S $filename1_short.sam --phred33 --very-sensitive-local --no-unal -p $threads_available");


	system ("samtools view -bh $filename1_short.sam > $filename1_short.bam");

	system ("samtools sort $filename1_short.bam -o $filename1_short.sorted.bam");

	system ("samtools index $filename1_short.sorted.bam");

	$timestamp = getLoggingTime();

	print HANDLE_LOG "\n\nPilon 1.24 started at $timestamp ";

	system ("pilon --genome $filename1_short.fasta --frags $filename1_short.sorted.bam --output $filename1_short.pilon --changes --threads $threads_available");

	$timestamp = getLoggingTime();

	print HANDLE_LOG "\n\nPilon finished at $timestamp ";

	print HANDLE_LOG "\n\nCorrected fasta file created: $filename1_short.pilon.fasta";
	
	chdir ("$my_dir");

	copy ("$path_assembly/$filename1_short.pilon.fasta", "$my_dir/all_assemblies/$filename1_short.pilon.fasta") or die "Copy failed: $!"; #copy fasta from SPAdes output to a folder to collect all assemblies at once. FIX this so that the uncorrected fasta is moved if nopilon. Take from black comp script.
}
}

#========================================================================================================
#				Metrics
#========================================================================================================

if ($skip_assembly eq 0) {

chdir ("$my_dir");

$path_assembly = "$filename1_short\_$wanted_coverage\X_spades";

}

sub metrics {

	$assembly_path = $_[0];

		print HANDLE_LOG "\n\nOpening $assembly_path...";
	
		open(ASSEMBLY_FASTA,"$assembly_path") || die("FATAL ERROR: could not open file: assembly_fasta");

		print HANDLE_LOG "and looking at metrics of assembly\n";
		print HANDLE_LOG "\n\n==========================================================================";
		print HANDLE_LOG "\n\tMetrics for $assembly_path";
		print HANDLE_LOG "\n==========================================================================";


	$number_of_contigs = 0;
	$bases_in_contig = 0;
	$total_consensus = 0;
	@contig_lengths = "";
	$N50 = 0;
	$non_base = 0;
	$number_AT = 0;
	$number_GC = 0;
	@GC = "";
	@sequence = "";
	$a = 0;
	$b = 0;
	@sorted_contig_lengths = "";
	$temp = 0;
	$N50_temp_var = 0;
	$contigs_over_1000 = 0;
	$total_number_bases = 0;
	
	while($line = <ASSEMBLY_FASTA>)  {


	if ($line =~ />/) {	
		
		if ($bases_in_contig ne 0) {

			push (@contig_lengths, $bases_in_contig);
			$total_consensus += $bases_in_contig;
			$bases_in_contig = 0;
					
		}
		$number_of_contigs += 1;		
	}
	
	else  { 
		chomp $line;			
		
		$bases_in_contig += length($line);
		
		@GC = split(//, $line);

		foreach $base (@GC) {
		
			$total_number_bases++;
		
			if ($base eq "A" || $base eq "T")  {	$number_AT++;  }
	
			elsif ($base eq "G" || $base eq "C")  {    $number_GC++;  }
	
			else { $non_base++; }
	
			@GC = "";
		}
		
	}
	}

	push (@contig_lengths, $bases_in_contig);

	close (ASSEMBLY_FASTA);

	print HANDLE_LOG "\n\nThe number of contigs is $number_of_contigs, the total nr of bases: $total_consensus\n";


	# ==============================================================================================
	#			 N50
	# ==============================================================================================

	@sorted_contig_lengths = sort { $b <=> $a } @contig_lengths;	#Sort array decreasing

	print HANDLE_LOG "\nLongest contig: $sorted_contig_lengths[0]\n";

	for $temp (@sorted_contig_lengths) {		
	

		$N50_temp_var += $temp;
		if ($N50_temp_var >= ($total_consensus / 2)) {
		
			$N50 = $temp;
			last;
		}			
	}	

	print HANDLE_LOG "\nN50: $N50\n";



	# =============================================================================================
	#			Contigs longer than 1 kb
	# =============================================================================================


	$contigs_over_1000 = 0;

	for $temp (@contig_lengths) {		
	
		if ($temp >= 1000) {		
			$contigs_over_1000 += 1;	

		}			
	}

	print HANDLE_LOG "\nThe number of contigs above 1 kb in length is $contigs_over_1000\n";


	# ================================================================================================
	#				GC
	# ================================================================================================


	if ($non_base > 0)  {			# function for non-ACTG
	
		print HANDLE_LOG "The GC-content of the sequence is ", sprintf "%.2f", ($number_GC/($number_GC + $number_AT))*100, "% but $non_base characters were not A, T, C or G and were excluded from GC-calculation\n";

	}

	else  {
		print HANDLE_LOG "\nThe GC-content of the sequence is ", sprintf "%.2f", ($number_GC/$total_number_bases)*100, "%\n";

	}

	print HANDLE_LOG "\n==========================================================================\n";

}

#==================== end of metrics-sub ================================================================
#========================================================================================================


#call metrics-sub:

if ($skip_assembly eq 0) {

	$do_metric = metrics("$path_assembly/$filename1_short.fasta");		


	if ($run_pilon eq "pilon") {
	
		$do_metric_pilon = metrics("$path_assembly/$filename1_short.pilon.fasta");
	
	}

}

#======================================================================================================
#		Build Consed-structure and map reads, if chosen 
#=======================================================================================================

#if ($skip_assembly eq 0) {

#if ($prepare_consed eq "consed") {

	

#	mkdir "$my_dir/$path_assembly/consed";
#	mkdir "$my_dir/$path_assembly/consed/chromat_dir";
#	mkdir "$my_dir/$path_assembly/consed/edit_dir";
#	mkdir "$my_dir/$path_assembly/consed/phdball_dir";
#	mkdir "$my_dir/$path_assembly/consed/phd_dir";
#	mkdir "$my_dir/$path_assembly/consed/sff_dir";
#	mkdir "$my_dir/$path_assembly/consed/solexa_dir";

#	copy ("$path_assembly/$filename1_short.fasta", "$path_assembly/consed/edit_dir/$filename1_short.fasta") or die "Copy failed: $!"; #copy fasta from SPAdes to Consed's edit folder


#	symlink ("$path_assembly/$read1_output", "$path_assembly/consed/edit_dir/$read1_output.link");	#link to reads R1

#	rename ("$path_assembly/consed/edit_dir/$read1_output.link", "$path_assembly/consed/edit_dir/$read1_output") or die "\nrename of link failed: $!"; #remove .link ending on files

#	symlink ("$path_assembly/$read1_output", "$path_assembly/consed/solexa_dir/$read1_output.link");

#	rename ("$path_assembly/consed/solexa_dir/$read1_output.link", "$path_assembly/consed/solexa_dir/$read1_output") or die "\nrename of link failed: $!"; 

#	open(HANDLE_FASTQFOF,">"."$path_assembly/consed/edit_dir/fastq.fof") || die("FATAL ERROR: could not open file:  fastq.fof ");

#	print HANDLE_FASTQFOF "$read1_output";

#	close HANDLE_FASTQFOF;

#	chdir ("$my_dir/$path_assembly/consed/edit_dir") || die "Cannot chdir to consed's edit folder ($!)";

#	$timestamp = getLoggingTime();

#	print HANDLE_LOG "\nfasta2Ace and AddSolexaReads started at $timestamp\n\n";

#	system ("fasta2Ace.perl $filename1_short.fasta");

#	system ("addSolexaReads.perl $filename1_short.ace fastq.fof $filename1_short.fasta");
	
	
#} #end if consed
#} #end if skip_assembly


#========================================================================================
#		Complete ARIBA results and move to assembly-folder
#========================================================================================

if ($run_ariba eq "ariba") {

#CARD Funkar inte längre, kolla om ARIBA fixat sökvägar till CARD

#	$path_ARIBA_CARD_results = "$filename1_short\_ARIBA\_CARD\_results";

#	move($path_ARIBA_CARD_results, "$path_assembly/$path_ARIBA_CARD_results") or die "The move operation failed: $!"; 

#	print HANDLE_LOG "\n\nThe complete ARIBA-report can be found in in folder $path_assembly/$filename1_short\_ARIBA\_CARD\_results/report.tsv\n\nA summary of the genes identified with CARD:\n\n";

#	close HANDLE_LOG;


#	system ("cut -f1,4,30,31 $path_assembly/$filename1_short\_ARIBA\_CARD\_results/report.tsv | uniq >> $filename1_short\_log_and_metrics.log");


#	open(HANDLE_LOG,">>"."$filename1_short\_log_and_metrics.log") || die("FATAL ERROR: could not open file:  log-file ");


#ResFinder

	$path_ARIBA_ResFinder_results = "$filename1_short\_ARIBA\_ResFinder\_results";

	move($path_ARIBA_ResFinder_results, "$path_assembly/$path_ARIBA_ResFinder_results") or die "The move operation failed: $!"; 

	print HANDLE_LOG "\n\nThe complete ARIBA-report can be found in in folder $path_assembly/$filename1_short\_ARIBA\_ResFinder\_results/report.tsv\n\nA summary of the genes identified with ResFinder:\n\n";

	close HANDLE_LOG;


	system ("cut -f1,4,30,31 $path_assembly/$filename1_short\_ARIBA\_ResFinder\_results/report.tsv | uniq >> $filename1_short\_log_and_metrics.log");


	open(HANDLE_LOG,">>"."$filename1_short\_log_and_metrics.log") || die("FATAL ERROR: could not open file:  log-file ");

}


#========================================================================================================
#			Cleanup - removes all the large files
#========================================================================================================


# Move Kraken-results to assembly-folder

if ($run_kraken eq "kraken") {

	move("$filename1_short.kraken.report", "$path_assembly/$filename1_short.kraken.report") or die "The move operation failed: $!";
	move("$filename1_short.kraken.out", "$path_assembly/$filename1_short.kraken.out") or die "The move operation failed: $!";
}

$cleanup = 1;		#set to something else than 1 if you want to skip cleanup

if ($cleanup eq 1) {

	print HANDLE_LOG "\n\nCleanup - unlinking files from trimmomatic, corrected reads from SPAdes and sam/bam's\n";	
	
	if ($trim eq "trim")  {
		

		($fastqfile1_temp) = $fastqfile1 =~ /(.*)\./;		#remove ".trimmed" in variable to be able to unlink unpaired etc.

		($fastqfile2_temp) = $fastqfile2 =~ /(.*)\./;

		chdir ("$my_dir");
		unlink "$fastqfile1_temp.trimmed";
		unlink "$fastqfile1_temp.unpaired_trimmed";
		unlink "$fastqfile2_temp.trimmed";
		unlink "$fastqfile2_temp.unpaired_trimmed";

	}

	chdir ("$path_assembly");
	
	unlink "$read1_output";		#skip this with "#" if you want to keep extracted reads
	unlink "$read2_output";		#skip this with "#" if you want to keep extracted reads


	unlink "$filename1_short.kraken.out";

	
	if ($skip_assembly eq 0) {

		if ($run_pilon eq "pilon") {
	
			unlink "$filename1_short.sam";
			unlink "$filename1_short.bam";
			unlink "$filename1_short.sorted.bam";		#keep this and the .bai below if you want to look at assembly alignment in Tablet or similar program
			unlink "$filename1_short.sorted.bam.bai";
		}
	
	
	chdir ("corrected");	
		
	unlink "$filename1_short\_R1_$wanted_coverage\X.00.0_0.cor.fastq.gz";	#removes corrected reads created by SPAdes
	unlink "$filename1_short\_R2_$wanted_coverage\X.00.0_0.cor.fastq.gz";
	unlink "$filename1_short\_R_unpaired.00.0_0.cor.fastq.gz";

	}

} #end if cleanup

#==========================================================================================================
#==========================================================================================================


chdir ("$my_dir");

$timestamp = getLoggingTime();

print HANDLE_LOG "\n\nPipe finished at $timestamp\n\n";
close HANDLE_LOG;


$log = "$filename1_short\_log_and_metrics.log";	

copy ($log, "all_assemblies/$log") or die "Copy failed: $!";  # copy log into folder with all assemblies

move($log, "$path_assembly/$log");	#move log into SPAdes-folder

print "\n";


exit;

