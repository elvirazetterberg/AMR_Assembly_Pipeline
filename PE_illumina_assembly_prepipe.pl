#!/usr/bin/perl
use File::Copy qw(copy);
use Cwd;

#   This script produces shell script that controls "PE_illumina_assembly_pipe.pl" which is an assembly pipe that uses SPAdes. Choose "yes" as last parameter to automatically start shell script.

#   1:  Start by putting all the fastq files in a folder. List their filenames with the terminal command "ls *.gz > input.txt" or "ls *.fastq > input.txt"
#   2:  Run script with: "perl PE_illumina_assembly_prepipe.pl input.txt trim/notrim 40 (wanted input coverage X) 1641481 (genome size Mbp) pilon/nopilon ariba/noariba kraken/nokraken	yes/no (run shell script)"
#   3:  If "no" was chosen as last parameter, the file "pipe_runs.sh" that is created by this script can be manually modified if different coverages or other settings needs to be changed for only a 
#       subset of the assemblies. If assembly work takes several hours it can be a good idea to check that the "pipe_runs.sh" has been created correctly before starting.
#   4:  Start shell script with "bash pipe_runs.sh"  


# All the parameters have to be given every time and in the same order as above, otherwise the script will return error message.
# add "perl " before script name if the script is not located in PATH.

# "trim" uses Trimmomatic. Pilon corrects assemblies. Ariba searches the reads for AMR-genes and mutations using CARD and ResFinder. Kraken classifies reads   
# according to taxon which can be used to quality check the reads, for species-identification and for metagenomics. If using ARIBA, also use Kraken to make sure AMR-genes are from correct genome.

# If no assembly is needed and for instance AMR-genes and/or species determination is the only goal, set the wanted coverage to "0". No SPAdes assembly will then be performed.

# Latest modification to script: 180302
# //Joakim Skarin



$inputfastqs = $ARGV[0];
$trim = $ARGV[1];
$wanted_coverage = $ARGV[2];
$genome_size = $ARGV[3];
$run_pilon = $ARGV[4];
$run_ariba = $ARGV[5];
$run_kraken = $ARGV[6];
$run_shell_script = $ARGV[7];


if ($trim ne "trim" && $trim ne "notrim") { die("\nFATAL ERROR: Trim ('trim/notrim') must be given\n");}

if ($run_pilon ne "pilon" && $run_pilon ne "nopilon") { die("\nFATAL ERROR: Pilon ('pilon/nopilon') must be given\n");}

#if ($prepare_consed ne "consed" && $prepare_consed ne "noconsed") { die("\nFATAL ERROR: consed ('consed/noconsed') must be given\n");}

if ($run_ariba ne "ariba" && $run_ariba ne "noariba") { die("\nFATAL ERROR: ARIBA ('ariba/noariba') must be given\n");}

if ($run_kraken ne "kraken" && $run_kraken ne "nokraken") { die("\nFATAL ERROR: Kraken ('kraken/nokraken') must be given\n");}

if ($run_shell_script ne "yes" && $run_shell_script ne "no") { print("\nWill assume \$run_shell_script = no\n"); $run_shell_script = 0;}

$fastq_counter = 0;

open(SHELL_SCRIPT,">"."pipe_runs.sh") || die("\nFATAL ERROR: could not open file:  pipe_runs.sh \n");
open(FASTQ_FILE,$inputfastqs) || die("FATAL ERROR: could not open file: $inputfastqs");

while (<FASTQ_FILE>) { 
	$fastq_counter += 1;

	if ($fastq_counter eq 1) {
		$fastq_1 = $_;
		chomp ($fastq_1);
} 
	if ($fastq_counter eq 2) {
		$fastq_2 = $_;
		chomp ($fastq_2);		
		print SHELL_SCRIPT "PE_illumina_assembly_pipe.pl $fastq_1 $fastq_2 $trim $wanted_coverage $genome_size $run_pilon $run_ariba $run_kraken\n"; 		
		$fastq_counter = 0;
  	}

}

close(SHELL_SCRIPT);
close(FASTQ_FILE);


system ("chmod +x pipe_runs.sh");


if ($run_shell_script eq "yes") {

	system ("bash pipe_runs.sh");
}


exit;

