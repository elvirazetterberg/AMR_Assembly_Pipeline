from curses import longname
import sys
import os
from datetime import datetime
import gzip
import re
# import pandas as pd

# Start by parsing the following command through the terminal, choosing only one option in each case:
# 'python assembly_pipeline_v4.py infile1/folder(???) infile2/none(???) here/there regular/parallel trim/notrim kraken/nokraken ariba/noariba wanted_coverage genome_size pilon/nopilon threads RAM'

# test run:
# python assembly_pipeline_v4.py SRR18825428_1.fastq.gz SRR18825428_2.fastq.gz here regular trim kraken ariba 40 124000000 pilon 40 0

'''OPTIONS'''
# infile1 / folder?? HOW DOES GENOME COVERAGE WORK THEN?
# infile2 / None ??
# - here/there: Where should all outputs be saved? If 'here' a new directory is created in 
# the current directory. If 'there' a path will be asked for.
# - regular/parallel: regular means running only one strain, parallel means running multiple strains
# - trim/notrim: trim means we run fastp, notrim means that we don't
# - kraken/nokraken: choose whether kraken should be run or not
# - ariba/noariba: choose whether to align AMR-genes with ariba
# - wanted_coverage: what coverage is requested? If 0, no assembly is performed
# - genome_size: what is the genome size of the input raw reads?
# - pilon/nopilon: choose whether to run pilon or not. Does not run if spades does not run (0 wanted coverage)
# - threads: maximum threads available
# - RAM: how much RAM that is available


def directory(date, time, there = False):
    
    ''' Function to create directory where all outputs from the pipeline are placed. 
    Date and time specific'''

    # Change the format of time from eg. 09:13:07.186006 to 09h13m07s
    stringtime = time[:2]+'h'+time[3:5]+'m'+time[6:8]+'s'

    # Choose path of directory
    if there:
        print('You requested to save all output files in another directory.')
        path = input('New path: ')
    else:
        path = os.getcwd()

    # Rename directory with date and time
    namedir = 'assembly_' + date + '_' + stringtime

    finalpath = os.path.join(path, namedir)

    os.mkdir(finalpath)
    
    return finalpath

def currenttime():
    '''Function that returns a string with the current time.'''
    time = str(datetime.time(datetime.now()))
    return time

def shortname(filename):
    '''Function that take a filename and returns a shorter version 
    including only the first continuous word-number sequence.'''
    short = re.search('[a-zA-Z1-9]+', filename).group()
    return short

def fastp_func(infile1, infile2, common_name):
    '''Function that takes two raw reads fastq files, one forward (1, right) and one reverse(2, left)
    and returns two trimmed fastq files as well as quality control documentation.'''

    loglines = f'Fastp started with {infile1} and {infile2}\n'

    outfile1 = f'out_fastp_{common_name}_1.fq.gz'
    outfile2 = f'out_fastp_{common_name}_2.fq.gz'

    fastpinput = 'fastp -i ' + infile1 + ' -I ' + infile2 + ' -o ' + outfile1 + ' -O ' + outfile2

    terminaltext = ' | tee fastp_func.txt'

    os.system(fastpinput+terminaltext)

    loglines += 'Fastp complete. Four output files returned:\n'
    loglines += f'{outfile1} \n{outfile2} \nfastp.html \nfastp.json \n\n'
    return outfile1, outfile2, loglines, 'fastp_func.txt'

def kraken_func(infile1, infile2, threads, common_name, path_kraken):
    ''' Function that runs Kraken on two raw reads fastq files, one forward (1, right) and one reverse(2, left), 
    in order to assign taxonomic labels to the sequences'''
                
    loglines = f'Kraken started with {infile1} and {infile2} as input with {threads} threads available \n' 

    kraken_output = f'out_kraken_{common_name}.out'
    kraken_report = f'report_kraken_{common_name}.report'

    krakeninput = f'kraken2 --db {path_kraken} --threads {threads} --output {kraken_output} --report {kraken_report} --paired {infile1} {infile2}'
    
    os.system(krakeninput)

    loglines += f'Kraken run finished. Two output files returned:\n'
    loglines += f'{kraken_output} \n{kraken_report}'
    return kraken_output, kraken_report, loglines

def reads_for_coverage(fastq_file, wanted_coverage, genome_size):
    '''Function that checks whether the requested coverage can be reached with the input
    files, returning the maximum coverage if this is not the case.'''
    loglines = (f'Running: reads_for_coverage')
    loglines += f'Checking if coverage can be achieved \n\n'

    bases_needed = int(wanted_coverage*genome_size/2)
    
    loglines = f'To achieve {wanted_coverage} X, {bases_needed} bases are needed from each fastq-file\n'

    total_bases = 0
    read_counter = 0
    row_counter = 1 # goes between 1 and 4
    
    with gzip.open(fastq_file, 'rt') as file:
        for line in file:
            if '@' in line:
                lenlist = re.findall('(length=[1-9]+)', line)
                if len(lenlist) > 0:
                    lenlist2 = lenlist[0].split('=')
                    readlength = int(lenlist2[1])
                    total_bases += readlength

            elif readlength == 0 and row_counter == 2:
                readlength = len(line)
                total_bases += readlength

            elif row_counter == 4:
                readlength = 0
                read_counter += 1

        # Give log output if the coverage can be achieved
            if total_bases >= bases_needed:
                loglines += f'Needed bases: {bases_needed} (per direction) which amounts to {read_counter} reads from fastq_1 which is {total_bases} bases\n\n'
                coverage = wanted_coverage
                break

            row_counter = row_counter%4 + 1 # makes the counter loop between 1 and 4


    # Give log output if the coverage CANNOT be achieved, and estimate new coverage
    if total_bases < bases_needed:
        loglines += f'There are not enough bases to achieve {wanted_coverage} X coverage.\n\n"'
        available_coverage = int(2*total_bases/genome_size)
        loglines += f'Using an estimated coverage of {available_coverage} X instead which amounts to {read_counter} reads from fastq_1 which is {total_bases} bases\n\n'
        coverage = available_coverage

    reads_needed = read_counter
    
    loglines += f'Function finished.\nOutputs: coverage {coverage}, reads needed {reads_needed}\n\n'

    return coverage, reads_needed, loglines

def shorten_fastq(fastq1_file, fastq2_file, reads_needed, common_name):
    '''Function that shortens the fastq files to only be long enough to reach 
    the requested coverage.'''

    loglines = f'shorten_fastq started to shorten {fastq1_file} and {fastq2_file} to only match wanted coverage.\n\n'

    lines_needed = reads_needed*4
    newname1 = f'X_{common_name}_1.fq.gz'
    newname2 = f'X_{common_name}_2.fq.gz'
    
    with gzip.open(fastq1_file, 'rt') as trim_me: # maybe change to 'rb'
        newfile = ''
        for i, line in enumerate(trim_me):
            newfile += line
            if i == lines_needed:
                break
    
    with gzip.open(newname1, 'wt') as one:
        one.write(newfile)

    with gzip.open(fastq2_file, 'rt') as trim_me:
        newfile = ''
        for i, line in enumerate(trim_me):
            newfile += line
            if i == lines_needed:
                break

    with gzip.open(newname2, 'wt') as one:
        one.write(newfile)

    loglines += f'Shortening complete.\nOutputs: {newname1}, {newname2}.\n\n'

    return newname1, newname2, loglines

def spades_func(file1, file2, path_spades, common_name, finalpath, threads): # threads, RAM
    '''Function that runs SPAdes to assemble contigs from short reads.'''

    loglines = 'SPAdes started\n'

    # To make sure X_spades output is in the correct output directory. 
    # Pilon output will also be added here
    assembly_path = f'{finalpath}/{common_name}_assembly'

    # commandline = '#SBATCH -p node -n 1 \n'
    commandline = f'python {path_spades}/spades.py --careful -o {assembly_path} --pe1-1 {file1} --pe1-2 {file2} -t {threads}'
    os.system(commandline)
    #"spades.py --careful -o $filename1_short\_$wanted_coverage\X_spades --pe1-1 $read1_output --pe1-2 $read2_output -t $threads_available -m $RAM_available"

    # rename from contigs.fasta to fasta
    os.system(f'cp {assembly_path}/contigs.fasta {assembly_path}/{common_name}.fasta')
    loglines += f'"contigs.fasta"-file copied and renamed to be called "{common_name}.fasta"'

    loglines += 'SPAdes finished.\n'
    loglines += f'All output files can be found here: {assembly_path}\n\n'

    return assembly_path, loglines

def pilon_func(fastafile, fasta1, fasta2, common_name, threads, assembly_path):
    '''Function that runs Pilon on contigs-file from SPAdes to 
    polish and assemble further.'''
    
    current = os.getcwd()
    
    os.chdir(assembly_path)
    
    loglines = 'Pilon started\n'
    loglines += f'Input files: {fastafile}, {fasta1}, {fasta2}\n'

    bowtie_build = f'bowtie2-build -f --threads {threads} --quiet {fastafile} {common_name}'
    os.system(bowtie_build)

    # inputs the two shortened fasta-files, if available
    bowtie = f'bowtie2 -x {common_name} -1 {fasta1} -2 {fasta2} -S {common_name}.sam --phred33 --very-sensitive-local --no-unal -p {threads}'
    os.system(bowtie)

    os.system(f'samtools view -bh {common_name}.sam > {common_name}.bam')
    os.system(f'samtools sort {common_name}.bam -o {common_name}.sorted.bam')
    os.system(f'samtools index {common_name}.sorted.bam')

    time = currenttime()+'\n'
    loglines += f'Pilon 1.24 started at {time}'
    
    os.system(f'pilon --genome {common_name}.fasta --frags {common_name}.sorted.bam --output {common_name}.pilon --changes --threads {threads}')
    
    time = currenttime()+'\n'
    loglines += f'Pilon finished at {time}\n'
    
    loglines += f'Corrected fasta file created: {common_name}.pilon.fasta'

    os.chdir(current)

    return loglines

def info(spades_assembly):
    '''Function that uses an assembly-file from SPAdes of Pilon 
    and returns the metrics of that assembly.'''

    # Output som pandas table for att satta ihop alla strains till en sammanstalld csv-fil, 
    # och varje strain till var sin csv-fil

    loglines = 'Looking at the metrics of assembly_fasta\n\n'

    number_of_contigs, bases_in_contig, total_bases = 0, 0, 0
    contig_lengths = []
    non_base, number_AT, number_GC = 0, 0, 0
    contigs_over_1000 = 0

    # Loop through and get metrics
    with open(spades_assembly, 'r') as s:
        for line in s:
            if '>' in line:
                number_of_contigs += 1
                total_bases += bases_in_contig

                if bases_in_contig != 0:
                    contig_lengths.append(bases_in_contig)
                    if bases_in_contig >= 1000:
                        contigs_over_1000 += 1

                bases_in_contig = 0 # resets

            else:
                bases_in_contig += len(line)
                for base in line:
                    if base == 'A' or base == 'T':
                        number_AT += 1
                    elif base == 'G' or base == 'C':
                        number_GC += 1
                    else:
                        non_base += 1
    
    loglines += f'The number of contigs: {number_of_contigs}, the total number of bases: {total_bases}\n\n'

    sorted_contig_lengths = contig_lengths.sort()
    longest = sorted_contig_lengths[0]
    loglines += f'Longest contig: {longest}\n'
    loglines += f'Contigs longer than 1 kb: {contigs_over_1000}'

    # N50
    temp = 0
    while temp <= total_bases/2:
        for length in sorted_contig_lengths:
            temp += length
            N_50 = length
    
    loglines += f'N50: {N_50}\n'

    # GC-content
    GC = number_GC/(number_GC + number_AT)*100
    loglines += f'The GC-content of the sequence is {GC}. {non_base} non-base characters were excluded from GC-calculation\n'

    loglines += f'-----------------------Metrics finished-----------------------'

    # PLACE ALL INFO IN PANDAS TABLE
    data = {'Total nr bases': total_bases, 'Nr contigs': number_of_contigs, 'Longest contig': longest, 
    'Nr contigs > 1kb': contigs_over_1000, 'N50':N_50, 'GC-content': GC}

    info_df = pd.DataFrame(data)

    return info_df, loglines

def log_parse(string):
    time = currenttime()
    log_name="unique_name"
    os.system(f"echo {time}: '{string}' >> {log_name}.txt")
    return


def ariba_fun(infile1,infile2,db_ariba):
    for db_name in db_ariba[1:-1].split(','): #klumpigt? as sysargv makes input a string, it is separated into a list here. Also parallell should do all db at same time?
        
        # OBS when making parallell the naming of files must take this into account. Right now Im deleting the privious runs
    
        # if there's allready an existing db, lite fult gjort?
        # In the pearlpipeline it says "$test_CARD = "CARD_reference_dataset_downloaded"; What does it mean? 

        log_parse(f' Starting ariba with {db_name}')
        if os.path.exists(f'out.{db_name}.fa'): 
            log_parse(f'Database {db_name} allready downloaded\n')
            os.system(f"rm -rf out.run.{db_name}") # är detta smart sätt att göra det? Ska det läggas till i log?

        else: # if database not downloaded. troligen onädigt i framtiden
            rm_db=input(f'Please note that running this will remove all existing files starting with "out.{db_name}". [y]/n?') 
            if rm_db.lower().startswith("y")==False:
                log_parse(f'Answered no: NOT Removing all existing files starting with "out.{db_name}". Exiting ariba.')
                exit() 
            else:
                log_parse(f'Answered yes: Removing all existing files starting with "out.{db_name}". Proceeding with ariba.')
                os.system(f"rm -rf out.{db_name}*")

            log_parse(f'Downloading database {db_name}')
            os.system(f"ariba getref {db_name} out.{db_name} >> {logname}.txt")

            log_parse(f'Preparing references with prefix out.{db_name}')
            os.system(f"ariba prepareref -f out.{db_name}.fa -m out.{db_name}.tsv out.{db_name}.prepareref >> {logname}.txt")

        log_parse(f'Running ariba on {db_name}')
        os.system(f"ariba run out.{db_name}.prepareref {infile1} {infile2} out.run.{db_name} >> {logname}.txt")

    log_parse(f'Ariba done.')
    return

def create_log(finalpath, time, date):

    lines = 15*'-' + 'LOGFILE' + 15*'-' + '\n\n'
    lines += f'New directory created with the adress {finalpath}\n'
    lines += f'Directory created at {time} on {date}\n'
    lines += 'All outputs will be saved in the new directory.\n\n'
    os.system(f"echo '{lines}' > {logname}.txt")
    
    return

# function that runs multiple strains in parallel. Inputs are all sys.argv[]
# Return lines for logfile?
def parallelize():
    pass

# function that runs everything for only one strain. Inputs are all sys.argv[]
# Return lines for logfile?
def regular():
    pass


def main():
    # python3 logwrite_fun.py /home/alma/Documents/kandidat/genomes/SRR18825428_1.fastq /home/alma/Documents/kandidat/genomes/SRR18825428_2.fastq noariba [vfdb_core]

    # os.system('SBATCH -p node -n 1')

    infile1 = sys.argv[1]
    infile2 = sys.argv[2]
    new_location = False # will ask for directory location if True
    ariba = sys.argv[3] == 'ariba'
    db_ariba=sys.argv[4]


# Let's start this pipeline!
    time = currenttime()
    date = str(datetime.date(datetime.now()))
    
# make directory for output
    finalpath = directory(date, time, new_location)

# Create log file
    global logname
    logname= 'unique_name' #need to make unique during parallell. Kanske likt directory? stringtime = time[:2]+'h'+time[3:5]+'m'+time[6:8]+'s'
    create_log(finalpath,time,date)

# Ariba 
    if ariba:
        header= '='*15 +'\n'+'ARIBA'+ '='*15 +'\n'
        log_parse(header)
        ariba_fun(infile1,infile2, db_ariba)
    
#testlog

# Run in parallel
    # if parallel:
    # parallelize()

    #os.system(f'mv {logname}.txt '+str(finalpath))
    

if __name__ == '__main__':
    main()  
