from datetime import datetime
import gzip
from typing import final
from numba import njit
import pandas as pd
import os
import re
import sys

from pytz import common_timezones


# Start by parsing the following command through the terminal, choosing only one option in each case:
# 'python assembly_pipeline_v6.py infile1/folder(???) infile2/None(???) here/there trim/notrim kraken/nokraken ariba/noariba wanted_coverage genome_size pilon/nopilon threads RAM'

# test run:
# python assembly_pipeline_v6.py SRR18825428_1.fastq.gz SRR18825428_2.fastq.gz here trim kraken noariba [vfdb_core] 40 1743985 nopilon 40 0
# 
# python assembly_pipeline_v6.py SRR18825428test_1.fastq.gz SRR18825428test_2.fastq.gz here trim kraken noariba [vfdb_core] 40 1743985 nopilon 40 0

# Lokal Alma:
# python Pipeline/assembly_pipeline_v6.py /home/alma/Documents/kandidat/genomes/SRR18825428_1.fastq /home/alma/Documents/kandidat/genomes/SRR18825428_2.fastq here ntrim nkraken ariba [vfdb_core] 40 1743985 npilon 40 0


'''OPTIONS'''
# - infile1 / directory: enter a directory to multiple files if parallelization is wanted. Will use same wanted coverage, however!
# - infile2 / None: enter None if parallelize.
# - here/there: Where should all outputs be saved? If 'here' a new directory is created in 
# the current directory. If 'there' a path will be asked for.
# - trim/notrim: trim means we run fastp, notrim means that we don't
# - kraken/nokraken: choose whether kraken should be run or not
# - ariba/noariba: choose whether to align AMR-genes with ariba
# - [vfdb_core]: list of AMR-databases for ariba, without spaces
# - wanted_coverage: what coverage is requested? If 0, no assembly is performed.
# - genome_size: what is the genome size of the organism?
# - pilon/nopilon: choose whether to run pilon or not. Does not run if spades does not run (0 wanted coverage)
# - threads: maximum threads available
# - RAM: how much RAM that is available


def script_log(finalpath, logname):
    os.system('script {logname}'.format(logname=logname)) # use logname???
    os.system('cd /proj/uppmax2022-2-14/private/campy_pipeline\n module load conda\nexport CONDA_ENVS_PATH=/proj/uppmax2022-2-14/private/campy_pipeline/environments\nconda activate pipeline_env\n')
    os.system('mv logfilescript.txt {finalpath}'.format(finalpath=finalpath))

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

def create_log(finalpath, time, date, logname):

    lines = 15*'-' + 'LOGFILE' + 15*'-' + '\n\n'
    lines += 'New directory created with the adress {finalpath}\n'.format(finalpath = finalpath)
    lines += 'Directory created at {time} on {date}\n'.format(time=time, date=date)
    lines += 'All outputs will be saved in the new directory.\n\n'
    os.system("echo '{lines}' > {logname}".format(lines=lines, logname=logname))
    
    return

def log_parse(string, logpath = ''):
    time = currenttime()
    os.system("echo {time}: '{string}\n' >> {logpath}{logname}".format(time=time, string=string, logpath=logpath, logname=logname))
    return

# @njit(parallel=True)
def ariba_fun(infile1,infile2,db_ariba):
    for db_name in db_ariba[1:-1].split(','): #klumpigt? as sysargv makes input a string, it is separated into a list here. Also parallell should do all db at same time?
        
        # OBS when making parallell the naming of files must take this into account. Right now Im deleting the privious runs
    
        # if there's allready an existing db, lite fult gjort?
        # In the pearlpipeline it says "$test_CARD = "CARD_reference_dataset_downloaded"; What does it mean? 

        log_parse(' Starting ariba with {db_name}'.format(db_name=db_name))
        if os.path.exists('out.{db_name}.fa'.format(db_name=db_name)): 
            log_parse('Database {db_name} allready downloaded\n'.format(db_name=db_name))
            os.system("rm -rf out.run.{db_name}".format(db_name=db_name)) # är detta smart sätt att göra det? Ska det läggas till i log?

        else: # if database not downloaded. troligen onädigt i framtiden
            rm_db=input('Please note that running this will remove all existing files starting with "out.{db_name}". [y]/n?'.format(db_name=db_name)) 
            if rm_db.lower().startswith("y")==False:
                log_parse('Answered no: NOT Removing all existing files starting with "out.{db_name}". Exiting ariba.'.format(db_name=db_name))
                exit() 
            else:
                log_parse('Answered yes: Removing all existing files starting with "out.{db_name}". Proceeding with ariba.'.format(db_name=db_name))
                os.system("rm -rf out.{db_name}*".format(db_name=db_name))

            log_parse('Downloading database {db_name}'.format(db_name=db_name))
            os.system("ariba getref {db_name} out.{db_name} >> {logname}".format(db_name=db_name, logname=logname))

            log_parse('Preparing references with prefix out.{db_name}'.format(db_name=db_name))
            os.system("ariba prepareref -f out.{db_name}.fa -m out.{db_name}.tsv out.{db_name}.prepareref >> {logname}".format(db_name=db_name, logname=logname))

        log_parse('Running ariba on {db_name}'.format(db_name=db_name))
        os.system("ariba run out.{db_name}.prepareref {infile1} {infile2} out.run.{db_name} >> {logname}".format(db_name=db_name, infile1=infile1, infile2=infile2, logname=logname))

    log_parse('Ariba done.\n')
    return

# @njit(parallel=True)
def fastp_func(infile1, infile2, common_name):
    '''Function that takes two raw reads fastq files, one forward (1, right) and one reverse(2, left)
    and returns two trimmed fastq files as well as quality control documentation.'''

    log_parse('Fastp started with {infile1} and {infile2}\n'.format(infile1=infile1, infile2=infile2))

    outfile1 = 'out_fastp_{common_name}_1.fq.gz'.format(common_name=common_name)
    outfile2 = 'out_fastp_{common_name}_2.fq.gz'.format(common_name=common_name)

    fastpinput = 'fastp -i {infile1} -I {infile2} -o {outfile1} -O {outfile2}'.format(infile1=infile1, infile2=infile2, outfile1=outfile1, outfile2=outfile2)

    # os.system(fastpinput) # I dont know if this generates outpu, but in that case I should be parsed into logfile like below
    # log_parse(fastpinput) 
    os.system('{fastpinput} >> {logname}'.format(fastpinput=fastpinput, logname=logname))
    
    log_parse('Fastp complete. Four output files returned:\n{outfile1} \n{outfile2} \nfastp.html \nfastp.json \n'.format(outfile1=outfile1, outfile2=outfile2))
    
    return outfile1, outfile2

# @njit(parallel=True)
def kraken_func(infile1, infile2, threads, common_name, path_kraken):
    ''' Function that runs Kraken on two raw reads fastq files, one forward (1, right) and one reverse(2, left), 
    in order to assign taxonomic labels to the sequences'''
                
    log_parse('Kraken started with {infile1} and {infile2} as input with {threads} threads available \n'.format(infile1=infile1, infile2=infile2, threads=threads))
    kraken_output = 'out_kraken_{common_name}.out'.format(common_name=common_name)
    kraken_report = 'report_kraken_{common_name}.report'.format(common_name=common_name)

    krakeninput = 'kraken2 --db {path_kraken} --threads {threads} --output {kraken_output} --report {kraken_report} --paired {infile1} {infile2}'.format(path_kraken=path_kraken, threads=threads,kraken_output=kraken_output, kraken_report=kraken_report, infile1=infile1, infile2=infile2)
    
    os.system(krakeninput)# I dont know if this generates outpu, but in that case I should be parsed into logfile like below
    #os.system(f'{krakeninput} >> {logname}') 

    log_parse('Kraken run finished. Two output files returned:\n')
    log_parse('{kraken_output} \n{kraken_report}'.format(kraken_output=kraken_output, kraken_report=kraken_report))
    return kraken_output, kraken_report

# @njit(parallel=True)
def reads_for_coverage(fastq_file, wanted_coverage, genome_size):
    '''Function that checks whether the requested coverage can be reached with the input
    files, returning the maximum coverage if this is not the case.'''
    log_parse('Running: reads_for_coverage')
    log_parse('Checking if coverage can be achieved \n\n')

    bases_needed = int(wanted_coverage*genome_size/2)
    
    log_parse('To achieve {wanted_coverage} X, {bases_needed} bases are needed from each fastq-file\n'.format(wanted_coverage=wanted_coverage, bases_needed=bases_needed))
    log_parse('Checking if wanted coverage can be achieved...')
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
                log_parse('Coverage can be reached! It amounts to {read_counter} reads from fastq_1 which is {total_bases} bases\n\n'.format(read_counter=read_counter, total_bases=total_bases))
                coverage = wanted_coverage
                break

            row_counter = row_counter%4 + 1 # makes the counter loop between 1 and 4


    # Give log output if the coverage CANNOT be achieved, and estimate new coverage
    if total_bases < bases_needed:
        log_parse('There are not enough bases to achieve {wanted_coverage} X coverage.\n"'.format(wanted_coverage=wanted_coverage))
        available_coverage = int((2*total_bases)/genome_size)
        log_parse('Using an estimated coverage of {available_coverage} X instead which amounts to {read_counter} reads and {total_bases} bases from fastq_1\n\n'.format(available_coverage=available_coverage, read_counter=read_counter, total_bases=total_bases))
        coverage = available_coverage

    reads_needed = read_counter
    
    log_parse('Function finished.\nOutputs: coverage {coverage}, reads needed {reads_needed}\n\n'.format(coverage=coverage, reads_needed=reads_needed))

    return coverage, reads_needed

# @njit(parallel=True)
def shorten_fastq(fastq1_file, fastq2_file, reads_needed, common_name):
    '''Function that shortens the fastq files to only be long enough to reach 
    the requested coverage.'''

    log_parse('shorten_fastq started to shorten {fastq1_file} and {fastq2_file} to only match wanted coverage.\n\n'.format(fastq1_file=fastq1_file, fastq2_file=fastq2_file))

    lines_needed = reads_needed*4
    newname1 = 'X_{common_name}_1.fastq.gz'.format(common_name=common_name)
    newname2 = 'X_{common_name}_2.fastq.gz'.format(common_name=common_name)
    
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

    log_parse('Shortening complete.\nOutputs: {newname1}, {newname2}.\n\n'.format(newname1=newname1, newname2=newname2))

    return newname1, newname2

# @njit(parallel=True)
def spades_func(file1, file2, path_spades, common_name, finalpath, threads): # threads, RAM
    '''Function that runs SPAdes to assemble contigs from short reads.'''

    log_parse('SPAdes started\n')

    # To make sure X_spades output is in the correct output directory. 
    # Pilon output will also be added here
    assembly_path = '{finalpath}/{common_name}_assembly'.format(finalpath=finalpath, common_name=common_name)

    # commandline = '#SBATCH -p node -n 1 \n'
    commandline = 'python {path_spades}/spades.py --careful -o {assembly_path} --pe1-1 {file1} --pe1-2 {file2} -t {threads}'.format(path_spades=path_spades, assembly_path=assembly_path, file1=file1, file2=file2, threads=threads)
    os.system(commandline)
    #"spades.py --careful -o $filename1_short\_$wanted_coverage\X_spades --pe1-1 $read1_output --pe1-2 $read2_output -t $threads_available -m $RAM_available"

    # rename from contigs.fasta to fasta to work with pilon
    os.system('cp {assembly_path}/contigs.fasta {assembly_path}/{common_name}.fasta').format(assembly_path=assembly_path, common_name=common_name)
    log_parse('"contigs.fasta"-file copied and renamed to be called "{common_name}.fasta"'.format(common_name=common_name))

    log_parse('SPAdes finished.\n')
    log_parse('All output files can be found here: {assembly_path}\n\n'.format(assembly_path=assembly_path))

    return assembly_path

# @njit(parallel=True)
def pilon_func(fastafile, fasta1, fasta2, common_name, threads, assembly_path):
    '''Function that runs Pilon on contigs-file from SPAdes to 
    polish and assemble further.'''
    
    current = os.getcwd()+'/'
    
    os.chdir(assembly_path)
    
    log_parse('Pilon started', current)
    log_parse('Input files: {fastafile}, {fasta1}, {fasta2}'.format(fastafile=fastafile, fasta1=fasta1, fasta2=fasta2), current)

    bowtie_build = 'bowtie2-build -f --threads {threads} --quiet {fastafile} {common_name}'.format(threads=threads, fastafile=fastafile, common_name=common_name)
    os.system(bowtie_build)

    # inputs the two shortened fasta-files, if available
    bowtie = 'bowtie2 -x {common_name} -1 {fasta1} -2 {fasta2} -S {common_name}.sam --phred33 --very-sensitive-local --no-unal -p {threads}'.format(common_name=common_name, fasta1=fasta1, fasta2=fasta2, threads=threads)
    os.system(bowtie)

    os.system('samtools view -bh {common_name}.sam > {common_name}.bam'.format(common_name=common_name))
    os.system('samtools sort {common_name}.bam -o {common_name}.sorted.bam'.format(common_name=common_name))
    os.system('samtools index {common_name}.sorted.bam'.format(common_name=common_name))

    time = currenttime()+'\n'
    log_parse('Pilon 1.24 started at {time}'.format(time=time), current)
    
    os.system('pilon --genome {common_name}.fasta --frags {common_name}.sorted.bam --output {common_name}.pilon --changes --threads {threads}'.format(common_name=common_name, threads=threads))
    
    log_parse('Pilon finished\n', current) #removed at time. keep? 
    
    log_parse('Corrected fasta file created: {common_name}.pilon.fasta'.format(common_name=common_name), current)

    os.chdir(current)

    return 

# @njit(parallel=True)
def info(spades_assembly):
    '''Function that uses an assembly-file from SPAdes of Pilon 
    and returns the metrics of that assembly.'''

    # Output som pandas table for att satta ihop alla strains till en sammanstalld csv-fil, 
    # och varje strain till var sin csv-fil

    log_parse( 'Looking at the metrics of assembly_fasta\n')

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
    
    log_parse('The number of contigs: {number_of_contigs}, the total number of bases: {total_bases}\n'.format(number_of_contigs=number_of_contigs, total_bases=total_bases))

    contig_lengths.sort()
    longest = contig_lengths[-1]
    log_parse('Longest contig: {longest}\n'.format(longest=longest))
    log_parse('Contigs longer than 1 kb: {contigs_over_1000}'.format(contigs_over_1000=contigs_over_1000))

    # N50
    temp = 0
    for length in contig_lengths:
        temp += length
        N_50 = length
        if temp >= (total_bases/2):
            break
    
    log_parse('N50: {N_50}\n'.format(N_50=N_50))

    # GC-content
    GC = round(number_GC*100/(number_GC + number_AT),2)
    log_parse('The GC-content of the sequence is {GC}%. {non_base} non-base characters were excluded from GC-calculation\n'.format(GC=GC, non_base=non_base))

    log_parse('-----------------------Metrics finished-----------------------')

    # PLACE ALL INFO IN PANDAS TABLE
    data = {'Total nr bases': total_bases, 'Nr contigs': number_of_contigs, 'Longest contig': longest, 
    'Nr contigs > 1kb': contigs_over_1000, 'N50':N_50, 'GC-content': GC}

    info_df = pd.DataFrame(data)

    return info_df

# function that runs everything for only one strain. Inputs are all sys.argv[]
def regular(path, infile1, infile2, run_fastp, kraken, ariba, db_ariba, run_spades, wanted_coverage, genome_size, pilon, threads, common_name):
    
    # script_log(path)

    time = currenttime()
    date = str(datetime.date(datetime.now()))

    path_tools = '/proj/uppmax2022-2-14/private/campy_pipeline/assembly/verktyg'
    path_spades = path_tools + '/SPAdes-3.15.4-Linux/bin'
    path_kraken = path_tools + '/minikraken2_v1_8GB'

    os.system('cp {infile1} {infile2} {path}'.format(infile1=infile1, infile2=infile2, path=path))
    os.chdir(path)


    # Create log file
    global logname
    logname = 'logfile.txt' # maybe not the smartest when parallell? but otherwise need to feed it into every function
    # Rename log file with date and time --- easier to just use logfile.txt globally?
    # stringtime = time[:2]+'h'+time[3:5]+'m'+time[6:8]+'s'
    # logname = 'LOGFILE' + date + '_' + stringtime
    create_log(path, time, date, logname)
    print('Pipeline started, please refer to logfile "{logname}" for updates.'.format(logname)) # add path to logfile later 

# Ariba 
    if ariba:
        header= '\n'+'='*15 +'ARIBA'+ '='*15 +'\n'
        log_parse(header)
        ariba_fun(infile1,infile2, db_ariba)
        os.system("ariba summary out_sum out.run.*/report.tsv")

# Fastp
    if run_fastp:
        header= '\n'+'='*15 +'FASTP'+ '='*15 +'\n'
        log_parse(header)
        outfile1_trim, outfile2_trim = fastp_func(infile1, infile2, common_name)

        infile1 = outfile1_trim
        infile2 = outfile2_trim

# Kraken
    if kraken:
        header= '\n'+'='*15 +'KRAKEN'+ '='*15 +'\n'
        log_parse(header)
        kraken_output, kraken_report = kraken_func(infile1, infile2, threads, common_name, path_kraken)

# Number of reads to match the wanted coverage
    if run_spades:
        header= '\n'+'='*15 +'READS FOR COVERAGE, SPADES'+ '='*15 +'\n'
        log_parse(header)
        coverage, reads_needed = reads_for_coverage(infile1, wanted_coverage, genome_size)
    else:
        coverage = 0

# Shortening fastq-files if the coverage can be reached with less
    if coverage > wanted_coverage:
        outfile1_shorten, outfile2_shorten = shorten_fastq(infile1, infile2, reads_needed, common_name)        
        infile1 = outfile1_shorten
        infile2 = outfile2_shorten
        os.system('mv {outfile1_shorten} {outfile2_shorten} {path}'.format(outfile1_shorten=outfile1_shorten, outfile2_shorten=outfile2_shorten, path=path))
        log_parse('Shortened fastq files from shorten_fastq function moved to directory\n\n')

# Spades
    if run_spades:
        header= '\n'+'='*15 +'SPADES'+ '='*15 +'\n'
        log_parse(header)
        assembly_path = spades_func(infile1, infile2, path_spades, common_name, path, threads)

# Pilon
    if pilon:
        header= '\n'+'='*15 +'PILON'+ '='*15 +'\n'
        log_parse(header)
        fastafile = '{assembly_path}/{common_name}.fasta'.format(assembly_path=assembly_path, common_name=common_name)
        pilon_func(fastafile, infile1, infile2, common_name, threads, assembly_path)
        
        # input file found here: assembly_path/SRR18825428.fasta

# Info/metrics
    if run_spades:
        header= '\n'+'='*15 +'INFO/METRICS, SPADES'+ '='*15 +'\n'
        log_parse(header)

        from_spades = '{assembly_path}/{common_name}.fasta'.format(assembly_path=assembly_path, common_name=common_name)
        
        info_df = info(from_spades)
        info_df.to_csv(header = True, path_or_buf = '{path}/{common_name}_info.csv'.format(path=path, common_name=common_name))
        
        # If we have multiple info_df then use pd.concat([info_df1, info_df2], axis=0) to stack the 2nd below the 1st.
        # This is useful when running in parallel.
        
        # Save info_df INSTEAD KEEP AS DF AND CONCAT WITH KRAKEN AND ALIGNMENT
        # info_df.to_csv(os.PathLike(f'{path}/{common_name}_metrics'))
    
    
    # return info_df to be able to concatenate if running in parallel.
    # 1. change ariba out.db.tsv and kraken report to csv
    # 2. change csv to df
    # all_df = pd.concat([info_df, kraken_df, alignment_df], axis=1) # stack in one row
    # all_df.to_csv(os.PathLike(f'{path}/{common_name}_metrics'))


# function that runs multiple strains in parallel. Inputs are all sys.argv[]
@njit(parallel=True)
def parallelize(finalpath, file_directory, run_fastp, kraken, ariba, db_ariba, run_spades, wanted_coverage, genome_size, pilon, threads, common_name):
    '''Function that takes a directory pf forward and reverse files to run the pipeline with in parallel.'''
    
    current = os.getcwd()
    os.chdir(file_directory)
    os.system("ls *.gz > input.txt")
    # go back
    os.chdir(current)

    dirlist = []
    with open('{file_directory}/input.txt'.format(file_directory=file_directory), 'r') as inp:
        linelist = inp.readlines()
        for i in range(0, len(linelist), 2):
            common_name = shortname(i)
            path = '{finalpath}/{common_name}'.format(finalpath=finalpath, common_name=common_name)
            os.mkdir(path)
            dirlist.append(path)
            regular(path, '{file_directory}/{linelist[i]}'.format(file_directory=file_directory, linelist=linelist), '{file_directory}/{linelist[i+1]}'.format(file_directory=file_directory, linelist=linelist), run_fastp, kraken, ariba, db_ariba, run_spades, wanted_coverage, genome_size, pilon, threads, common_name)

    # os.system(f'cd {finalpath}') # change back to finalpath to place all info in

    pass

def main():
    # os.system('SBATCH -p node -n 1')
    """
    path/to/file1 path/to/file2 here nopar notrim nokraken ariba [db1, db2] 0 size nopilon thr ram
    """
    infile1 = sys.argv[1] # 
    infile2 = sys.argv[2]
    new_location = sys.argv[3] == 'there' # will ask for directory location if True
    run_fastp = sys.argv[4] == 'trim' # will run fastp if True
    kraken = sys.argv[5] == 'kraken'
    ariba = sys.argv[6] == 'ariba'
    db_ariba = sys.argv[7] 
    wanted_coverage = int(sys.argv[8]) # if wanted coverage == 0, then don't run spades
    genome_size = int(sys.argv[9])
    pilon = sys.argv[10] == 'pilon'
    threads = sys.argv[11]
    RAM = sys.argv[12] # this has not been implemented

    run_spades = wanted_coverage != 0
    common_name = shortname(infile1) # until here only work if not parallel

    if pilon and run_spades == False: # Since pilon requires spades output, this 
        pilon = False
        pilon_lines = 'Pilon not run since SPAdes was not run (!)\n\n'

    ''' -------------------CHANGE !?!?!?!----------------------- '''
# Hardcoded, location of non-conda tools


# Let's start this pipeline!
    time = currenttime()
    date = str(datetime.date(datetime.now()))
    
# make directory for output
    finalpath = directory(date, time, new_location)

    '''CONTINUE FROM HERE'''
    if os.path.isdir(infile1):
        parallelize(finalpath, infile1, run_fastp, kraken, ariba, db_ariba, run_spades, wanted_coverage, genome_size, pilon, threads, common_name)
    else:
        regular(finalpath, infile1, infile2, run_fastp, kraken, ariba, db_ariba, run_spades, wanted_coverage, genome_size, pilon, threads, common_name)
    # if infile1 == directory (hur)?
        # parallelize()
    # else:
    #     regular()


    # OBS OBS convert to csv later instead
    # if run_spades:
    # os.system(f'mv {info_df} {finalpath}')
    # log.writelines(f'Dataframe moved')

# Close and move the logfile to the correct directory

#    os.system('mv {logname} '+str(finalpath))
    

if __name__ == '__main__':
    main()  
