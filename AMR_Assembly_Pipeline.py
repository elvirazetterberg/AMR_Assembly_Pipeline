from datetime import datetime
import gzip
import pandas as pd
import os
import re
import sys
import concurrent.futures as future
import glob


"""

                    HI  AND  WELCOME  TO  THE _____
             ____  _ ____  _____  _     _ __   _  _____  
            |  _ \| |  _ \|  ___|| |   | |  \ | ||  ___| 
            | |_) | | |_) | |__  | |   | |   \| || |___  
            |  __/| |  __/| |___ | |___| | |\ \ || |___   
            |_|   |_|_|   |_____||_____|_| | \__||_____| 

                            UU x SVA


Before running:
    * Make sure all tools mentioned in the sections 'The pipeline tools and where to find them' are installed properly.
    * Hard code the paths to the SPAdes-tool and Minikraken2 database in the function 'regular'

*************Now you should be ready to test the pipeline!*************

Start by parsing the following command through the terminal, choosing only one option in each case:
'python AMR_Assembly_Pipeline.py infile1/directory infile2/None here/path trim/notrim kraken/nokraken ariba/noariba wanted_coverage genome_size pilon/nopilon threads RAM'



REFERENCE:
The work of third year students at the Master's Programme in Molecular Biotechnology Engineering (X) at Uppsala University 2022. 
GitHub: https://github.com/elvirazetterberg/Pipeline.git  #change name 
Contributors: Alma Nilsson, Corinne Olivero, Elvira Zetterberg, Evelina Andersson, Julia Sulyaeva, Moa Qvarnlof

"""




'''OPTIONS'''
# - infile1 / directory: enter a directory to multiple short read files to run the pipeline in 
# parallel. Will use same wanted coverage, however!
# - infile2 / None: enter None if a directory was entered as infile1.
# - here/path: Where should all outputs be saved? If 'here' a new directory is created in
# the current directory. If a path is given the output will be saved there.
# - trim/notrim: trim means we run fastp, notrim means that we don't.
# - kraken/nokraken: choose whether kraken should be run or not.
# - ariba/noariba: choose whether to align AMR-genes with ariba.
# - [vfdb_core]: list of AMR-databases for ariba, without spaces.
# - wanted_coverage: what coverage is requested? If 0, no assembly is performed.
# - genome_size: what is the genome size of the organism?
# - pilon/nopilon: choose whether to run pilon or not. Does not run if no spades (0 wanted coverage).
''' THIS VERSION CANNOT RUN PILON.'''
# - threads: maximum threads available.
# - RAM: maximum RAM available.

def directory(date, time, here):
    
    ''' Function to create directory where all outputs from the pipeline are placed. 
    Date and time specific'''

    # Change the format of time from eg. 09:13:07.186006 to 09h13m07s
    stringtime = time[:2]+'h'+time[3:5]+'m'+time[6:8]+'s'

    # Choose path of directory
    if not here:
        log_parse('You requested to save all output files in another directory.')
        path = new_location
    else:
        path = os.getcwd()

    # Rename directory with date and time
    namedir = 'pipeline_output_' + date + '_' + stringtime

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
    splitit = filename.split('/')
    name = splitit[-1]
    short = re.search('[a-zA-Z0-9]+', name).group()
    return short

def create_log(path, time, date, infile1, infile2):

    lines = f'                                              \n\
                    HI  AND  WELCOME  TO  THE _____         \n\
             ____  _ ____  _____  _     _ __   _  _____     \n\
            |  _ \| |  _ \|  ___|| |   | |  \ | ||  ___|    \n\
            | |_) | | |_) | |__  | |   | |   \| || |___     \n\
            |  __/| |  __/| |___ | |___| | |\ \ || |___     \n\
            |_|   |_|_|   |_____||_____|_| | \__||_____|    \n\
                                                            \n\
                            UU x SVA                        \n\n\n'

    lines += 15*'-' + 'LOGFILE' + 15*'-' + '\n\n'
    lines += f'Pipeline called with following arguments:\n'
    lines += f'\t infile1: {infile1}, infile2: {infile2}, new_location: {new_location},\n'
    lines += f'\t run_fastp: {run_fastp}, kraken: {kraken}, ariba: {ariba}, db_ariba: {db_ariba}, wanted_coverage: {wanted_coverage},\n'
    lines += f'\t genome_size: {genome_size}, pilon: {pilon}, threads: {threads}\n\n'
    lines += f'The following packages and versions used:\n'
    lines += os.popen("conda list | awk '/^python /{print $1\"\t\"$2}'").read().strip() + '\n'
    lines += os.popen("conda list | awk '/^fastp /{print $1\"\t\"$2}'").read().strip() + '\n'
    lines += os.popen("conda list | awk '/^kraken2 /{print $1\"\t\"$2}'").read().strip() + '\n'
    lines += os.popen("conda list | awk '/^spades /{print $1\"\t\"$2}'").read().strip() + '\n'
    lines += os.popen("conda list | awk '/^bowtie2 /{print $1\"\t\"$2}'").read().strip() + '\n'
    lines += os.popen("conda list | awk '/^ariba /{print $1\"\t\"$2}'").read().strip() + '\n\n'
    lines += f'New directory created with the adress {path}\n'
    lines += f'Directory created at {time} on {date}\n'
    lines += 'All outputs will be saved in the new directory.\n\n'
    log_parse(lines, path)
    
    return

def log_parse(string, logpath):
    time = currenttime()
    os.system(f"echo {time}: '{string}\n' >> {logpath}/{logname}")
    
    return

def ariba_fun(path, infile1, infile2, db_ariba):
    print(base_dir)
    
    for db_name in db_ariba: 
        log_parse(f' Starting ariba with {db_name}', path)
        if os.path.exists(f'{base_dir}/out.{db_name}.fa'): #if databases already downloaded
            log_parse(f'Database {db_name} already downloaded', path)
            os.system(f"rm -rf {base_dir}/out.run.{db_name}") # OBS need warning?

        else: # if database not downloaded.
            os.system(f"rm -rf {base_dir}/out.{db_name}*")
            log_parse(f'Downloading database {db_name}', path)
            os.system(f"ariba getref {db_name} {base_dir}/out.{db_name} >> {path}/{logname}") # >> kommer inte funka, finns inte logname i finalpath som vi står i
            # Provade >> {path}/{logname} istället för >> {logname}. Kanske funkar

            log_parse(f'Preparing references with prefix out.{db_name}', path)
            os.system(f"ariba prepareref -f {base_dir}/out.{db_name}.fa -m {base_dir}/out.{db_name}.tsv {base_dir}/out.{db_name}.prepareref >> {path}/{logname}") # samma här m >>

        os.chdir(path) # go to output path

        log_parse(f'Running ariba on {db_name}', path)
        os.system(f"ariba run {base_dir}/out.{db_name}.prepareref {infile1} {infile2} out.run.{db_name} >> {path}/{logname}")

    os.system(f"mkdir {path}/Ariba_output")     # Making dir to ensure output is not in main dir
    os.system(f"mv out.* {path}/Ariba_output")  # No other names start with out. , right?
    log_parse(f'Ariba done.\n', path)
    return

def fastp_func(path, infile1, infile2, common_name):
    '''Function that takes two raw reads fastq files, one forward (1, right) and one reverse(2, left)
    and returns two trimmed fastq files as well as quality control documentation.'''

    log_parse(f'Fastp started with {infile1} and {infile2}\n', path)

    outfile1 = f'out_fastp_{common_name}_1.fq.gz'
    outfile2 = f'out_fastp_{common_name}_2.fq.gz'
    html = f'out_fastp_{common_name}.html'
    json = f'out_fastp_{common_name}.json'

    fastpinput = f'fastp -i {infile1} -I {infile2} -o {outfile1} -O {outfile2} -h {html} -j {json}'

    os.system(fastpinput)
    
    log_parse(f'Fastp complete. Four output files returned:\n{outfile1} \n{outfile2} \n{html} \n{json} \n', path)
    
    # Move to output folder
    if path != finalpath:
        os.system(f'mv {outfile1} {outfile2} {html} {json} {path}')
    
    outpath1 = f'{path}/{outfile1}'
    outpath2 = f'{path}/{outfile2}'

    return outpath1, outpath2

def kraken_func(path, infile1, infile2, threads, common_name, path_kraken):
    ''' Function that runs Kraken on two raw reads fastq files, one forward (1, right) and one reverse(2, left), 
    in order to assign taxonomic labels to the sequences'''

    log_parse(f'Kraken started with {infile1} and {infile2} as input with {threads} threads available \n', path)
    kraken_output = f'out_kraken_{common_name}.out'
    kraken_report = f'report_kraken_{common_name}.report'

    krakeninput = f'kraken2 --db {path_kraken} --threads {threads} --output {kraken_output} --report {kraken_report} --paired {infile1} {infile2}'
    
    os.system(krakeninput)

    lines = f'Kraken run finished. Two output files returned:\n'
    lines += f'{kraken_output} \n{kraken_report}\n'
    log_parse(lines, path)

    # Move to output folder
    if path != finalpath:
        os.system(f'mv {kraken_output} {kraken_report} {path}')

def reads_for_coverage(path, fastq_file, wanted_coverage, genome_size):
    '''Function that checks whether the requested coverage can be reached with the input
    files, returning the maximum coverage if this is not the case.'''

    log_parse(f'Running: reads_for_coverage', path)

    bases_needed = int(wanted_coverage*genome_size/2)
    
    log_parse(f'To achieve {wanted_coverage} X, {bases_needed} bases are needed from each fastq-file\n', path)
    log_parse(f'Checking if wanted coverage can be achieved...\n', path)
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
                log_parse(f'Coverage can be reached! It amounts to {read_counter} reads from fastq_1 which is {total_bases} bases\n\n', path)
                coverage = wanted_coverage
                break

            row_counter = row_counter%4 + 1 # makes the counter loop between 1 and 4


    # Give log output if the coverage CANNOT be achieved, and estimate new coverage
    if total_bases < bases_needed:
        log_parse(f'There are not enough bases to achieve {wanted_coverage} X coverage.\n"', path)
        available_coverage = int((2*total_bases)/genome_size)
        log_parse( f'Using an estimated coverage of {available_coverage} X instead which amounts to {read_counter} reads and {total_bases} bases from fastq_1\n\n', path)
        coverage = available_coverage

    reads_needed = read_counter
    
    log_parse( f'Function finished.\nOutputs: coverage {coverage}, reads needed {reads_needed}\n\n', path)

    return coverage, reads_needed

def shorten_fastq(path, fastq1_file, fastq2_file, reads_needed, common_name):
    '''Function that shortens the fastq files to only be long enough to reach 
    the requested coverage.'''

    log_parse(f'Shorten_fastq started.\n', path)
    log_parse(f'Shortening {fastq1_file} and {fastq2_file} to only match wanted coverage.\n\n', path)

    lines_needed = reads_needed*4
    newname1 = f'X_{common_name}_1.fastq.gz'
    newname2 = f'X_{common_name}_2.fastq.gz'

    newpath1 = f'{path}/{newname1}'
    newpath2 = f'{path}/{newname2}'
    
    with gzip.open(fastq1_file, 'rt') as trim_me: # maybe change to 'rb'
        newfile = ''
        for i, line in enumerate(trim_me):
            newfile += line
            if i == lines_needed:
                break
    
    with gzip.open(newpath1, 'wt') as one:
        one.write(newfile)

    with gzip.open(fastq2_file, 'rt') as trim_me:
        newfile = ''
        for i, line in enumerate(trim_me):
            newfile += line
            if i == lines_needed:
                break

    with gzip.open(newpath2, 'wt') as one:
        one.write(newfile)

    log_parse( f'Shortening complete.\nOutputs: {newname1}, {newname2}.\n\n', path)

    return newpath1, newpath2

def spades_func(path, file1, file2, path_spades, common_name, threads, RAM):
    '''Function that runs SPAdes to assemble contigs from short reads.'''

    log_parse('SPAdes started\n', path)

    # To make sure X_spades output is in the correct output directory. 
    # Pilon output will also be added here
    assembly_path = f'{path}/{common_name}_assembly'

    # commandline = '#SBATCH -p node -n 1 \n'
    commandline = f'python {path_spades}/spades.py --careful -o {assembly_path} --pe1-1 {file1} --pe1-2 {file2} -t {threads} -m {RAM}'
    os.system(commandline)
    #"spades.py --careful -o $filename1_short\_$wanted_coverage\X_spades --pe1-1 $read1_output --pe1-2 $read2_output -t $threads_available -m $RAM_available"

    # rename from contigs.fasta to fasta to work with pilon
    os.system(f'cp {assembly_path}/contigs.fasta {assembly_path}/{common_name}.fasta')
    log_parse( f'"contigs.fasta"-file copied and renamed to be called "{common_name}.fasta"', path)

    log_parse('SPAdes finished.\n', path)
    log_parse(f'All output files can be found here: {assembly_path}\n\n', path)

    return assembly_path

def pilon_func(path, fastafile, fasta1, fasta2, common_name, threads, assembly_path):
    '''Function that runs Pilon on contigs-file from SPAdes to 
    polish and assemble further.'''
    
    os.chdir(assembly_path) # Cannot change directory when parallel.
    
    log_parse('Pilon started', path)
    log_parse(f'Input files: {fastafile}, {fasta1}, {fasta2}', path)

    bowtie_build = f'bowtie2-build -f --threads {threads} --quiet {fastafile} {common_name}'
    os.system(bowtie_build)

    # inputs the two shortened fasta-files, if available
    bowtie = f'bowtie2 -x {common_name} -1 {fasta1} -2 {fasta2} -S {common_name}.sam --phred33 --very-sensitive-local --no-unal -p {threads}'
    os.system(bowtie)

    os.system(f'samtools view -bh {common_name}.sam > {common_name}.bam')
    os.system(f'samtools sort {common_name}.bam -o {common_name}.sorted.bam')
    os.system(f'samtools index {common_name}.sorted.bam')

    time = currenttime()+'\n'
    log_parse( f'Pilon 1.24 started at {time}', path)
    
    os.system(f'pilon --genome {common_name}.fasta --frags {common_name}.sorted.bam --output {common_name}.pilon --changes --threads {threads}')
    
    log_parse(f'Pilon finished\n', path)
    
    log_parse( f'Corrected fasta file created: {common_name}.pilon.fasta', path)

    return 

def info(path, spades_assembly):
    '''Function that uses an assembly-file from SPAdes of Pilon 
    and returns the metrics of that assembly.'''

    # Output som pandas table for att satta ihop alla strains till en sammanstalld csv-fil, 
    # och varje strain till var sin csv-fil
    os.chdir(path)

    log_parse( 'Looking at the metrics of assembly_fasta\n', path)

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
    
    log_parse( f'The number of contigs: {number_of_contigs}, the total number of bases: {total_bases}\n', path)

    contig_lengths.sort()
    longest = contig_lengths[-1]
    log_parse(  f'Longest contig: {longest}\n', path)
    log_parse(  f'Contigs longer than 1 kb: {contigs_over_1000}', path)


    # N50
    temp = 0
    for length in contig_lengths:
        temp += length
        N_50 = length
        if temp >= (total_bases/2):
            break
    
    log_parse( f'N50: {N_50}\n', path)

    # GC-content
    GC = round(number_GC*100/(number_GC + number_AT),2)
    log_parse( f'The GC-content of the sequence is {GC}%. {non_base} non-base characters were excluded from GC-calculation\n', path)

    log_parse( f'-----------------------Metrics finished-----------------------', path)

    # PLACE ALL INFO IN PANDAS TABLE
    data = {'Total nr bases': total_bases, 'Nr contigs': number_of_contigs, 'Longest contig': longest, 
    'Nr contigs > 1kb': contigs_over_1000, 'N50':N_50, 'GC-content': GC}

    info_df = pd.DataFrame(data=data, index=[0])

    return info_df

def regular(path, infile1, infile2, run_fastp, kraken, ariba, db_ariba, run_spades, wanted_coverage, genome_size, pilon, threads, RAM):
    '''Function that runs the regular pipeline. This function is called from the parallelize
    function in the case of calling the pipeline with a directory of multiple reads files.
    Requires an output path, a forward file (1), a reverse file (2) as well as other predetermined 
    parameters'''

    time = currenttime()
    date = str(datetime.date(datetime.now()))

#                   |  |
#                   |  |
#                 __|  |__
#                 \      /
#                  \    /
#                   \  /   
#                    \/
################################################
# CHANGE THESE PATHS TO FIT YOUR DOWNLOADS! :D
################################################

    path_tools = '/proj/uppmax2022-2-14/private/campy_pipeline/assembly/verktyg'
    path_spades = path_tools + '/SPAdes-3.15.4-Linux/bin'
    path_kraken = path_tools + '/minikraken2_v1_8GB'

################################################
# NO MORE MODIFICATIONS NEEDED (HOPEFULLY)
################################################

    os.system(f'cp {infile1} {infile2} {path}')

# if path for infiles has been sent in, then shorten the names. Otherwise there will be no change.
    
    f1 = infile1.split('/')[-1]
    f2 = infile2.split('/')[-1]

    common_name = shortname(f1)

    infile1 = f'{path}/{f1}'
    infile2 = f'{path}/{f2}'

# Create log file
    create_log(path, time, date, infile1, infile2)
    
# Ariba 
    if ariba:
        header= '\n'+'='*15 +'ARIBA'+ '='*15 +'\n'
        log_parse(header, path)
        ariba_fun(path, infile1,infile2, db_ariba)
        #os.system("ariba summary out_sum out.run.*/report.tsv") #change from v5

# Fastp
    if run_fastp:
        header= '\n'+'='*15 +'FASTP'+ '='*15 +'\n'
        log_parse(header, path)
        outfile1_trim, outfile2_trim = fastp_func(path, infile1, infile2, common_name)

        infile1 = outfile1_trim
        infile2 = outfile2_trim

# Kraken
    if kraken:
        header= '\n'+'='*15 +'KRAKEN'+ '='*15 +'\n'
        log_parse(header, path)
        kraken_func(path, infile1, infile2, threads, common_name, path_kraken)

# Number of reads to match the wanted coverage
    if run_spades:
        header= '\n'+'='*15 +'READS FOR COVERAGE, SPADES'+ '='*15 +'\n'
        log_parse(header, path)
        coverage, reads_needed = reads_for_coverage(path, infile1, wanted_coverage, genome_size)
    else:
        coverage = 0

# Shortening fastq-files if the coverage can be reached with less
    if coverage > wanted_coverage:
        outfile1_shorten, outfile2_shorten = shorten_fastq(path, infile1, infile2, reads_needed, common_name)        
        infile1 = outfile1_shorten
        infile2 = outfile2_shorten

# Spades
    if run_spades:
        header= '\n'+'='*15 +'SPADES'+ '='*15 +'\n'
        log_parse(header, path)
        assembly_path = spades_func(path, infile1, infile2, path_spades, common_name, threads, RAM)

# Pilon
    # if pilon:
    #     header= '\n'+'='*15 +'PILON'+ '='*15 +'\n'
    #     log_parse(header, path)
    #     fastafile = f'{assembly_path}/{common_name}.fasta'
    #     pilon_func(path, fastafile, infile1, infile2, common_name, threads, assembly_path)
        
# Info/metrics
    if run_spades:
        header= '\n'+'='*15 +'INFO/METRICS, SPADES'+ '='*15 +'\n'
        log_parse(header, path)

        from_spades = f'{assembly_path}/{common_name}.fasta'
        
        info_df = info(path, from_spades)
        info_df.to_csv(header = True, path_or_buf = f'{path}/{common_name}_info.csv')

    return 1
        
def map_func(dir, f):
    '''Function to map regular to files and directory when running in parallel'''
    return regular(dir, f[0], f[1], run_fastp, kraken, ariba, db_ariba, run_spades, wanted_coverage, genome_size, pilon, threads, RAM)

def parallelize(finalpath, file_directory):
    '''Function that takes a directory of forward and reverse files to run the 
    pipeline with in parallel. Also takes a final path to the collective directory'''
    

    os.chdir(file_directory)
    os.system(f"ls *.gz > input.txt")
    # go back
    os.chdir(base_dir)

    dirlist = []
    files = []
    with open(f'{file_directory}/input.txt', 'r') as inp:
        linelist = inp.readlines()
        for i in range(0, len(linelist), 2):
            common_name = shortname(linelist[i])
            path = f'{finalpath}/{common_name}'
            os.mkdir(path)
            dirlist.append(path)
            f1 = linelist[i].strip('\n')
            f2 = linelist[i+1].strip('\n')
            files.append((f'{file_directory}/{f1}', f'{file_directory}/{f2}'))
    
    with future.ThreadPoolExecutor() as ex:
        results = ex.map(map_func, dirlist, files)

    os.chdir(finalpath)
 
    # Creating combined info-files for parallellized genomes, currently names are last but works.
    finalname="sum_info"
    all_filenames = [i for i in glob.glob(f'{finalpath}/*/*info.csv')]
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ], axis=0)
    just_names = [name.split('/')[-1].split('_')[0] for name in all_filenames]
    combined_csv["Genome Name"] = just_names
    combined_csv.to_csv(path_or_buf= f'{finalpath}/{finalname}.csv', index=False, encoding='utf-8-sig')
    
def main():
    """
    path/to/file1 path/to/file2 here nopar notrim nokraken ariba [db1, db2] 0 size nopilon thr ram
    """
    global new_location, run_fastp, kraken, ariba, db_ariba, wanted_coverage, genome_size, pilon, threads, RAM, run_spades, finalpath, logname, base_dir 
    infile1 = sys.argv[1] # 
    infile2 = sys.argv[2]
    here = sys.argv[3].lower() == 'here' # will ask for directory location if False
    new_location = sys.argv[3]
    run_fastp = sys.argv[4].lower() == 'trim' # will run fastp if True
    kraken = sys.argv[5].lower() == 'kraken'
    ariba = sys.argv[6].lower() == 'ariba'
    db_ariba = sys.argv[7][1:-1].strip(" ").split(',')
    wanted_coverage = int(sys.argv[8]) # if wanted coverage == 0, then don't run spades
    genome_size = int(sys.argv[9])
    pilon = sys.argv[10].lower() == 'pilon'
    threads = sys.argv[11]
    RAM = sys.argv[12]

    run_spades = wanted_coverage != 0

    if pilon and run_spades == False: # Since pilon requires spades output, this 
        pilon = False

# Let's start this pipeline!

    logname = 'logfile.txt'
    time = currenttime()
    date = str(datetime.date(datetime.now()))
    base_dir = os.getcwd()

# make directory for output
    finalpath = directory(date, time, here)

    if os.path.isdir(infile1):
        parallelize(finalpath, infile1)
        #os.system("ariba summary out.sum out.run.*/report.tsv")

    else:
        regular(finalpath, infile1, infile2, run_fastp, kraken, ariba, db_ariba, run_spades, wanted_coverage, genome_size, pilon, threads, RAM)


if __name__ == '__main__':
    main()  
