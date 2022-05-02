import sys
import os
from datetime import datetime
import gzip
import re
import csv
import pandas as pd


# Start by parsing the following command through the terminal, choosing only one option in each case:
# 'python assembly_pipeline_v2.py infile1/folder(???) infile2/none(???) here/there regular/parallel trim/notrim kraken/nokraken ariba/noariba wanted_coverage genome_size pilon/nopilon threads RAM'
# go-to:
# python assembly_pipeline_v2.py SRR18825428_1.fastq.gz SRR18825428_2.fastq.gz here regular trim kraken noariba 0 0 nopilon

'''OPTIONS'''
# infile1
# - here/there: Where should all outputs be saved? If 'here' a new directory is created in 
# the current directory. If 'there' a path will be asked for.
# - regular/parallel: regular means running only one strain, parallel means running multiple strains
# - trim/notrim: trim means we run fastp, notrim means that we don't
# - kraken/nokraken
# - ariba/noariba
# - wanted_coverage
# - genome_size
# - pilon/nopilon
# - threads
# - RAM


''' Function to create directory where all outputs from the pipeline are placed. 
Date and time specific'''
def directory(date, time, there = False):

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
    time = str(datetime.time(datetime.now()))
    return time

def shortname(filename):
    short = re.search('[a-zA-Z1-9]+', filename).group()
    return short

def fastp_func(infile1, infile2, common_name):

    loglines = f'Fastp started with {infile1} and {infile2}\n'

    outfile1 = f'out_fastp_{common_name}_1.fq.gz'
    outfile2 = f'out_fastp_{common_name}_2.fq.gz'

    fastpinput = 'fastp -i ' + infile1 + ' -I ' + infile2 + ' -o ' + outfile1 + ' -O ' + outfile2

    os.system(fastpinput)

    loglines += 'Fastp complete. Four output files returned:\n'
    loglines += f'{outfile1} \n{outfile2} \nfastp.html \nfastp.json \n\n'
    return outfile1, outfile2, loglines

def reads_for_coverage(wanted_coverage, genome_size):
    bases_needed = int(wanted_coverage*genome_size/2)
    
    loglines = f'To achieve {wanted_coverage} X, {bases_needed} bases are needed from each fastq-file\n'


    return loglines

def spades_func(file1, file2, path_spades, common_name, finalpath):

    loglines = 'SPAdes started\n'
    # To make sure X_spades output is in the correct output directory
    assembly_path = f'{finalpath}/{common_name}_spades'
    # commandline = '#SBATCH -p node -n 1 \n'
    commandline += f'python {path_spades}/spades.py --careful -o {assembly_path} --pe1-1 {file1} --pe1-2 {file2}'
    os.system(commandline)
    #"spades.py --careful -o $filename1_short\_$wanted_coverage\X_spades --pe1-1 $read1_output --pe1-2 $read2_output -t $threads_available -m $RAM_available"

    # rename from contigs.fasta to fasta
    os.system(f'cp {assembly_path}/contigs.fasta {assembly_path}/{common_name}.fasta')
    loglines += f'"contigs.fasta"-file copied and renamed to be called "{common_name}.fasta"'

    loglines += 'SPAdes finished.\n'
    loglines += f'All output files can be found here: {assembly_path}\n\n'

    return loglines

def pilon_func():
    pass

def info(assembly_fasta):
    # Output som pandas table for att satta ihop alla strains till en sammanstalld csv-fil, 
    # och varje strain till var sin csv-fil

    loglines = 'Looking at the metrics of assembly_fasta'

    number_of_contigs, bases_in_contig, total_consensus = 0, 0, 0
	#@contig_lengths = ""
    N50, non_base, number_AT, number_GC = 0, 0, 0, 0
	#@GC = ""
	#@sequence = "";
    a, b = 0, 0
	#@sorted_contig_lengths = ""
    temp, N50_temp_var, contigs_over_1000, total_number_bases = 0, 0, 0, 0

    assembly = open(assembly_fasta, 'r')
    for line in assembly:
        if '>' in line:
            # bla...
            pass

    return loglines


# function that runs multiple strains in parallel. Inputs are all sys.argv[]
# Return lines for logfile?
def parallelize():
    pass

# function that runs everything for only one strain. Inputs are all sys.argv[]
# Return lines for logfile?
def regular():
    pass


def main():
    infile1 = sys.argv[1]
    infile2 = sys.argv[2]
    common_name = shortname(infile1) # until here only work if not parallel
    new_location = sys.argv[3] == 'there' # will ask for directory location if True
    parallel = sys.argv[4] == 'parallel' # TO BE MODIFIED to run strains in parallel if True
    run_fastp = sys.argv[5] == 'trim' # will run fastp if True
    kraken = sys.argv[6] == 'kraken'
    ariba = sys.argv[7] == 'ariba'
    spades = sys.argv[8] != '0' # if wanted coverage == 0, then don't run spades
    genome_size = sys.argv[9]
    pilon = sys.argv[10] == 'pilon'
    threads = sys.argv[11]
    RAM = sys.argv[12]




    

    path_tools = '/proj/uppmax2022-2-14/private/campy_pipeline/assembly/verktyg'
    path_spades = path_tools + '/SPAdes-3.15.4-Linux/bin'
    path_kraken = path_tools + '/minikraken2_v1_8GB'
    
    time = currenttime()
    date = str(datetime.date(datetime.now()))
    
    # make directory for output
    finalpath = directory(date, time, new_location)
    

# Create log file
    logname = 'logfile'
    log = open(logname, 'w')

    lines = 15*'-' + 'LOGFILE' + 15*'-' + '\n\n'

    lines += f'New directory created with the adress {finalpath}\n'
    lines += f'Directory created at {time} on {date}\n'
    lines += 'All outputs will be saved in the new directory.\n\n'
    log.writelines(lines)

# Parallel
    # if parallel:
    # parallelize()

# Fastp
    if run_fastp:
        time = currenttime()+'\n'
        log.writelines(time)

        # infile1 = input('Give the fastq, gzipped forward file you want for fastp: ')
        # infile2 = input('Give the fastq, gzipped reverse file you want for fastp: ')

        outfile1, outfile2, fastp_lines = fastp_func(infile1, infile2, common_name)
        log.writelines(fastp_lines)

        infile1 = outfile1
        infile2 = outfile2

# Kraken
    if kraken:

        #chdir "$path_assembly";

        time = currenttime()+'\n'
        log.writelines(time)
            
        if path_kraken: # If path_kraken exists do thing, maybe it should be os.path.exists('kraken_db') instead? Does this even has to be an if statement at all?
        
            loglines = f'Kraken started with {infile1} and {infile2}\n' #unclear if "filename" is right
            #system ("kraken2 --db $kraken_DB --threads $threads_available --output $filename1_short.kraken.out --report $filename1_short.kraken.report --paired $fastqfile1 $fastqfile2");
            # Add later lol --threads {threads_available}
            # in the old pipeline they have two names for the input (fastqfile1 and filename1_short) and I do not know why
            #krakeninput = f'kraken2 --db {path_kraken} --output output_test.kraken.out --report output_test.kraken.report --paired {infile1} {infile2}'
            krakeninput = f'kraken2 --db {path_kraken} --output out_kraken.out --report report_kraken.report --paired {infile1} {infile2}'
            os.system(krakeninput)

            time = currenttime()
            log.writelines(time)

            #csv_kraken = r"kraken_report.csv"

            with open('report_kraken.report', 'rb') as kraken_report: # add "or die" or similar? see old pipeline
            #    result_reader = csv.reader(kraken_report, delimiter=' ')
            #   # result_reader.next()
            #    for row in result_reader:
            #        for (i,v) in enumerate(row):
            #            columns[i].append(v)
            #print(columns[0])

            #with open(txt_file, "r") as in_text:
                out_filename= 'kraken_csv.csv'
                
                df = pd.read_csv(kraken_report, sep="/t")
                df.to_csv(out_filename, index=False)

                df = pd.read_csv(out_filename, header=None)
                print(df.shape)
                #Due to python indexing, the enumerating of the columns starts at 0
                short_report = df.iloc(:, [0, 2, 5])
                print(short_report.shape)
                #column1 = df.iloc[:,0]
                #column3 = df.iloc[:,2]
                #column6 = df.iloc[:,5]

                #short_report = pd.concat([column1, column3, column6], axis=1)

                print(short_report)

            #    kraken_report.readlines()

# Number of reads to match the wanted coverage
    if spades:
        pass
# Spades
    if spades:
        time = currenttime()+'\n'
        log.writelines(time)

        spades_lines = spades_func(infile1, infile2, path_spades, common_name, finalpath)
        log.writelines(spades_lines)

# Pilon
    time = currenttime()
    log.writelines(time)


# Move files to correct folder
    os.system('mv ' + outfile1 + ' ' + outfile2 + ' fastp.html fastp.json ' + str(finalpath))
    log.writelines('Fastp output files moved to directory\n\n')

# Close and move the logfile to the correct directory
    log.close()
    os.system('mv logfile '+str(finalpath))
    

if __name__ == '__main__':
    main()  
