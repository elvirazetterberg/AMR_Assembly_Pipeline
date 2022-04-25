import sys
import os
from datetime import datetime
import gzip

# Parsing through the terminal:
# Start with python assembly_pipeline_v1.py here/there regular/parallel trim/notrim kraken/nokraken ariba/noariba pilon/nopilon
# - here/there: Where should all outputs be saved? If 'here' a new directory is created in 
# the current directory. If 'there' a path will be asked for.
# - regular/parallel: regular means running only one strain, parallel means running multiple strains
# - trim/notrim: trim means we run fastp, notrim means that we don't

# pipeline tar in short reads
# tar in huruvida vi vill assembla


#fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

#os.system("ls")

# OBS ge fastq-info sa den kan sattas in i pandas-tabellen (info) ocksa

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

def fastp_func(infile1, infile2):

    loglines = f'Fastp started with {infile1} and {infile2}\n'

    ''' for the interleaved file '''
    # with gzip.open(rawreads, 'rt') as r:
    #     read = ''
    #     for i, line in enumerate(r):
    #         if line[0] == '@' and i != 0:
    #             with gzip.open('read1.fq.gz', 'wt') as one:
    #             # read1 = open('read1.fq.gz', 'w')
    #                 one.write(read)
    #             read = ''
    #         read += line
    
        # # read2 = open('read2.fq.gz', 'w')
        # with gzip.open('read2.fq.gz', 'wt') as two:
        #     two.writelines(read)

    outfile1 = 'out.R1.fq.gz'
    outfile2 = 'out.R2.fq.gz'

    fastpinput = 'fastp -i ' + infile1 + ' -I ' + infile2 + ' -o ' + outfile1 + ' -O ' + outfile2

    os.system(fastpinput)

    loglines += 'Fastp complete. Four output files returned:\n'
    loglines += f'{outfile1} \n{outfile2} \nfastp.html \nfastp.json \n'
    return outfile1, outfile2, loglines

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
    time = str(datetime.time(datetime.now()))
    date = str(datetime.date(datetime.now()))
    
    if sys.argv[1] == 'there':
        there = True
        finalpath = directory(date, time, there)
    else:
        finalpath = directory(date, time)

# Create log file
    filename = 'logfile'
    log = open(filename, 'w')

    lines = 15*'-' + 'LOGFILE' + 15*'-' + '\n\n'

    lines += f'New directory created with the adress {finalpath}\n'
    lines += f'Directory created at {time} on {date}\n'
    lines += 'All outputs will be saved in the new directory.\n\n'
    log.writelines(lines)

# Parallel
    # if sys.argv[1] == 'parallel':
    #     parallelize()
    # else:
    #     continue

# Fastp
    infile1 = input('Give the fastq, gzipped forward file you want for fastp: ')
    infile2 = input('Give the fastq, gzipped reverse file you want for fastp: ')

    outfile1, outfile2, loglines = fastp_func(infile1, infile2)
    log.writelines(loglines)

    os.system('mv ' + outfile1 + ' ' + outfile2 + ' fastp.html fastp.json ' + str(finalpath))
    log.writelines('Fastp output files moved to directory\n')



# Close and move the logfile to the correct directory
    log.close()
    os.system('mv logfile '+str(finalpath))

if __name__ == '__main__':
    main()  
