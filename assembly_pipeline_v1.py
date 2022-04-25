import sys
import os
from datetime import datetime
import gzip

#import fastp

# pipeline tar in short reads
# tar in huruvida vi vill assembla


#fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

#os.system("ls")

# OBS ge fastq-info sa den kan sattas in i pandas-tabellen (info) ocksa

''' Function to create directory where all outputs from the pipeline are placed. 
Date and time specific'''
def directory(date, time):

    # Change the format of time from eg. 09:13:07.186006 to 09h13m07s
    stringtime = time[:2]+'h'+time[3:5]+'m'+time[6:8]+'s'

    # Choose path of directory
    print('The default path to all output files is a new directory in your current path.')
    print('Do you want to save all files in another directory?')
    ans = input('[y/n]: ')
    while ans not in ['y', 'n']:
        ans = input(['Please write y or n: [y/n]: '])
    if ans == 'y':
        path = input('New path: ')
    elif ans != 'y':
        path = os.getcwd()

    # Rename directory with date and time
    namedir = 'assembly_' + date + '_' + stringtime

    finalpath = os.path.join(path, namedir)

    os.mkdir(finalpath)
    
    return finalpath

def fastp_func(infile1, infile2):
    # format = input('Interleaved or two files?: [i/2]')
    # if format == 'i':


    loglines = f'Fastp started with {infile1} and {infile2}'

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


    loglines += 'Fastp complete\n\n'
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


def main():
    time = str(datetime.time(datetime.now()))
    date = str(datetime.date(datetime.now()))
    
    finalpath = directory(date, time)

    # Create log file
    filename = 'logfile'
    log = open(filename, 'w')

    lines = 'LOGFILE' + '\n\n'
    lines += 'New directory created with the adress ' + finalpath +'\n'
    lines += 'Directory created at '+ time+ ' on '+ date
    lines += 'All outputs will be saved in the new directory.'
    log.writelines(lines)

# Fastp
    infile1 = input('Give the fastq, gzipped forward file you want for fastp: ')
    infile2 = input('Give the fastq, gzipped reverse file you want for fastp: ')
    outfile1, outfile2, loglines = fastp_func(infile1, infile2)
    log.writelines(loglines)
    os.system('mv ' + outfile1 + ' ' + outfile2 + ' fastp.html fastp.json ' + str(finalpath))
    log.writelines('Fastp output files moved to directory')
    # os.system('mv' + outfile2 +str(finalpath))
    # os.system('mv fastp.html' +str(finalpath))
    # os.system('mv fastp.json' +str(finalpath))



# Close and move the logfile to the correct directory
    log.close()
    os.system('mv logfile '+str(finalpath))

if __name__ == '__main__':
    main()  
