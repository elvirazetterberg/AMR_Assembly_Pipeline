import sys
import os
from datetime import datetime

#import fastp

# pipeline tar in short reads
# tar in huruvida vi vill assembla


#fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

#os.system("ls")

# OBS ge fastq-info så den kan sättas in i pandas-tabellen (info) också

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
    else:
        path = os.getcwd()

    # Rename directory with date and time
    namedir = 'assembly_' + date + '_' + stringtime

    finalpath = os.path.join(path, namedir)

    os.mkdir(finalpath)
    
    return finalpath

def fastp_func(rawreads):
    loglines = 'Fastp started'

    read = ''
    for i, line in enumerate(rawreads):
        if line[0] == '@' and i != 0:
            read1 = read
            read = ''
        read += line
    
    read2 = read

    var = os.system(f'fastp -i {read1} -I {read2} -o out.R1.fq.gz -O out.R2.fq.gz')

    loglines += 'Fastp complete'
    return var, loglines

def info(assembly_fasta):
    # Output som pandas table för att sätta ihop alla strains till en sammanställd csv-fil, 
    # och varje strain till var sin csv-fil

    loglines = f'Looking at the metrics of {assembly_fasta}'

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
    lines += f'Directory created at {time} on {date}'
    lines += 'All outputs will be saved in the new directory.'
    log.writelines(lines)

    rawreads = input('Give the fastq, gzipped file you want for fastp: ')
    output, loglines = fastp_func(rawreads)
    


# Close and move the logfile to the correct directory
    log.close()
    os.system(f"mv logfile {finalpath}")

if __name__ == '__main__':
    main()  