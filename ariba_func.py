import sys
import os
from datetime import datetime
import gzip
import re

# Start by parsing the following command through the terminal, choosing only one option in each case:
# 'python assembly_pipeline_v2.py infile1/folder(???) infile2/none(???) here/there regular/parallel trim/notrim kraken/nokraken ariba/noariba wanted_coverage genome_size pilon/nopilon threads RAM'

# test run:
# python assembly_pipeline_v2.py SRR18825428_1.fastq.gz SRR18825428_2.fastq.gz here regular trim kraken ariba 40 124000000 pilon 40 0

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


def ariba(infile1,infile2,db_ariba):
    for db_name in db_ariba[1:-1].split(','): #klumpigt? as sysargv makes input a string, it is separated into a list here.
        
        #OBS when making parallell the naming of files must take this into account. Right now Im deleting the privious runs

        #if there's allready an existing db, lite fult gjort?
        # In the pearlpipeline it says "$test_CARD = "CARD_reference_dataset_downloaded"; What does it mean? 

        if os.path.exists(f'out.{db_name}.fa'): 
            loglines=f'{currenttime()}: Database {db_name} allready downloaded\n'
            os.system(f"rm -rf out.run.{db_name}") # är detta smart sätt att göra det?

        else: # if database not downloaded. troligen onädigt i framtiden
            rm_db=input(f'Please note that running this will remove all existing files starting with "out.{db_name}". [y]/n?') 
            if rm_db.lower().startswith("y")==False:
                
                exit() 
            else:
                os.system(f"rm -rf out.{db_name}*")

            loglines= f'{currenttime()} Downloading database {db_name}\n'
            os.system(f"ariba getref {db_name} out.{db_name}")

            loglines+=f'{currenttime()} Preparing references with prefix out.{db_name} \n'
            os.system(f"ariba prepareref -f out.{db_name}.fa -m out.{db_name}.tsv out.{db_name}.prepareref")

        loglines +=f'{currenttime()}: Running ariba on {db_name}\n'
        os.system(f"ariba run out.{db_name}.prepareref {infile1} {infile2} out.run.{db_name}")

    loglines+=f'{currenttime()}: Ariba done.\n'
    return loglines
#-------------------------------------------------------------------------------------------------------
 #Pipeline
#-------------------------------------------------------------------------------------------------------

#call : python3 ariba_test.py /home/alma/Documents/kandidat/genomes/SRR18825428_1.fastq /home/alma/Documents/kandidat/genomes/SRR18825428_2.fastq ariba [vfdb_core]
def main():
    infile1 = sys.argv[1]
    infile2 = sys.argv[2]
    ariba = sys.argv[3] == 'ariba'
    db_ariba =sys.argv[4]
    new_location =False # will ask for directory location if True


    # Hardcoded, location of non-conda tools
    #path_tools = '/proj/uppmax2022-2-14/private/campy_pipeline/assembly/verktyg'
    #path_spades = path_tools + '/SPAdes-3.15.4-Linux/bin'

    # Let's start this pipeline!
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

#Ariba----------------------------------------------------------------------
#("Available databases: argannot, card, ncbi, megares, plasmidfinder, resfinder, srst2_argannot, vfdb_core, vfdb_full, virulencefinder")

    if ariba:
        log.writelines(f'=================================\n{currenttime()}: Ariba\n=================================\n')
        log.writelines(ariba(infile1,infile2, db_ariba))



if __name__ == '__main__':
    main()  

