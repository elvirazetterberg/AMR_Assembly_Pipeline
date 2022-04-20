import sys
import os
from datetime import datetime
#import fastp

# pipeline tar in short reads
# tar in huruvida vi vill assembla


#fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

#os.system("ls")

def directory():

    date = str(datetime.date(datetime.now()))

    print('The default path to all output files is a new directory in your current path.')
    ans = input('Do you want to save all files in another directory? [y/n]')
    while ans not in ['y', 'n']:
        ans = input(['Please write y or n: [y/n]: '])
    if ans == 'y':
        path = input('New path: ')
    else:
        path = os.getcwd()

    newdir = 'assembly_' + date

    finalpath = os.path.join(path, newdir)

    os.mkdir(finalpath)
    
    return finalpath





def __main__():
    finalpath = directory()

    filename = 'logfile'

    log = open(filename, 'w')

    lines = 'LOGFILE' + '\n\n'
    lines += 'New directory created:' + finalpath

    log.writelines(lines)







# Close and move the logfile to the correct directory
    log.close()
    os.system(f"mv logfile {finalpath}")
