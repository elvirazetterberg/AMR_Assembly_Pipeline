from datetime import datetime
import pandas as pd
import os
import concurrent.futures as future
import subprocess as sp

"""

                    HI  AND  WELCOME  TO  THE 
     ____  ____  _____    ____  _ ____  _____  _     _ __   _  _____  
    |  _ \|  _ \| ____|  |  _ \| |  _ \|  ___|| |   | |  \ | ||  ___| 
    | |_) | |_) | |__    | |_) | | |_) | |__  | |   | |   \| || |___  
    |  __/|    /| |___   |  __/| |  __/| |___ | |___| | |\ \ || |___  
    |_|   |_|\_\|_____|  |_|   |_|_|   |_____||_____|_| | \__||_____| 

                        COMMAND GENERATOR!

This script generates a command in the following format:

'python assembly_pipeline_v8.py infile1/directory infile2/None here/there trim/notrim kraken/nokraken ariba/noariba wanted_coverage genome_size pilon/nopilon threads'

And will start the pipeline with the generated command if wanted.


Reference:

"""


def command_generator(frst,scnd,one_or_sev):

    #Check path for input and check if files exist:
    path=input(f"If not in current directory, please input absolute path for {one_or_sev}. Otherwise press ENTER: ") or sp.getoutput("pwd") #os.curdir not working?
    print(path)
    file_ex = sp.getoutput(f' [ -f {path}/{frst} ] && echo "File exist" || echo "File does not exist"')
    if file_ex!="File exist":
        exit(f'Could not find "{frst}". Exiting.')
    # Check if output should be in current directory
    pos=input("If you want to change the path for the output directory, please write 'there'. Otherwise press ENTER: ") or "here"

    # Check what tools are wanted, default is nothing
    for tool in tools:
        t=input(f"Do you want to run {tool}? y/[n]:\t") or "no"
        if t.lower().startswith("y"):
            chosen_tools.append(tool.lower()) 
        else:
            chosen_tools.append(f"no{tool}")
    if "fastp" in chosen_tools: # special case: fastq is called using the word "trim"
        chosen_db[0]="trim"

    # If ariba, Check what databases are wanted. OBS ugly
    no="cont"
    while no!="stop":
        no="stop"
        if "ariba" in chosen_tools:
            db=list(input("On what databases do you wish to run ariba? \t").strip("[] ").split(","))
            for d in db:
                if d not in databases:
                    print(f"One or more databases not recognized. Please choose one of the following:\n {databases}")
                    no="cont"
            print(f"Chosen databases: {chosen_db}")
    
    # Checking parameters for assebly, if none given default is 0
    cov=input("What coverage is requested? If 0, no assembly is performed.\t") or "0"
    gen_siz=input("What is the genome size of your organism?\t") or "0"
    if cov!=0:
        pilon=input("Do you wish to run pilon? y/[n]\t")+"pilon" or "no"
        if pilon.lower().startswith("y"): # not working for some reason
            pilon="pilon"
    threads=input("Please write the maximum threads available:\t") or "0"

    # Compiling the info into the command, one for the case of several genomes and one for a singular gneome
    if one_or_sev=="directory":
        command=f"python assembly_pipeline_v8.py {path}/{frst} {scnd} {pos} {chosen_tools[0]}"+\
            f"{chosen_tools[1]} {chosen_tools[2]} {chosen_db} {cov} {gen_siz} {pilon} {threads}"
    else:
        command=f"python assembly_pipeline_v8.py {path}/{frst}_1.fastq.gz {path}/{scnd}_2.fastq.gz {pos} {chosen_tools[0]}"+\
            f"{chosen_tools[1]} {chosen_tools[2]} {chosen_db} {cov} {gen_siz} {pilon} {threads}" #if they enter the filename this is innacurate
    return  command

def main():  
    
    # print(f" \n\n\
    #  ____  ____   _____    ____  _ ____  _____  _     _ __   _  _____  \n\
    # |  _ \|    \ | ____|  |  _ \| |  _ \|  ___|| |   | |  \ | ||  ___| \n\
    # | |_) |    / | |__    | |_) | | |_) | |__  | |   | |   \| || |___  \n\
    # |  __/| |\ \ | |___   |  __/| |  __/| |___ | |___| | |\ \ || |___  \n\
    # |_|   |_| \_\|_____|  |_|   |_|_|   |_____||_____|_| | \__||_____| \n\n\
    # COMMAND  GENERATOR  FOR  ASSEMBLY, & AMR  &  VIRULENCE  DETECTION\n\n")

    print(f" \n\n\
     ____  ____  _____    ____  _ ____  _____  _     _ __   _  _____  \n\
    |  _ \|    \| ____|  |  _ \| |  _ \|  ___|| |   | |  \ | ||  ___| \n\
    | |_) | |_) | |__    | |_) | | |_) | |__  | |   | |   \| || |___  \n\
    |  __/|    /| |___   |  __/| |  __/| |___ | |___| | |\ \ || |___  \n\
    |_|   |_|\_\|_____|  |_|   |_|_|   |_____||_____|_| | \__||_____| \n\n\
        \
    COMMAND  GENERATOR  FOR  ASSEMBLY, & AMR  &  VIRULENCE  DETECTION\n\n")

    # VARIABLES:
    global tools, chosen_tools, databases, chosen_db

    tools=["fastq","Kraken", "Ariba"]
    chosen_tools=[]

    databases=["vfdb_core","vfdb_full", "ncbi", "megares","plasmidfinder","virulencefinder","card", "srst2_argannot","argannot","resfinder"]
    chosen_db=[]

    # START OF COMMAND GENERATION:
    sev=input("Do you wish to run several genomes? [y]/n:\t") or "yes" 
    if sev.lower().startswith("y"):
        directory=input("Please input name of directory containing genomes: ")
        command=command_generator(directory,"None","directory")
    else:
        name=input("please enter the name (ID) of your genome: ") or exit("No name given. Exiting.") #should we ask for name of file or name? 

        command=command_generator(name,name,"genome")


    a=input(f"Your command is: \n {command}\nDo you wish to run the pipeline imemdiately? y/[n]") or "no"
    if a.lower is "y":
        os.system(command)
    else:
        print("Command generation finished. Bye!")
    
if __name__ == '__main__':
    main()  