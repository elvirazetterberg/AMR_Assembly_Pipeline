import os
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

'python AMR_Assembly_Pipeline.py infile1/directory infile2/None here/path trim/notrim kraken/nokraken ariba/noariba [databases] wanted_coverage genome_size pilon/nopilon threads RAM'

And will start the pipeline with the generated command if wanted.

REFERENCE:
The work of third year students at the Master's Programme in Molecular Biotechnology Engineering (X) at Uppsala University 2022. 
GitHub: https://github.com/elvirazetterberg/AMR_Assembly_Pipeline.git
Contributors: Alma Nilsson, Corinne Olivero, Elvira Zetterberg, Evelina Andersson, Julia Sulyaeva, Moa Qvarnlöf

"""

def while_yn(string, default_ans):
    ans=""
    yes_or_no = False
    while yes_or_no==False:
        ans=input(string) or default_ans
        yes_or_no = ans.lower().startswith("y") or ans.lower().startswith("n")
        if yes_or_no==False:
            print(f"Please answer 'yes' or 'no'. Giver answer '{ans}' does not satisfy the criteria.")
    return ans

def while_not_nr(string, default_ans):
    ans=""
    is_nr = False
    while is_nr==False:
        ans=input(string) or default_ans
        is_nr = ans.isnumeric()
        if is_nr==False:
            print(f"Giver answer '{ans}' is not numeric. Please try again.")
    return ans

def command_generator(frst,scnd,one_or_sev):

    #Check path for input and check if files exist:
    path=input(f"If not in current directory, please input absolute path for {one_or_sev}. Otherwise press ENTER: ") or sp.getoutput("pwd") 
    print("Current directory:",path)

    file_ex = sp.getoutput(f' [ {path}/{frst} ] && echo "File exist" || echo "File does not exist"')
    if file_ex!="File exist":
        exit(f'Could not find "{frst}". Exiting.')
    # Check if output should be in current directory
    pos=input("If you want to change the path for the output directory, please enter output path. Otherwise press ENTER: ") or "here"

    # Check what tools are wanted, default is nothing
    for tool in tools:
        t=while_yn(f"Do you want to run {tool}? y/[n]:\t", "no")
        if t.lower().startswith("y"):
            chosen_tools.append(tool.lower()) 
        else:
            chosen_tools.append(f"no{tool.lower()}")
    if "fastp" in chosen_tools: # special case: fastq is called using the word "trim"
        chosen_tools[0]="trim"
        

    global chosen_db    
    # If ariba, Check what databases are wanted
    if "ariba" in chosen_tools: 
        print("On what databases do you wish to run Ariba?")
        for db in databases:
            t=while_yn(f"Do you want to run {db}? y/[n]:\t", "no")
            if t.lower().startswith("y"):
                chosen_db+=db.lower()+','
        chosen_db=chosen_db[:-1]
    chosen_db+=']'

    # Checking parameters for assebly, if none given default is 0
    cov=input
    cov=while_not_nr("What coverage is requested? If 0, no assembly is performed.\t", "0")
    gen_siz=while_not_nr("What is the genome size of your organism?\t", "0")
    if cov!=0:
        pilon=while_yn("Do you wish to run pilon? y/[n]\t", "no")+ "pilon"
        if pilon.lower().startswith("y"):
            pilon="pilon"
    threads=while_not_nr("Please write the maximum threads available:\t", "0")

    # Compiling the info into the command, one for the case of several genomes and one for a singular gneome
    if one_or_sev=="directory": 
        command=f"AMR_Assembly_Pipeline.py {path}/{frst} {scnd} {pos} {chosen_tools[0]} "+\
            f"{chosen_tools[1]} {chosen_tools[2]} {chosen_db} {cov} {gen_siz} {pilon} {threads}"
    else:
        command=f"AMR_Assembly_Pipeline.py {path}/{frst} {path}/{scnd} {pos} {chosen_tools[0]} "+\
            f"{chosen_tools[1]} {chosen_tools[2]} {chosen_db} {cov} {gen_siz} {pilon} {threads}" 

    command += " " + while_not_nr("What amount of RAM is needed?","0")
    return  "python " + command

def main():  
    
    print(f" \n\n\
     ____  ____  _____    ____  _ ____  _____  _     _ __   _  _____  \n\
    |  _ \|  _ \| ____|  |  _ \| |  _ \|  ___|| |   | |  \ | ||  ___| \n\
    | |_) | |_) | |__    | |_) | | |_) | |__  | |   | |   \| || |___  \n\
    |  __/|    /| |___   |  __/| |  __/| |___ | |___| | |\ \ || |___  \n\
    |_|   |_|\_\|_____|  |_|   |_|_|   |_____||_____|_| | \__||_____| \n\n\
        \
    COMMAND  GENERATOR  FOR  ASSEMBLY, & AMR  &  VIRULENCE  DETECTION\n\n")

    # VARIABLES:
    global tools, chosen_tools, databases, chosen_db

    tools=["fastp","Kraken", "Ariba"]
    chosen_tools=[]

    databases=["ncbi","megares","plasmidfinder","card", "srst2_argannot","argannot","resfinder","virulencefinder","vfdb_core","vfdb_full"]
    chosen_db='['

    # START OF COMMAND GENERATION:
    sev=while_yn("Do you wish to run several genomes? [y]/n:\t", "yes")
    if sev.lower().startswith("y"):
        directory=input("Please input name of directory containing genomes: ") or exit("No name given. Exiting.") 
        command=command_generator(directory,"None","directory")
    else:
        name=input("please enter the name (ID) of your strain: ").replace("_2.fastq.gz","").replace("_1.fastq.gz","") or exit("No name given. Exiting.") 
        command=command_generator(name+"_1.fastq.gz",name+"_2.fastq.gz","genome")


    a=while_yn(f"Your command is: \n {command}\n\nDo you wish to run the pipeline imemdiately? y/[n]","no")
    if a.lower is "y":
        os.system(command)
    else:
        print("Command generation finished. Bye!")
    
if __name__ == '__main__':
    main()  