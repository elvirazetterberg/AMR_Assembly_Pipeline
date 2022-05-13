import os
import pandas as pd
import glob
import csv
import sys
"""
concatenates csv files
"""
inputfil="SRR18825428_1.fastq.gz\nSRR18825428_2.fastq.gz\nSRR19139723_1.fastq.gz\nSRR19139723_2.fastq.gz"
every_other_elements = [inputfil.split('\n')[index].strip('_2.fastq.gz') for index in range(1, len(inputfil.split('\n')), 2)]

finalname="finallittlebaby"
    
infopath= os.getcwd() # correct? where are we standing?
all_filenames = [i for i in glob.glob(f'assembly*/*assembly/info.csv')]  

combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ], axis=0) 
# combined_csv.to_csv( f'{infopath}/{finalname}.csv', index=False, encoding='utf-8-sig')

# df = pd.read_csv(f"{finalname}.csv")
print(combined_csv)
combined_csv["Genome Name"] = every_other_elements
print(combined_csv)
combined_csv.to_csv( f'{infopath}/{finalname}.csv', index=False, encoding='utf-8-sig')
