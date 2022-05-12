import os
import pandas as pd
import glob
import csv
import sys
"""
Takes 1 dimensional csvfiles with the suffix "commonsuffix" and adds them together horisontally.
"""


path='/home/alma/Documents/kandidat/2d_sum_sev_gen' # hardcoded path to all assembly directories. change as needed

#Ariba 2 col per db, i 1a antal hits, i andra gennamn

extension = 'csv'
common_suffix='addsomethinggood'


all_filenames = [i for i in glob.glob(f'assembly*/*/info.csv')]             # all csvfiles are inside the assembly folder and then their folder?

print(all_filenames)

exit()
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ], axis=0)              # axis =0 makes horisontal 1 makes vert
combined_csv.to_csv( f'{path}Finallittlebaby.csv', index=False, encoding='utf-8-sig')   # creates final file 

