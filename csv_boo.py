import os
import pandas as pd
import glob
import csv
import sys

path='/home/alma/Documents/kandidat/2d_sum_sev_gen' # hardcoded path to all assembly directories. change as needed

extension = 'csv'
common_suffix='addsomethinggood'

all_filenames = [i for i in glob.glob(f'*{common_suffix}.csv')]                         # all csvfiles are 
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ], axis=0)              # axis =0 makes horisontal
combined_csv.to_csv( f'{path}Finallittlebaby.csv', index=False, encoding='utf-8-sig')   # creates final file 