import os
import pandas as pd
import glob
import csv
import sys
sys.path.append('/home/alma/Documents/kandidat/2d_sum_sev_gen')
# print(os.getcwd())

path='/home/alma/Documents/kandidat/2d_sum_sev_gen' # hardcoded path to all assembly directories. change as needed


extension = 'csv'
all_filenames = [i for i in glob.glob('ass*/out.sum.csv')]
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ], axis=1)
combined_csv.to_csv( f'{path}combined_csv.csv', index=False, encoding='utf-8-sig')

# with open(path/'combined_csv.csv', 'rb') as f:
#     mycsv = csv.reader(f)
#     for row in mycsv:
#         for col in row:

print(combined_csv)