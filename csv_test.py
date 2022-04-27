import pandas as pd

fil = 'test.csv'

df = pd.read_csv(fil, header=None)
ny = df.iloc[:,1]
ny2 = df.iloc[:,3]

df2 = pd.concat([ny,ny2], axis=1)
print(df2)