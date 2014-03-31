import pandas as pd

df = pd.read_csv('plasmid_insert_data.csv')

for p in df.index:
    print df['plasmid_name'][p]

