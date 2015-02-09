import pandas as pd

df = pd.read_excel('results/470_Report_Expression_Test_Results.xls')

with open('results.txt', 'w') as tablefile:
    tablefile.write(df.to_string())
df.to_pickle('results.p')
