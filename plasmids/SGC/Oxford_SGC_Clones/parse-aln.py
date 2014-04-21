import pandas as pd
df = pd.DataFrame.from_csv('aln.csv')

nconflicts = df['nconflicts_target_domain_region']

for i in range(0,20):
    print 'Number of plasmids with nconflicts == %d:' % i, sum(nconflicts == i)

print 'Number of plasmids with nconflicts < 100:', sum(nconflicts < 100)
