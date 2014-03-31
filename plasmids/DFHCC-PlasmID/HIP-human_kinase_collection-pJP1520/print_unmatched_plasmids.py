import pandas as pd
plasmid_df = pd.DataFrame.from_csv('plasmid_insert_data.csv')

unmatched_plasmids = plasmid_df[ plasmid_df['UniProtAC'] == 'None' ]
print unmatched_plasmids.to_string()

