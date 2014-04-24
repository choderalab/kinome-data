import pandas as pd

output_columns = ['targetID', 'top_plasmid_nconflicts', 'nmatching_PDB_structures', 'top_PDB_chain_ID', 'top_cnstrct_expr_tag', 'top_cnstrct_auth_score', 'top_cnstrct_expr_sys', 'top_cnstrct_taxname', 'DB_target_rank', 'DB_target_score']

df = pd.read_csv('PDB_constructs-data.csv')

# print df.sort(columns=['DB_target_rank']).to_string(columns=output_columns)

selection = df[df['top_cnstrct_auth_score'] > 0]
print selection.sort(columns=['DB_target_rank']).to_string(columns=output_columns)
