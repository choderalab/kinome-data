import pandas as pd
import numpy as np

# ========
# Read in data
# ========
sgc_plasmids = pd.DataFrame.from_csv('../plasmids/SGC/Oxford_SGC_Clones/aln.csv')

selected_sgc_plasmids = sgc_plasmids[ sgc_plasmids['nconflicts_target_domain_region'] < 10 ]
print 'Number of SGC plasmids with < 10 conflicts in the target domain region:', len(selected_sgc_plasmids)
selected_sgc_plasmids.sort('DB_target_rank', inplace=True)
selected_sgc_plasmids.reset_index(inplace=True) # add numerical index and move 'cloneID' to a column

pdb_constructs = pd.DataFrame.from_csv('PDB_constructs-data.csv')


# ========
# Filter out PDB constructs with no matching UniProtAC, and those with an authenticity score <= 0
# ========
pdb_constructs_with_matching_UniProtAC = pdb_constructs['UniProtAC'].notnull()
pdb_constructs = pdb_constructs[ pdb_constructs_with_matching_UniProtAC ]

pdb_constructs_with_positive_auth_score = pdb_constructs['top_cnstrct_auth_score'] > 0
pdb_constructs = pdb_constructs[ pdb_constructs_with_positive_auth_score ]

pdb_constructs_sorted = pdb_constructs.sort(columns='DB_target_rank')
nselected_pdb_constructs = 96 - len(selected_sgc_plasmids)
print '\nNumber of PDB constructs with a matching UniProt entry:', len(pdb_constructs_sorted)

intersecting_UniProt_entry_names = [ x for x in pdb_constructs_sorted['UniProt_entry_name'].values if x in list(selected_sgc_plasmids['UniProt_entry_name']) ]
print '\nThe following %d UniProt entries are present in both SGC Oxford and HIP pJP1520 plasmid libraries. Only the SGC Oxford plasmids will be selected for expression testing.' % len(intersecting_UniProt_entry_names)
print intersecting_UniProt_entry_names

nonintersecting_UniProtACs = [ True if x not in list(selected_sgc_plasmids['UniProtAC']) else False for x in pdb_constructs_sorted['UniProtAC'].values ]
print '\nNumber of PDB constructs with matching HIP pJP1520 plasmids, which do not intersect with the SGC Oxford kinases:', sum(nonintersecting_UniProtACs)


selected_pdb_constructs = pdb_constructs_sorted[ nonintersecting_UniProtACs ]
selected_pdb_constructs = selected_pdb_constructs.head(nselected_pdb_constructs)

selected_pdb_constructs.reset_index(inplace=True)
selected_pdb_constructs.set_index(np.arange(len(selected_sgc_plasmids), 96), inplace=True)

# ========
# Rename columns to be more readable, and so that sgc and pdb column names match
# ========
pdb_rename_dict = {
'top_cloneID': 'cloneID',
'top_plasmid_nconflicts': 'nconflicts_target_domain_region',
'nmatching_PDB_structures':'nPDBs',
'top_cnstrct_expr_tag':'expr_tag',
'top_cnstrct_auth_score':'auth_score',
}
selected_pdb_constructs.rename(columns=pdb_rename_dict, inplace=True)

sgc_rename_dict = {
'matching_targetID': 'targetID',
}
selected_sgc_plasmids.rename(columns=sgc_rename_dict, inplace=True)

# ========
# Add column defining plasmid library source
# ========
selected_sgc_plasmids['plasmid_lib'] = pd.Series(['SGC Oxford'] * len(selected_sgc_plasmids), index=selected_sgc_plasmids.index)
selected_pdb_constructs['plasmid_lib'] = pd.Series(['HIP pJP1520'] * len(selected_pdb_constructs), index=selected_pdb_constructs.index)

# ========
# Output data
# ========

merged = pd.concat([selected_sgc_plasmids, selected_pdb_constructs])

merged.to_csv('96-kinases-sgc_and_hip.csv')

