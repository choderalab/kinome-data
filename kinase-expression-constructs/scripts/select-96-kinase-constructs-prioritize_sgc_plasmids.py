import os, copy
from lxml import etree
import pandas as pd
import numpy as np

ofilename = '96-kinases-sgc_and_hip'

# ========
# Read in data
# ========
sgc_plasmids = pd.DataFrame.from_csv('../plasmids/SGC/Oxford_SGC_Clones/aln.csv')

# remove SGC plasmids with 100 or more conflicts in the target domain region
selected_sgc_plasmids = sgc_plasmids[ sgc_plasmids['nconflicts_target_domain_region'] < 100 ]
selected_sgc_plasmids.reset_index(inplace=True) # add numerical index and move 'cloneID' to a column
print 'Number of SGC plasmids with < 100 conflicts in the target domain region:', len(selected_sgc_plasmids)

hip_plasmids = pd.DataFrame.from_csv('../plasmids/DFHCC-PlasmID/HIP-human_kinase_collection-pJP1520/aln.csv')
hip_plasmids.reset_index(inplace=True) # add numerical index and move 'cloneID' to a column
print 'Number of HIP plasmids:', len(hip_plasmids)

parser = etree.XMLParser(remove_blank_text=True)
pdbconstructs_xml = etree.parse('../PDB-constructs/PDB_constructs-data.xml', parser).getroot()
print 'Number of targets:', len(pdbconstructs_xml.findall('target'))
print 'Number of PDB constructs:', len(pdbconstructs_xml.findall('target/PDB_construct'))

# ========
# Iterate through targets
# ========
targets_results = {
'targetID':[],
'DB_target_rank':[],
'plasmid_source':[],
'plasmid_ID':[],
'plasmid_nconflicts_in_domain':[],
'plasmid_nextraneous_residues':[],
'nPDBs':[],
'top_pdb_ID':[],
'top_pdb_expr_tag':[],
'top_pdb_auth_score':[],
'top_pdb_nextraneous_residues':[],
'top_pdb_taxon':[],
}

hip_results = copy.deepcopy(targets_results)

# firstly look for targets with SGC plasmids
for target in pdbconstructs_xml:
    targetID = target.get('targetID')
    matching_sgc_plasmids = selected_sgc_plasmids[ selected_sgc_plasmids['matching_targetID'] == targetID ]
    # if found, rank by nextraneous residues and nconflicts, and take the top-ranked plasmid
    if len(matching_sgc_plasmids) > 0:
        matching_sgc_plasmids.sort(('nextraneous_plasmid_residues', 'nconflicts_target_domain_region'), inplace=True)
        chosen_sgc_plasmid = matching_sgc_plasmids.head(1)
        targets_results['targetID'].append(targetID)
        targets_results['plasmid_source'].append('SGC Oxford')
        targets_results['plasmid_ID'].append(chosen_sgc_plasmid['cloneID'].values[0])
        targets_results['plasmid_nconflicts_in_domain'].append(chosen_sgc_plasmid['nconflicts_target_domain_region'].values[0])
        targets_results['plasmid_nextraneous_residues'].append(chosen_sgc_plasmid['nextraneous_plasmid_residues'].values[0])

        targets_results['DB_target_rank'].append(target.get('DB_target_rank'))
        targets_results['nPDBs'].append(target.get('nPDBs'))

        # get pdbconstruct data
        pdbconstructs = target.findall('PDB_construct')
        if len(pdbconstructs) > 0:
            targets_results['top_pdb_ID'].append(pdbconstructs[0].get('PDBconstructID'))
            targets_results['top_pdb_expr_tag'].append(pdbconstructs[0].get('expr_tag_string'))
            targets_results['top_pdb_auth_score'].append(pdbconstructs[0].get('auth_score'))
            targets_results['top_pdb_nextraneous_residues'].append(pdbconstructs[0].get('nextraneous_residues'))
            targets_results['top_pdb_taxon'].append(pdbconstructs[0].get('taxname'))
        else:
            targets_results['top_pdb_ID'].append('-')
            targets_results['top_pdb_expr_tag'].append('-')
            targets_results['top_pdb_auth_score'].append('-')
            targets_results['top_pdb_nextraneous_residues'].append('-')
            targets_results['top_pdb_taxon'].append('-')

targets_results = pd.DataFrame(targets_results)

# for targets without SGC plasmids, look for HIP plasmids
for target in pdbconstructs_xml:
    targetID = target.get('targetID')
    # skip if this target has a matching SGC plasmid
    if targetID in list(targets_results['targetID']):
        continue

    matching_hip_plasmids = hip_plasmids[ hip_plasmids['matching_targetID'] == targetID ]
    if len(matching_hip_plasmids) > 0:
        matching_hip_plasmids.sort(('nextraneous_plasmid_residues', 'nconflicts_target_domain_region'), inplace=True)
        chosen_hip_plasmid = matching_hip_plasmids.head(1)
        hip_results['targetID'].append(targetID)
        hip_results['plasmid_source'].append('HIP pJP1520')
        hip_results['plasmid_ID'].append(chosen_hip_plasmid['cloneID'].values[0])
        hip_results['plasmid_nconflicts_in_domain'].append(chosen_hip_plasmid['nconflicts_target_domain_region'].values[0])
        hip_results['plasmid_nextraneous_residues'].append(chosen_hip_plasmid['nextraneous_plasmid_residues'].values[0])

        hip_results['DB_target_rank'].append(target.get('DB_target_rank'))
        hip_results['nPDBs'].append(target.get('nPDBs'))

        # get pdbconstruct data
        pdbconstructs = target.findall('PDB_construct')
        if len(pdbconstructs) > 0:
            hip_results['top_pdb_ID'].append(pdbconstructs[0].get('PDBconstructID'))
            hip_results['top_pdb_expr_tag'].append(pdbconstructs[0].get('expr_tag_string'))
            hip_results['top_pdb_auth_score'].append(pdbconstructs[0].get('auth_score'))
            hip_results['top_pdb_nextraneous_residues'].append(pdbconstructs[0].get('nextraneous_residues'))
            hip_results['top_pdb_taxon'].append(pdbconstructs[0].get('taxname'))
        else:
            hip_results['top_pdb_ID'].append('-')
            hip_results['top_pdb_expr_tag'].append('-')
            hip_results['top_pdb_auth_score'].append('-')
            hip_results['top_pdb_nextraneous_residues'].append('-')
            hip_results['top_pdb_taxon'].append('-')

hip_results = pd.DataFrame(hip_results)
# filter targets to keep only those with nPDBs > 0
hip_results = hip_results[ hip_results['nPDBs'].astype(int) > 0 ]
print '%d/%d targets have one or matching HIP pJP1520 plasmids, and > 0 matching PDB constructs' % (len(hip_results), len(pdbconstructs_xml))

# filter targets to keep only those with top_pdb_auth_score > 0
hip_results = hip_results[ hip_results['top_pdb_auth_score'].astype(int) > 0 ]
print '%d/%d targets have one or matching HIP pJP1520 plasmids, > 0 matching PDB constructs, and a positive authenticity score' % (len(hip_results), len(pdbconstructs_xml))

# convert 'DB_target_rank' column to int dtype (from str), ready for sorting
hip_results['DB_target_rank'] = hip_results['DB_target_rank'].astype(int)
# sort targets by 'DB_target_rank'
hip_results.sort('DB_target_rank', inplace=True)

# Take the required number of targets with HIP plasmids, and add to the main list of selected targets
ntargets_w_hip_plasmid_required = 96 - len(targets_results)
selected_hip_results = hip_results.head(ntargets_w_hip_plasmid_required)
selected_hip_results.set_index(np.arange(len(selected_hip_results)) + len(targets_results), inplace=True)
targets_results = pd.concat([targets_results, selected_hip_results])


# Write csv file
#targets_results.set_index('targetID', inplace=True)
targets_results.to_csv(ofilename + '.csv')

# Write text file
output_columns=['targetID', 'DB_target_rank', 'plasmid_source', 'plasmid_ID', 'plasmid_nconflicts_in_domain', 'plasmid_nextraneous_residues', 'nPDBs', 'top_pdb_ID', 'top_pdb_expr_tag', 'top_pdb_auth_score', 'top_pdb_nextraneous_residues', 'top_pdb_taxon']
with open(ofilename + '.txt', 'w') as otxtfile:
    otxtfile.write(targets_results.to_string(columns=output_columns))

