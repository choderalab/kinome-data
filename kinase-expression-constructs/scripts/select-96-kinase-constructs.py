import os, copy
from lxml import etree
import pandas as pd
import numpy as np

ofilename = '96-kinases-sgc_and_hip'

output_columns=['targetID', 'DB_target_rank', 'plasmid_source', 'plasmid_ID', 'plasmid_nconflicts_in_domain', 'plasmid_nextraneous_residues', 'nPDBs', 'top_pdb_ID', 'top_pdb_expr_tag', 'top_pdb_auth_score', 'top_pdb_nextraneous_residues', 'top_pdb_taxon']

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

def elem_all_true(inputarrays):
    assert len(inputarrays[0]) == len(inputarrays[-1])
    inputarrays = [list(inputarray) for inputarray in inputarrays]
    result = [True] * len(inputarrays[0])
    result = np.array(result)
    for i in range(len(inputarrays[0])):
        for inputarray in inputarrays:
            if inputarray[i] == False:
                result[i] = False
    return result

def elem_any_true(inputarrays):
    assert len(inputarrays[0]) == len(inputarrays[-1])
    inputarrays = [list(inputarray) for inputarray in inputarrays]
    result = [False] * len(inputarrays[0])
    result = np.array(result)
    for i in range(len(inputarrays[0])):
        for inputarray in inputarrays:
            if inputarray[i] == True:
                result[i] = True
    return result

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
            targets_results['top_pdb_ID'].append(np.nan)
            targets_results['top_pdb_expr_tag'].append(np.nan)
            targets_results['top_pdb_auth_score'].append(np.nan)
            targets_results['top_pdb_nextraneous_residues'].append(np.nan)
            targets_results['top_pdb_taxon'].append(np.nan)

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
        targets_results['targetID'].append(targetID)
        targets_results['plasmid_source'].append('HIP pJP1520')
        targets_results['plasmid_ID'].append(chosen_hip_plasmid['cloneID'].values[0])
        targets_results['plasmid_nconflicts_in_domain'].append(chosen_hip_plasmid['nconflicts_target_domain_region'].values[0])
        targets_results['plasmid_nextraneous_residues'].append(chosen_hip_plasmid['nextraneous_plasmid_residues'].values[0])

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
            targets_results['top_pdb_ID'].append(np.nan)
            targets_results['top_pdb_expr_tag'].append(np.nan)
            targets_results['top_pdb_auth_score'].append(np.nan)
            targets_results['top_pdb_nextraneous_residues'].append(np.nan)
            targets_results['top_pdb_taxon'].append(np.nan)

targets_results = pd.DataFrame(targets_results)
# convert 'DB_target_rank' column to int dtype (from str), ready for sorting
targets_results['DB_target_rank'] = targets_results['DB_target_rank'].astype(int)
# the following columns have to be converted to 'float64', since they can contain np.NaN
targets_results['plasmid_nextraneous_residues'] = targets_results['plasmid_nextraneous_residues'].astype(float)
targets_results['top_pdb_auth_score'] = targets_results['top_pdb_auth_score'].astype(float)
targets_results['top_pdb_nextraneous_residues'] = targets_results['top_pdb_nextraneous_residues'].astype(float)

print '%d/%d targets have a plasmid with < 40 extraneous residues' % (len(targets_results[ targets_results['plasmid_nextraneous_residues'] < 40 ]), len(pdbconstructs_xml))
print '%d/%d targets have nPDBs > 0 and top_pdb_auth_score > 0' % (len( targets_results[ elem_all_true( [targets_results['nPDBs'] > 0, targets_results['top_pdb_auth_score'] > 0] ) ] ), len(pdbconstructs_xml))

# filter targets
# print targets_results.sort('plasmid_nextraneous_residues').reset_index().to_string(columns=output_columns)
selected_targets = targets_results[ targets_results['plasmid_nconflicts_in_domain'] < 40 ]
selected_targets = selected_targets[ selected_targets['plasmid_nextraneous_residues'] < 40 ]
selected_targets.sort('DB_target_rank', inplace=True)
selected_targets = pd.concat((selected_targets, targets_results[ elem_all_true( [targets_results['nPDBs'] > 0, targets_results['top_pdb_auth_score'] > 0, targets_results['plasmid_nextraneous_residues'] >= 40 ] ) ].sort('top_pdb_nextraneous_residues')))
# selected_targets = selected_targets[ selected_targets['DB_target_rank'] < 400 ]

# sort targets by 'DB_target_rank'
# selected_targets.sort('DB_target_rank', inplace=True)

selected_targets.reset_index(inplace=True)


# Write csv file
selected_targets.to_csv(ofilename + '.csv')

# Write text file
with open(ofilename + '.txt', 'w') as otxtfile:
    otxtfile.write(targets_results.to_string(columns=output_columns))

