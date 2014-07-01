import sys, os, copy
from lxml import etree
import pandas as pd
import numpy as np

ofilename = 'selected-kinases-addgene_and_hip'

output_columns=['targetID', 'DB_target_rank', 'plasmid_source', 'plasmid_ID', 'plasmid_nconflicts', 'plasmid_nextraneous_residues', 'nPDBs', 'top_pdb_ID', 'top_pdb_expr_tag', 'top_pdb_auth_score', 'top_pdb_nextraneous_residues', 'family', 'top_pdb_taxname', 'selected_construct_source', 'selected_construct_nextraneous_residues']

# ========
# Read in data
# ========
addgene_plasmids = pd.DataFrame.from_csv('../../plasmids/addgene/human-kinase-ORF-collection/aln.csv')

# remove addgene plasmids with 100 or more conflicts in the target domain region
selected_addgene_plasmids = addgene_plasmids[ addgene_plasmids['nconflicts_target_domain_region'] < 100 ]
selected_addgene_plasmids.reset_index(inplace=True) # add numerical index and move 'cloneID' to a column
print 'Number of addgene plasmids with < 100 conflicts in the target domain region:', len(selected_addgene_plasmids)

hip_plasmids = pd.DataFrame.from_csv('../../plasmids/DFHCC-PlasmID/HIP-human_kinase_collection-pJP1520/aln.csv')
hip_plasmids.reset_index(inplace=True) # add numerical index and move 'cloneID' to a column
print 'Number of HIP plasmids:', len(hip_plasmids)

all_plasmids = pd.concat((selected_addgene_plasmids, hip_plasmids))
all_plasmids.reset_index(inplace=True)

plasmid_source = ['addgene'] * len(selected_addgene_plasmids) + ['HIP pJP1520'] * len(hip_plasmids)
all_plasmids['plasmid_source'] = pd.Series(plasmid_source)

parser = etree.XMLParser(remove_blank_text=True)
pdbconstructs_xml = etree.parse('../../PDB-constructs/PDB_constructs-data.xml', parser).getroot()
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
'target_UniProt_seq':[],
'DB_target_rank':[],
'plasmid_source':[],
'plasmid_ID':[],
'plasmid_nconflicts':[],
'plasmid_nextraneous_residues':[],
'plasmid_aa_seq':[],
# 'plasmid_dna_orf_seq':[],
'plasmid_dna_seq':[],
'nPDBs':[],
'top_pdb_ID':[],
'top_pdb_expr_tag':[],
'top_pdb_auth_score':[],
'top_pdb_nextraneous_residues':[],
'top_pdb_aa_seq':[],
'family':[],
'top_pdb_taxname':[],
'selected_construct_source':[],
'selected_construct_nextraneous_residues':[],
# 'selected_construct_aa_seq':[],
# 'selected_construct_aa_start':[],
# 'selected_construct_aa_end':[],
# 'selected_construct_dna_seq':[],
# 'selected_construct_dna_start':[],
# 'selected_construct_dna_end':[],
}



for target in pdbconstructs_xml:
    targetID = target.get('targetID')
    target_uniprot_ac = target.get('UniProtAC')
    target_uniprot_seq = target.find('seq').text

    # get plasmid data
    matching_plasmids = all_plasmids[ all_plasmids['matching_targetID'] == targetID ]
    matching_addgene_plasmids = matching_plasmids[ matching_plasmids['plasmid_source'] == 'addgene' ]
    matching_hip_plasmids = matching_plasmids[ matching_plasmids['plasmid_source'] == 'HIP pJP1520' ]

    matching_plasmids.sort(('nconflicts_target_domain_region', 'nextraneous_plasmid_residues'), inplace=True)
    matching_addgene_plasmids.sort(('nconflicts_target_domain_region', 'nextraneous_plasmid_residues'), inplace=True)
    matching_hip_plasmids.sort(('nconflicts_target_domain_region', 'nextraneous_plasmid_residues',), inplace=True)

    # get pdbconstruct data - this is already sorted
    pdbconstructs = target.findall('PDB_construct')
    top_pdb_id = pdbconstructs[0].get('PDBconstructID') if len(pdbconstructs) > 0 else None
    top_pdb_expr_tag = pdbconstructs[0].get('expr_tag_string') if len(pdbconstructs) > 0 else None
    top_pdb_auth_score = int(pdbconstructs[0].get('auth_score')) if len(pdbconstructs) > 0 else None
    top_pdb_nextraneous_residues = int(pdbconstructs[0].get('nextraneous_residues')) if len(pdbconstructs) > 0 else None
    top_pdb_taxname = pdbconstructs[0].get('taxname') if len(pdbconstructs) > 0 else None
    top_pdb_aa_seq = pdbconstructs[0].find('seq').text if len(pdbconstructs) > 0 else None

    top_addgene_plasmid_nextraneous_residues = matching_addgene_plasmids['nextraneous_plasmid_residues'].values[0] if len(matching_addgene_plasmids) > 0 else None
    top_addgene_plasmid_nconflicts = matching_addgene_plasmids['nconflicts_target_domain_region'].values[0] if len(matching_addgene_plasmids) > 0 else None
    top_addgene_plasmid_aa_seq = matching_addgene_plasmids['construct_aa_seq'].values[0] if len(matching_addgene_plasmids) > 0 else None
   #  top_addgene_plasmid_dna_orf_seq = matching_addgene_plasmids['construct_dna_orf_seq'].values[0] if len(matching_addgene_plasmids) > 0 else None
    top_addgene_plasmid_dna_seq = matching_addgene_plasmids['construct_dna_seq'].values[0] if len(matching_addgene_plasmids) > 0 else None
   #  top_addgene_expr_tag = matching_addgene_plasmids['Expression tag'].values[0] if len(matching_addgene_plasmids) > 0 else None
    top_hip_plasmid_nextraneous_residues = matching_hip_plasmids['nextraneous_plasmid_residues'].values[0] if len(matching_hip_plasmids) > 0 else None
    top_hip_plasmid_nconflicts = matching_hip_plasmids['nconflicts_target_domain_region'].values[0] if len(matching_hip_plasmids) > 0 else None
    top_hip_plasmid_aa_seq = matching_hip_plasmids['construct_aa_seq'].values[0] if len(matching_hip_plasmids) > 0 else None
    top_hip_plasmid_dna_orf_seq = matching_hip_plasmids['construct_dna_orf_seq'].values[0] if len(matching_hip_plasmids) > 0 else None
    top_hip_plasmid_dna_seq = matching_hip_plasmids['construct_dna_seq'].values[0] if len(matching_hip_plasmids) > 0 else None

    # if targetID == 'RN5A_HUMAN_D0':
    #     import ipdb; ipdb.set_trace()

    if len(matching_plasmids) == 0:
        continue   # skip if no plasmids at all

    elif top_pdb_auth_score <= 0:
        continue   # skip if no authentic PDB constructs

    else:
        selected_construct_source = 'PDB'
        top_plasmid = matching_plasmids.head(1)

    top_plasmid_nextraneous_residues = top_plasmid['nextraneous_plasmid_residues'].values[0]
    top_plasmid_nconflicts = top_plasmid['nconflicts_target_domain_region'].values[0]
    top_plasmid_aa_seq = top_plasmid['construct_aa_seq'].values[0]
    # top_plasmid_dna_orf_seq = top_plasmid['construct_dna_orf_seq'].values[0]
    top_plasmid_dna_seq = top_plasmid['construct_dna_seq'].values[0]


    targets_results['targetID'].append(targetID)
    targets_results['target_UniProt_seq'].append(target_uniprot_seq)
    targets_results['selected_construct_source'].append(selected_construct_source)
    targets_results['plasmid_source'].append(top_plasmid['plasmid_source'].values[0])
    targets_results['plasmid_ID'].append(top_plasmid['cloneID'].values[0])
    targets_results['plasmid_nconflicts'].append(top_plasmid['nconflicts_target_domain_region'].values[0])
    targets_results['plasmid_nextraneous_residues'].append(top_plasmid['nextraneous_plasmid_residues'].values[0])
    targets_results['family'].append(top_plasmid['UniProt_family'].values[0])

    targets_results['DB_target_rank'].append(target.get('DB_target_rank'))
    targets_results['nPDBs'].append(target.get('nPDBs'))

    targets_results['top_pdb_ID'].append(top_pdb_id)
    targets_results['top_pdb_expr_tag'].append(top_pdb_expr_tag)
    targets_results['top_pdb_auth_score'].append(top_pdb_auth_score)
    targets_results['top_pdb_nextraneous_residues'].append(top_pdb_nextraneous_residues)
    targets_results['top_pdb_taxname'].append(top_pdb_taxname)
    targets_results['top_pdb_aa_seq'].append(top_pdb_aa_seq)

    # if targetID == 'CDKL1_HUMAN_D0':
    #     import ipdb; ipdb.set_trace()

    targets_results['plasmid_aa_seq'].append(top_plasmid_aa_seq)
    # targets_results['plasmid_dna_orf_seq'].append(top_plasmid_dna_orf_seq)
    targets_results['plasmid_dna_seq'].append(top_plasmid_dna_seq)

    targets_results['selected_construct_nextraneous_residues'].append(top_plasmid_nextraneous_residues)



    if targets_results['selected_construct_nextraneous_residues'][-1] > targets_results['plasmid_nextraneous_residues'][-1]:
        if targetID in ['CDK2_HUMAN_D0', 'KAPCA_HUMAN_D0', 'MK01_HUMAN_D0', 'MK03_HUMAN_D0', 'MK13_HUMAN_D0', 'MK14_HUMAN_D0', 'NEK7_HUMAN_D0']:
            additional_message = 'NOTE, this target is ok (checked manually).'
        else:
            additional_message = ''
        print 'WARNING: selected construct is shorter than plasmid for target %s - %s' % (targetID, additional_message)




targets_results = pd.DataFrame(targets_results)
# convert 'DB_target_rank' column to int dtype (from str), ready for sorting
targets_results['DB_target_rank'] = targets_results['DB_target_rank'].astype(int)
# the following columns have to be converted to 'float64', since they can contain np.NaN
targets_results['plasmid_nextraneous_residues'] = targets_results['plasmid_nextraneous_residues'].astype(float)
targets_results['top_pdb_auth_score'] = targets_results['top_pdb_auth_score'].astype(float)
targets_results['top_pdb_nextraneous_residues'] = targets_results['top_pdb_nextraneous_residues'].astype(float)
targets_results['selected_construct_nextraneous_residues'] = targets_results['selected_construct_nextraneous_residues'].astype(float)

print '%d/%d targets have a plasmid with < 40 extraneous residues' % (len(targets_results[ targets_results['plasmid_nextraneous_residues'] < 40 ]), len(pdbconstructs_xml))
print '%d/%d targets have nPDBs > 0 and top_pdb_auth_score > 0' % (len( targets_results[ elem_all_true( [targets_results['nPDBs'] > 0, targets_results['top_pdb_auth_score'] > 0] ) ] ), len(pdbconstructs_xml))

# filter targets
# print targets_results.sort('plasmid_nextraneous_residues').reset_index().to_string(columns=output_columns)
# selected_targets = targets_results[ targets_results['plasmid_nconflicts'] < 40 ]
# selected_targets = selected_targets[ selected_targets['plasmid_nextraneous_residues'] < 40 ]
# selected_targets = selected_targets[ selected_targets['DB_target_rank'] < 450 ]
# selected_targets.sort('DB_target_rank', inplace=True)
# selected_targets = pd.concat((selected_targets, targets_results[ elem_all_true( [targets_results['nPDBs'] > 0, targets_results['top_pdb_auth_score'] > 0, targets_results['plasmid_nextraneous_residues'] >= 40 ] ) ].sort('top_pdb_nextraneous_residues')))

# nextraneous_residues_cutoff = 10000
# plasmid_targets = targets_results[ elem_all_true( [targets_results['plasmid_source'] == 'SGC Oxford', targets_results['plasmid_nextraneous_residues'] < nextraneous_residues_cutoff] ) ]
# plasmid_targets.sort('DB_target_rank', inplace=True)
# pdb_targets = targets_results[ elem_all_true( [targets_results['top_pdb_auth_score'] > 0, targets_results['top_pdb_nextraneous_residues'] < nextraneous_residues_cutoff] ) ]
# pdb_targets.sort('DB_target_rank', inplace=True)
# selected_targets = pd.concat((plasmid_targets, pdb_targets))
# selected_targets = targets_results[ elem_any_true( [ elem_all_true( [targets_results['plasmid_source'] == 'SGC Oxford', targets_results['plasmid_nextraneous_residues'] < nextraneous_residues_cutoff] ), elem_all_true( [targets_results['top_pdb_auth_score'] > 0, targets_results['top_pdb_nextraneous_residues'] < nextraneous_residues_cutoff] ) ] ) ]

# selected_targets = targets_results[ elem_any_true( [ targets_results['plasmid_source'] == 'SGC Oxford', targets_results['top_pdb_auth_score'] > 0 ] ) ]
selected_targets = targets_results[ [True] * len(targets_results) ]
selected_targets.sort(('plasmid_source', 'selected_construct_nextraneous_residues'), inplace=True)

# sort targets by 'DB_target_rank'
# selected_targets.sort('DB_target_rank', inplace=True)

selected_targets.reset_index(inplace=True)


# Write csv file
selected_targets.to_csv(ofilename + '.csv')

# unfiltered_targets = targets_results[ targets_results['plasmid_nconflicts'] < 25 ]
# unfiltered_targets.to_csv(ofilename + '-unfiltered.csv')

# Write text file
with open(ofilename + '.txt', 'w') as otxtfile:
    otxtfile.write(targets_results.to_string(columns=output_columns))

