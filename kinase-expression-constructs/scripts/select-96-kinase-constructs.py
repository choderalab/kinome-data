import sys, os, copy
from lxml import etree
import pandas as pd
import numpy as np

ofilename = '96-kinases-sgc_and_hip'

output_columns=['targetID', 'DB_target_rank', 'plasmid_source', 'plasmid_ID', 'plasmid_nconflicts', 'plasmid_nextraneous_residues', 'nPDBs', 'top_pdb_ID', 'top_pdb_expr_tag', 'top_pdb_auth_score', 'top_pdb_nextraneous_residues', 'family', 'top_pdb_taxname', 'selected_construct_source', 'selected_construct_nextraneous_residues']

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

all_plasmids = pd.concat((selected_sgc_plasmids, hip_plasmids))
all_plasmids.reset_index(inplace=True)

plasmid_source = ['SGC Oxford'] * len(selected_sgc_plasmids) + ['HIP pJP1520'] * len(hip_plasmids)
all_plasmids['plasmid_source'] = pd.Series(plasmid_source)

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
'plasmid_nconflicts':[],
'plasmid_nextraneous_residues':[],
'nPDBs':[],
'top_pdb_ID':[],
'top_pdb_expr_tag':[],
'top_pdb_auth_score':[],
'top_pdb_nextraneous_residues':[],
'family':[],
'top_pdb_taxname':[],
'selected_construct_source':[],
'selected_construct_nextraneous_residues':[],
}



for target in pdbconstructs_xml:
    targetID = target.get('targetID')

    # get plasmid data
    matching_plasmids = all_plasmids[ all_plasmids['matching_targetID'] == targetID ]
    matching_sgc_plasmids = matching_plasmids[ matching_plasmids['plasmid_source'] == 'SGC Oxford' ]
    matching_hip_plasmids = matching_plasmids[ matching_plasmids['plasmid_source'] == 'HIP pJP1520' ]

    matching_plasmids.sort(('nextraneous_plasmid_residues', 'nconflicts_target_domain_region'), inplace=True)
    matching_sgc_plasmids.sort(('nextraneous_plasmid_residues', 'nconflicts_target_domain_region'), inplace=True)
    matching_hip_plasmids.sort(('nextraneous_plasmid_residues', 'nconflicts_target_domain_region'), inplace=True)

    # get pdbconstruct data - this is already sorted
    pdbconstructs = target.findall('PDB_construct')
    top_pdb_id = pdbconstructs[0].get('PDBconstructID') if len(pdbconstructs) > 0 else None
    top_pdb_expr_tag = pdbconstructs[0].get('expr_tag_string') if len(pdbconstructs) > 0 else None
    top_pdb_auth_score = int(pdbconstructs[0].get('auth_score')) if len(pdbconstructs) > 0 else None
    top_pdb_nextraneous_residues = int(pdbconstructs[0].get('nextraneous_residues')) if len(pdbconstructs) > 0 else None
    top_pdb_taxname = pdbconstructs[0].get('taxname') if len(pdbconstructs) > 0 else None

    top_sgc_plasmid_nextraneous_residues = matching_sgc_plasmids['nextraneous_plasmid_residues'].values[0] if len(matching_sgc_plasmids) > 0 else None
    top_sgc_plasmid_nconflicts = matching_sgc_plasmids['nconflicts_target_domain_region'].values[0] if len(matching_sgc_plasmids) > 0 else None
    top_hip_plasmid_nextraneous_residues = matching_hip_plasmids['nextraneous_plasmid_residues'].values[0] if len(matching_hip_plasmids) > 0 else None
    top_hip_plasmid_nconflicts = matching_hip_plasmids['nconflicts_target_domain_region'].values[0] if len(matching_hip_plasmids) > 0 else None

    # if targetID == 'RN5A_HUMAN_D0':
    #     import ipdb; ipdb.set_trace()

    if len(matching_plasmids) == 0:
        continue   # skip if no plasmids at all

    elif len(matching_hip_plasmids) > 0 and (len(pdbconstructs) == 0 or top_pdb_auth_score <= 0):
        continue   # skip if HIP plasmids, but no authentic PDB constructs

    elif len(matching_sgc_plasmids) > 0 and len(matching_hip_plasmids) == 0 and (len(pdbconstructs) == 0 or top_pdb_auth_score <= 0):
        # if SGC plasmids, no HIP plasmids, no authentic PDB constructs - select SGC
        top_plasmid = matching_sgc_plasmids.head(1)
        selected_construct_source = 'SGC'

    elif len(matching_sgc_plasmids) > 0 and len(matching_hip_plasmids) == 0 and len(pdbconstructs) > 0 and top_pdb_auth_score > 0:
        # if SGC plasmids, no HIP plasmids, authentic PDB constructs - select shortest of SGC and PDB
        top_plasmid = matching_sgc_plasmids.head(1)
        if top_pdb_nextraneous_residues < top_sgc_plasmid_nextraneous_residues:
            selected_construct_source = 'PDB'
        else:
            selected_construct_source = 'SGC'

    elif len(matching_sgc_plasmids) == 0 and len(matching_hip_plasmids) > 0 and len(pdbconstructs) > 0 and top_pdb_auth_score > 0:
        # if no SGC plasmids, HIP plasmids, authentic PDB constructs - select PDB (HIP plasmid cannot be trusted to express well)
        top_plasmid = matching_hip_plasmids.head(1)
        selected_construct_source = 'PDB'

    elif len(matching_sgc_plasmids) > 0 and len(matching_hip_plasmids) > 0 and (len(pdbconstructs) == 0 or top_pdb_auth_score <= 0):
        # if SGC and HIP plasmids, no authentic PDB constructs - select SGC plasmid (HIP plasmid cannot be trusted to express well)
        top_plasmid = matching_sgc_plasmids.head(1)
        selected_construct_source = 'SGC'

    elif len(pdbconstructs) > 0 and len(matching_sgc_plasmids) > 0 and len(matching_hip_plasmids) > 0 and top_pdb_auth_score > 0:
        # if SGC and HIP plasmids and authentic PDB constructs are available - further disambiguation required
        if top_pdb_nextraneous_residues < top_sgc_plasmid_nextraneous_residues and top_pdb_nextraneous_residues < top_hip_plasmid_nextraneous_residues:
            # if PDB is shortest - select PDB
            selected_construct_source = 'PDB'
            top_plasmid = matching_plasmids.head(1)
        elif top_hip_plasmid_nextraneous_residues < top_sgc_plasmid_nextraneous_residues and top_hip_plasmid_nextraneous_residues < top_pdb_nextraneous_residues:
            # if HIP is shortest - select PDB or SGC, whichever is shortest (HIP plasmid cannot be trusted to express well)
            top_plasmid = matching_sgc_plasmids.head(1)
            if top_pdb_nextraneous_residues < top_sgc_plasmid_nextraneous_residues:
                selected_construct_source = 'PDB'
            else:
                selected_construct_source = 'SGC'
        else:
            # if SGC is shortest - select SGC
            top_plasmid = matching_sgc_plasmids.head(1)
            selected_construct_source = 'SGC'

    else:
        print 'WARNING for target %s: shouldn\'t get here.' % targetID
        import ipdb; ipdb.set_trace()




    targets_results['targetID'].append(targetID)
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

    if selected_construct_source == 'SGC':
        targets_results['selected_construct_nextraneous_residues'].append(top_sgc_plasmid_nextraneous_residues)
    elif selected_construct_source == 'HIP':
        targets_results['selected_construct_nextraneous_residues'].append(top_hip_plasmid_nextraneous_residues)
    elif selected_construct_source == 'PDB':
        targets_results['selected_construct_nextraneous_residues'].append(top_pdb_nextraneous_residues)



    if targets_results['selected_construct_nextraneous_residues'][-1] > targets_results['plasmid_nextraneous_residues'][-1]:
        if targetID in ['CDK2_HUMAN_D0', 'KAPCA_HUMAN_D0', 'MK01_HUMAN_D0', 'MK03_HUMAN_D0', 'MK13_HUMAN_D0', 'MK14_HUMAN_D0', 'NEK7_HUMAN_D0']:
            additional_message = 'NOTE, this target is ok (checked manually).'
        else:
            additional_message = ''
        print 'WARNING: selected construct is shorter than plasmid for target %s - %s' % (targetID, additional_message)


# # firstly look for targets with SGC plasmids
# for target in pdbconstructs_xml:
#     targetID = target.get('targetID')
#     matching_sgc_plasmids = selected_sgc_plasmids[ selected_sgc_plasmids['matching_targetID'] == targetID ]
#     # if found, rank by nextraneous residues and nconflicts, and take the top-ranked plasmid
#     if len(matching_sgc_plasmids) > 0:
#         matching_sgc_plasmids.sort(('nextraneous_plasmid_residues', 'nconflicts_target_domain_region'), inplace=True)
#         chosen_sgc_plasmid = matching_sgc_plasmids.head(1)
#         targets_results['targetID'].append(targetID)
#         targets_results['plasmid_source'].append('SGC Oxford')
#         targets_results['plasmid_ID'].append(chosen_sgc_plasmid['cloneID'].values[0])
#         targets_results['plasmid_nconflicts'].append(chosen_sgc_plasmid['nconflicts_target_domain_region'].values[0])
#         targets_results['plasmid_nextraneous_residues'].append(chosen_sgc_plasmid['nextraneous_plasmid_residues'].values[0])
#         targets_results['family'].append(chosen_sgc_plasmid['UniProt_family'].values[0])
#
#         targets_results['DB_target_rank'].append(target.get('DB_target_rank'))
#         targets_results['nPDBs'].append(target.get('nPDBs'))
#
#         # get pdbconstruct data
#         pdbconstructs = target.findall('PDB_construct')
#         if len(pdbconstructs) > 0:
#             targets_results['top_pdb_ID'].append(pdbconstructs[0].get('PDBconstructID'))
#             targets_results['top_pdb_expr_tag'].append(pdbconstructs[0].get('expr_tag_string'))
#             targets_results['top_pdb_auth_score'].append(pdbconstructs[0].get('auth_score'))
#             targets_results['top_pdb_nextraneous_residues'].append(pdbconstructs[0].get('nextraneous_residues'))
#             targets_results['top_pdb_taxname'].append(pdbconstructs[0].get('taxname'))
#         else:
#             targets_results['top_pdb_ID'].append(np.nan)
#             targets_results['top_pdb_expr_tag'].append(np.nan)
#             targets_results['top_pdb_auth_score'].append(np.nan)
#             targets_results['top_pdb_nextraneous_residues'].append(np.nan)
#             targets_results['top_pdb_taxname'].append(np.nan)
#
#         # compute lowest nextraneous residues
#         if len(pdbconstructs) > 0:
#             if int(pdbconstructs[0].get('nextraneous_residues')) < chosen_sgc_plasmid['nextraneous_plasmid_residues'].values[0]:
#                 targets_results['lowest_nextraneous_residues'].append(pdbconstructs[0].get('nextraneous_residues'))
#             else:
#                 targets_results['lowest_nextraneous_residues'].append(chosen_sgc_plasmid['nextraneous_plasmid_residues'].values[0])
#         else:
#             targets_results['lowest_nextraneous_residues'].append(chosen_sgc_plasmid['nextraneous_plasmid_residues'].values[0])
#
#
# # for targets without SGC plasmids, look for HIP plasmids
# for target in pdbconstructs_xml:
#     targetID = target.get('targetID')
#     # skip if this target has a matching SGC plasmid
#     if targetID in list(targets_results['targetID']):
#         continue
#
#     matching_hip_plasmids = hip_plasmids[ hip_plasmids['matching_targetID'] == targetID ]
#     if len(matching_hip_plasmids) > 0:
#         matching_hip_plasmids.sort(('nextraneous_plasmid_residues', 'nconflicts_target_domain_region'), inplace=True)
#         chosen_hip_plasmid = matching_hip_plasmids.head(1)
#         targets_results['targetID'].append(targetID)
#         targets_results['plasmid_source'].append('HIP pJP1520')
#         targets_results['plasmid_ID'].append(chosen_hip_plasmid['cloneID'].values[0])
#         targets_results['plasmid_nconflicts'].append(chosen_hip_plasmid['nconflicts_target_domain_region'].values[0])
#         targets_results['plasmid_nextraneous_residues'].append(chosen_hip_plasmid['nextraneous_plasmid_residues'].values[0])
#         targets_results['family'].append(chosen_hip_plasmid['UniProt_family'].values[0])
#
#         targets_results['DB_target_rank'].append(target.get('DB_target_rank'))
#         targets_results['nPDBs'].append(target.get('nPDBs'))
#
#         # get pdbconstruct data
#         pdbconstructs = target.findall('PDB_construct')
#         if len(pdbconstructs) > 0:
#             targets_results['top_pdb_ID'].append(pdbconstructs[0].get('PDBconstructID'))
#             targets_results['top_pdb_expr_tag'].append(pdbconstructs[0].get('expr_tag_string'))
#             targets_results['top_pdb_auth_score'].append(pdbconstructs[0].get('auth_score'))
#             targets_results['top_pdb_nextraneous_residues'].append(pdbconstructs[0].get('nextraneous_residues'))
#             targets_results['top_pdb_taxname'].append(pdbconstructs[0].get('taxname'))
#         else:
#             targets_results['top_pdb_ID'].append(np.nan)
#             targets_results['top_pdb_expr_tag'].append(np.nan)
#             targets_results['top_pdb_auth_score'].append(np.nan)
#             targets_results['top_pdb_nextraneous_residues'].append(np.nan)
#             targets_results['top_pdb_taxname'].append(np.nan)
#
#         # compute lowest nextraneous residues
#         if len(pdbconstructs) > 0:
#             if int(pdbconstructs[0].get('nextraneous_residues')) < chosen_hip_plasmid['nextraneous_plasmid_residues'].values[0]:
#                 targets_results['lowest_nextraneous_residues'].append(pdbconstructs[0].get('nextraneous_residues'))
#             else:
#                 targets_results['lowest_nextraneous_residues'].append(chosen_hip_plasmid['nextraneous_plasmid_residues'].values[0])
#         else:
#             targets_results['lowest_nextraneous_residues'].append(chosen_hip_plasmid['nextraneous_plasmid_residues'].values[0])




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
selected_targets.sort(('selected_construct_source', 'selected_construct_nextraneous_residues'), inplace=True)

# sort targets by 'DB_target_rank'
# selected_targets.sort('DB_target_rank', inplace=True)

selected_targets.reset_index(inplace=True)


# Write csv file
selected_targets.to_csv(ofilename + '.csv')

unfiltered_targets = targets_results[ targets_results['plasmid_nconflicts'] < 25 ]
unfiltered_targets.to_csv(ofilename + '-unfiltered.csv')

# Write text file
with open(ofilename + '.txt', 'w') as otxtfile:
    otxtfile.write(targets_results.to_string(columns=output_columns))

