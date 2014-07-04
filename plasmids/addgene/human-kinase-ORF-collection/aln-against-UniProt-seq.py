import sys, os, re
from operator import itemgetter
import openpyxl
from openpyxl import Workbook
import Bio, Bio.Seq, Bio.Alphabet, Bio.SubsMat.MatrixInfo, Bio.pairwise2
from lxml import etree
from lxml.builder import E
import TargetExplorer as clab
import pandas as pd
from numpy import nan

# TODO sort out plasmids with no or v few aas - may need to search for the correct ORF

DB_path = os.path.join('/Users/partond', 'kinomeDB', 'database', 'database.xml')
DB_root = etree.parse(DB_path).getroot()

plasmid_df = pd.read_csv('plasmid_insert_data.csv')


DB_targets = DB_root.findall('entry/UniProt/domains/domain[@targetID]')

# To be used to construct a pandas DataFrame
output_data = {
'matching_targetID':[],
'cloneID':[],
'DB_target_rank':[],
'UniProt_family':[],
'UniProt_seq_aln':[],
'clone_seq_aln':[],
'nconflicts_target_domain_region':[],
'nextraneous_plasmid_residues':[],
'len_seq_aln':[],
'pctidentity_target_domain_region':[],
'construct_aa_seq':[],
'construct_dna_seq':[],
'construct_dna_orf_seq':[],
}

# ===========
# Iterate through plasmids
# ===========

for p in plasmid_df.index:
    cloneID = plasmid_df['cloneID'][p]
    UniProtAC = plasmid_df['UniProtAC'][p]
    UniProt_entry_name = plasmid_df['UniProt_entry_name'][p]
    insert_dna_seq = plasmid_df['insert_dna_seq'][p]
   #  print UniProt_entry_name, cloneID
    if type(insert_dna_seq) != str or insert_dna_seq == 'None':
       #  print 'No DNA seq found - skipping.'
        continue

    if UniProtAC == 'None':
        continue

    DB_UniProt_node = DB_root.find('entry/UniProt[@AC="%s"]' % UniProtAC)
    DB_entry = DB_UniProt_node.getparent()
    UniProt_family = DB_UniProt_node.get('family')
    UniProt_seq = ''.join(DB_UniProt_node.findtext('isoforms/canonical_isoform/sequence').strip().split('\n'))

    domains = DB_UniProt_node.findall('domains/domain[@targetID]')

    # =====
    # Find ORF and translate to aa sequence
    # =====

    dna_seq_frames = [insert_dna_seq[i:] for i in range(3)]
    orf_dna_seqs = [''] * 3
    for frame in range(3):
        in_frame = False
        for i in range(0, len(dna_seq_frames[frame]), 3):
            codon = dna_seq_frames[frame][i:i+3]
            if codon == 'ATG':
                in_frame = True
            elif codon in ['TAG', 'TAA', 'TGA']:
                in_frame = False
            if in_frame:
                orf_dna_seqs[frame] += codon

    orf_dna_seq_lens = [len(x) for x in orf_dna_seqs]

    orf_aa_seqs = [ Bio.Seq.translate(orf_dna_seqs[i], to_stop=True) for i in range(3) ]

    longest_orf_index = max(enumerate(orf_dna_seq_lens), key=lambda x: x[1])[0]

    # set the selected sequence as the longest aa seq derived from the three DNA seq frames
    insert_dna_orf_seq = orf_dna_seqs[longest_orf_index]
    insert_aa_seq = orf_aa_seqs[longest_orf_index]

    # ignore if the translated aa sequence is particularly short
    if len(insert_aa_seq) < 30:
        # print 'aa seq < 30 residues: %s' % insert_aa_seq
        continue

    # skip if DNA sequence contains non-standard letters
    non_standard_dna_letter = False
    for nucleotide in set(insert_dna_orf_seq):
        if nucleotide not in ['A', 'C', 'T', 'G']:
            non_standard_dna_letter = True
            break
    if non_standard_dna_letter:
        continue

    # orf_regex = '^ATG[A-Z]*(TAG|TAA|TGA)'
    # for i in range(0, len(insert_dna_orf_seq), 3):
    #     codon = insert_dna_orf_seq[i:i+3]
    #     if codon in ['TAG', 'TAA', 'TGA']:
    #         insert_dna_orf_seq = insert_dna_orf_seq[:i]
    #         break
        # re_search = re.search(orf_regex, insert_dna_seq)
        # if re_search != None:
        #     insert_dna_orf_seq = insert_dna_seq[re_search.start() : re_search.end()]
        #     break

    # print len(insert_dna_seq), len(insert_dna_orf_seq)
    if abs(len(insert_dna_seq) - len(insert_dna_orf_seq)) not in [0,3]:
        print UniProt_entry_name, cloneID
        print longest_orf_index
        print '+   ' + '\n+   '.join(orf_aa_seqs)
        print len(insert_dna_seq), len(insert_dna_orf_seq)
        print insert_dna_seq
        print ' ' * len(insert_dna_orf_seq) + '^'
        print insert_dna_orf_seq

    print ''

    # =====
    # Run alignment
    # =====
    matrix = Bio.SubsMat.MatrixInfo.gonnet
    gap_open = -10
    gap_extend = -0.5
    try:
        aln = Bio.pairwise2.align.globalds(UniProt_seq, insert_aa_seq, matrix, gap_open, gap_extend)
    except:
        print insert_aa_seq
        import ipdb; ipdb.set_trace()
    aln = [aln[0][0], aln[0][1]]


    # alnIDs = [targetID, cloneID]
    # pre_aln_seqs = [domain_seq, insert_aa_seq]
    # if nan in pre_aln_seqs:
    #     continue
    #
    # try:
    #     aligned_seqs = clab.align.run_clustalo(alnIDs, pre_aln_seqs)
    # except Exception as e:
    #     print alnIDs
    #     print pre_aln_seqs
    #     print [type(x) for x in pre_aln_seqs]
    #     raise e
    # alignment = [ [alnIDs[i], seq] for i, seq in enumerate(aligned_seqs) ]




    # Calculate the number of plasmid residues outside the target domain region (excluding expression tags)
    # Also use this to determine which target domain the plasmid is most likely to represent
    nextraneous_plasmid_residues = []
    domain_seqs = [''.join(domain.findtext('sequence').strip().split('\n')) for domain in domains]
    for d in range(len(domains)):
        nextraneous_plasmid_residues.append(0)
        domain_seq_regex = re.compile( ''.join( [ aa + '-*' for aa in domain_seqs[d] ] ) )
        UniProt_domain_aln_coords = re.search(domain_seq_regex, aln[0])
        for a in range(len(aln[0])):
            # print a, aln[1][a], UniProt_domain_aln_coords.start(), UniProt_domain_aln_coords.end()
            if (a < UniProt_domain_aln_coords.start() and aln[1][a] != '-') or (a >= UniProt_domain_aln_coords.end() and aln[1][a] != '-'):
                nextraneous_plasmid_residues[-1] += 1
    domainID, nextraneous_plasmid_residues = min(enumerate(nextraneous_plasmid_residues), key=itemgetter(1))
    targetID = UniProt_entry_name + '_D' + str(domainID)
    domain_seq = domain_seqs[domainID]

    # Make mismatching residues in plasmid sequence lower case
    plasmid_aln_list = list(aln[1])
    for a in range(len(aln[0])):
        if aln[1][a] != aln[0][a]:
            plasmid_aln_list[a] = aln[1][a].lower()
    aln[1] = ''.join(plasmid_aln_list)



    UniProt_domain_aln_coords = re.search(domain_seq_regex, aln[0])

    # Count conflicting residues within the target domain region
    nconflicts = []
    pctidentity = []
    nconflicts.append(0)
    for a in range(UniProt_domain_aln_coords.start(), UniProt_domain_aln_coords.end()):
        if aln[0][a] == '-' or aln[1][a] == '-':
            nconflicts[-1] += 1
        elif aln[0][a].upper() != aln[1][a].upper():
            nconflicts[-1] += 1

    # Calculate percent identity
    pctidentity.append( (len(domain_seq) - nconflicts[-1]) * 100. / float(len(domain_seq)) )

    nconflicts = nconflicts[domainID]
    pctidentity = pctidentity[domainID]

    DB_target_score_node = DB_root.find('entry/target_score/domain[@targetID="%s"]' % targetID)
    DB_target_rank = DB_target_score_node.get('target_rank')





    # # ===========
    # # convert conflicts in plasmid seq to lower case and find the target domain region in the alignment coords
    # # ===========
    # seq = list(alignment[1][1]) # sequence needs to be mutable so convert to list
    # for aa_iter in range(len(seq)):
    #     if seq[aa_iter] != alignment[0][1][aa_iter]:
    #         seq[aa_iter] = seq[aa_iter].lower()
    # alignment[1][1] = ''.join(seq)
    #
    # #regex_string = ''.join([ aa + '-+' for aa in domain_seq ])
    # #regex_obj = re.search(domain_seq, alignment[0][1])
    # #if regex_obj == None:
    # #    print domain_seq
    # #    print alignment[0][1]
    # #    print alignment[1][1]
    # #target_domain_start_aln_coords = regex_obj.start()
    # #target_domain_end_aln_coords = regex_obj.end()
    #
    # regex_start = re.search('[A-Z]', alignment[0][1])
    # regex_end = re.search('[A-Z]', alignment[0][1][::-1])
    # target_domain_start_aln_coords = regex_start.start()
    # target_domain_end_aln_coords = len(alignment[0][1]) - regex_end.start() - 1  # 0-based index of the last residue
    #
    # # ===========
    # # calculate number of conflicts and percentage
    # # ===========
    # aln_seqs_target_domain_region = [aln_seq[1][target_domain_start_aln_coords : target_domain_end_aln_coords + 1] for aln_seq in alignment]
    # lcase_regex = re.compile('[a-z]')
    # nconflicts = len( [ aa_iter for aa_iter in range(len(aln_seqs_target_domain_region[0])) if aln_seqs_target_domain_region[0][aa_iter] == '-' or aln_seqs_target_domain_region[1][aa_iter] == '-' or re.match(lcase_regex, aln_seqs_target_domain_region[1][aa_iter]) ] )
    # pctidentity = (len(domain_seq) - nconflicts) * 100. / float(len(domain_seq))

    # ===========
    # append to data dict, to be used later to construct pandas DataFrame
    # ===========

    output_data['matching_targetID'].append(targetID)
    output_data['cloneID'].append(cloneID)
    output_data['DB_target_rank'].append(DB_target_rank)
    output_data['UniProt_family'].append(UniProt_family)
    output_data['UniProt_seq_aln'].append(aln[0])
    output_data['clone_seq_aln'].append(aln[1])
    output_data['nconflicts_target_domain_region'].append(nconflicts)
    output_data['nextraneous_plasmid_residues'].append(nextraneous_plasmid_residues)
    output_data['len_seq_aln'].append(str(len(domain_seq)))
    output_data['pctidentity_target_domain_region'].append('%.2f' % pctidentity)
    output_data['construct_aa_seq'].append(insert_aa_seq)
    output_data['construct_dna_seq'].append(insert_dna_seq)
    output_data['construct_dna_orf_seq'].append(insert_dna_orf_seq)

# construct pandas DataFrame and write to csv
output_data = pd.DataFrame(output_data)
#df.to_csv('aln-against-UniProt-seq.csv', index=False)

output_data = output_data.sort('nconflicts_target_domain_region', ascending=True)
output_data.index = range(len(output_data)) # redo indices, as these will have been sorted in the previous step

# write data txt file
with open('aln.txt', 'w') as otxtfile:
    columnwidths = {}
    for key in output_data.keys():
        columnwidths[key] = max( [ len(str(x)) for x in output_data[key] ] )

    IDcolumnwidth = max( [columnwidths['matching_targetID'], columnwidths['cloneID']] )

    for i in range(len(output_data)):
        otxtfile.write('%-*s' % (IDcolumnwidth, output_data['matching_targetID'][i]))
        otxtfile.write('  ')
        otxtfile.write('%*s' % (columnwidths['pctidentity_target_domain_region'], '') )
        otxtfile.write('  ')
        otxtfile.write('%*s' % (columnwidths['nconflicts_target_domain_region'], '') )
        otxtfile.write(' ')
        otxtfile.write('%-*s' % (columnwidths['len_seq_aln'], '') )
        otxtfile.write('  ')
        otxtfile.write('%s' % output_data['UniProt_seq_aln'][i])
        otxtfile.write('\n')

        otxtfile.write('%-*s' % (IDcolumnwidth, output_data['cloneID'][i]))
        otxtfile.write('  ')
        otxtfile.write('%*s' % (columnwidths['pctidentity_target_domain_region'], output_data['pctidentity_target_domain_region'][i]) )
        otxtfile.write('  ')
        otxtfile.write('%*s' % (columnwidths['nconflicts_target_domain_region'], str(output_data['nconflicts_target_domain_region'][i])) )
        otxtfile.write('/')
        otxtfile.write('%-*s' % (columnwidths['len_seq_aln'], output_data['len_seq_aln'][i]) )
        otxtfile.write('  ')
        otxtfile.write('%s' % output_data['clone_seq_aln'][i])
        otxtfile.write('\n\n')

#with open('aln.html', 'w') as ohtmlfile:
#    ohtmlfile.write( etree.tostring(output_html_tree, pretty_print=True) )

output_data.set_index('cloneID', inplace=True)
output_data.to_csv('aln.csv')

