import sys, os, re
from lxml import etree
from lxml.builder import E
import pandas as pd
import TargetExplorer

# ========
# input, params etc.
# ========

# get plasmid data
plasmid_filepaths = ['../../plasmids/addgene/human-kinase-ORF-collection/aln.csv', '../../plasmids/DFHCC-PlasmID/HIP-human_kinase_collection-pJP1520/aln.csv']

plasmid_data = pd.DataFrame.from_csv(plasmid_filepaths[0])
for p in range(1, len(plasmid_filepaths)):
    df = pd.DataFrame.from_csv(plasmid_filepaths[p])
    plasmid_data = pd.concat([plasmid_data, df])

# get PDB construct data
pdbconstruct_data = etree.parse('../../PDB-constructs/PDB_constructs-data.xml').getroot()

alignments_dir = 'alignments'

# ========
# defs
# ========

def gen_html(title, alignment, alignment_IDs, additional_data_fields=[], aa_css_class_list=None, sections={}):
    '''
    additional_data_fields structure: [ { [PDB_ID]_[PDB_CHAIN_ID] : data }, { [PDB_ID]_[PDB_CHAIN_ID] : data }, ... ] with a separate dict for each data type, and where data is of len(alignment). Column headers not required.
    aa_css_class_list can be used to override the CSS classes assigned to each residue. Should be given as a list of lists, with shape: (len(alignment), len(alignment[0])).
    sections structure: {int(position): [section_title, additional_field_header1, additional_field_header2, ...]}
    '''
    output_html_tree = E.html(
        E.head(
            E.link()
        ),
        E.body(
            E.h2(title),
            E.table()
        )
    )

    link = output_html_tree.find('head/link')
    link.set('type', 'text/css')
    link.set('href', 'seqlib.css')
    link.set('rel', 'stylesheet')
    html_body = output_html_tree.find('body')
    html_table = html_body.find('table')

    for i in range(len(alignment)):
        if i in sections:
            #html_table.append( E.tr( E.td( E.div( E.strong(E.u(sections[i])) ), nowrap='' ) ) )
            tr = etree.SubElement(html_table, 'tr')
            for col in range(len(sections[i])):
                tr.append( E.td( E.div( E.strong(E.u(sections[i][col])) ), nowrap='' ) )
        # row
        row = E.tr()

        # alignment ID div
        row.append( E.td( E.div( str(alignment_IDs[i]) ,CLASS='ali') ) )

        # add any additional data fields (also in divs)
        for d in range(len(additional_data_fields)):
            data = additional_data_fields[d][alignment_IDs[i]]
            if data != None:
                row.append( E.td( E.div(data, CLASS='ali'), nowrap='') )
            else:
                row.append( E.td( E.div('', CLASS='ali'), nowrap='') )

        # format sequence with css classes. Returned as a list of span objects
        if aa_css_class_list != None:
            if aa_css_class_list[i] != None:
                prettyseq = TargetExplorer.core.seq2pretty_html(alignment[i], aa_css_class_list=aa_css_class_list[i])
            else:
                prettyseq = TargetExplorer.core.seq2pretty_html(alignment[i])
        else:
            prettyseq = TargetExplorer.core.seq2pretty_html(alignment[i])

        # set up sequence div
        seq_div = E.div(id='sequence',CLASS='ali')
        seq_div.set('style','background-color:#dddddd;letter-spacing:-5px')

        # add sequence to div
        for span in prettyseq:
            seq_div.append(span)

        row.append( E.td( seq_div, nowrap='' ) )

        # add the row to the table
        html_table.append(row)

    return output_html_tree


def process_target(t):
    target = pdbconstruct_data[t]
    targetID = target.get('targetID')
    UniProt_seq = target.findtext('seq')
    target_domain_seq = target.findtext('domain_seq')

    print 'Working on target:', targetID

    plasmids = plasmid_data[ plasmid_data['matching_targetID'] == targetID ]
   #  print plasmids['ID']
   #  print [p for p in plasmids['matching_targetID'].values if p[0] != 'H']

    if len(plasmids) == 0:
        return None

    # sort plasmids
    plasmids.sort(('nconflicts_target_domain_region', 'nextraneous_plasmid_residues'), inplace=True)

    cloneIDs = list(plasmids.index)
    plasmid_seqs = list(plasmids['construct_aa_seq'])

    pdbconstructs = target.findall('PDB_construct')
    pdbconstructIDs = [x.get('PDBconstructID') for x in pdbconstructs]
    pdbconstruct_seqs = [x.findtext('seq') for x in pdbconstructs]

    alignment_IDs = [targetID] + cloneIDs + pdbconstructIDs
    pre_alignment_seqs = [UniProt_seq] + plasmid_seqs + pdbconstruct_seqs

    # run alignment
    aln = TargetExplorer.align.run_clustalo(alignment_IDs, pre_alignment_seqs)
    UniProt_seq_aln = aln[0]

    # find domain start and end points in the aln coords
    target_domain_seq_regex = ''.join([ aa + '-*' for aa in target_domain_seq ])[:-2] # ignore last '-*'
    aligned_UniProt_seq_target_domain_match = re.search(target_domain_seq_regex, UniProt_seq_aln)
    # len_aligned_UniProt_seq_target_domain_match = aligned_UniProt_seq_target_domain_match.end() - aligned_UniProt_seq_target_domain_match.start()
    domain_start = aligned_UniProt_seq_target_domain_match.start()
    domain_end = aligned_UniProt_seq_target_domain_match.end()

    # for plasmid and PDB seqs, put non-matching aas in lower case
    for seq_iter in range(len(aln)):
        if seq_iter == 0:
            continue
        seq = list(aln[seq_iter]) # sequence needs to be mutable so convert to list
        for aa_iter in range(len(seq)):
            if seq[aa_iter] != UniProt_seq_aln[aa_iter]:
                seq[aa_iter] = seq[aa_iter].lower()
        aln[seq_iter] = ''.join(seq)

    # generate custom set of css classes which will be used to highlight target domain of UniProt sequence in red ('c4')
    aa_css_class_list_UniProt_seq = [None] * len(aln)
    aa_css_class_list_UniProt_seq[0] = ['bl'] * len(aln[0])
    aa_css_class_list_UniProt_seq[0][ domain_start : domain_end ] = ['c4'] * (domain_end - domain_start)

    # generate section titles/positions for the alignment html
    sections_dict = {0: ['UniProt seq'], 1: ['Plasmids', 'nconf', 'nextran'], len(plasmids) + 1: ['PDB constructs', 'nconf', 'nextran', 'expr_tag', 'organism']}

    # generate additional data fields to be added to the alignment html
    additional_data = [ {targetID: None} for i in range(4) ]
    for i in range(len(cloneIDs)):
        print plasmids['nconflicts_target_domain_region']
        additional_data[0][cloneIDs[i]] = str(int(plasmids['nconflicts_target_domain_region'].values[i]))
        additional_data[1][cloneIDs[i]] = str(int(plasmids['nextraneous_plasmid_residues'].values[i]))
        additional_data[2][cloneIDs[i]] = None
        additional_data[3][cloneIDs[i]] = None
    for i in range(len(pdbconstructIDs)):
        additional_data[0][pdbconstructIDs[i]] = pdbconstructs[i].get('nconflicts_target_domain_region')
        additional_data[1][pdbconstructIDs[i]] = pdbconstructs[i].get('nextraneous_residues')
        additional_data[2][pdbconstructIDs[i]] = pdbconstructs[i].get('expr_tag_string')
        additional_data[3][pdbconstructIDs[i]] = pdbconstructs[i].get('taxname')

    # generate and write html
    aln_html = gen_html(title=targetID, alignment=aln, alignment_IDs=alignment_IDs, aa_css_class_list=aa_css_class_list_UniProt_seq, additional_data_fields=additional_data, sections=sections_dict)
    ohtml_filename = os.path.join(alignments_dir, targetID + '.html')
    with open(ohtml_filename, 'w') as ohtml_file:
        ohtml_file.write( etree.tostring(aln_html, pretty_print=True) )

    return aln


# ========
# main
# ========

if __name__ == '__main__':
    import multiprocessing
    pool =  multiprocessing.Pool()
    targets_results = pool.map(process_target, range(len(pdbconstruct_data)))

