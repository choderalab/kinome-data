import os, argparse, re
import TargetExplorer
from lxml import etree
from lxml.builder import E
import pandas as pd
import Bio.pairwise2
import Bio.SubsMat.MatrixInfo

css_path = 'seqlib.css'

# ========
# Command-line args
# ========

argparser = argparse.ArgumentParser()
argparser.add_argument('--database_path', type=str, help='Path to a TargetExplorer database XML file', required=True)
args = argparser.parse_args()

# ========
# Read in plasmid data
# ========

df = pd.read_csv('plasmid-data.csv')

# ========
# Read in database
# ========

DB_root = etree.parse(args.database_path)

# ========
# HTML layout
# ========

output_html_tree = E.html(
    E.head(
        E.link()
    ),
    E.body(
    )
)

output_html_body = output_html_tree.find('body')
css_link = output_html_tree.find('head/link')
css_link.set('type','text/css')
css_link.set('href',css_path)
css_link.set('rel','stylesheet')

def gen_html(aln, alnIDs, aa_css_class_list=None):
    '''
    additional_data_fields structure: [ { [PDB_ID]_[PDB_CHAIN_ID] : data }, { [PDB_ID]_[PDB_CHAIN_ID] : data }, ... ] with a separate dict for each data type, and where data is of len(alignment). Column headers not required.
    aa_css_class_list can be used to override the CSS classes assigned to each residue. Should be given as a list of lists, with shape: (len(alignment), len(alignment[0])).
    '''
    html_table = E.table(STYLE='margin-bottom: 1cm')

    for r in range(len(aln)):
        row = E.tr()
        row.append( E.td( E.div(alnIDs[r], CLASS='ali'), nowrap='' ) )
        if aa_css_class_list != None:
            if aa_css_class_list[r] != None:
                prettyseq = TargetExplorer.core.seq2pretty_html(aln[r], aa_css_class_list=aa_css_class_list[r])
            else:
                prettyseq = TargetExplorer.core.seq2pretty_html(aln[r])
        else:
            prettyseq = TargetExplorer.core.seq2pretty_html(aln[r])
        seq_div = E.div(id='sequence', CLASS='ali')
        seq_div.set('style','background-color:#dddddd;letter-spacing:-5px')
        for span in prettyseq:
            seq_div.append(span)
        row.append( E.td( seq_div, nowrap='' ) )
        html_table.append(row)

    return html_table

# ========
# Iterate through plasmids
# ========

ofile = open('aln.txt', 'w')

for i in range(len(df)):
# for i in range(2):
    plasmid_data = df.loc[i]
    cloneID = plasmid_data['cloneID']
    print 'Working on cloneID %s' % cloneID
    UniProtAC = plasmid_data['UniProtAC']
    DB_entry = DB_root.find('entry/UniProt[@AC="%s"]/..' % UniProtAC)
    UniProt_entry_name = DB_entry.find('UniProt').get('entry_name')

    UniProt_seq = ''.join(DB_entry.findtext('UniProt/isoforms/canonical_isoform/sequence').strip().split('\n'))
    plasmid_aa_seq = plasmid_data['aa_seq']

    # Separate expression tag from the plasmid insert sequence

    expr_tag_regex = '(^MG{0,1}HHHHHHSSGVD[A-Z]*GTENLYFQSM)|(^MGSSHHHHHHSSGRENLYFQGHM)|(^MHHHHHHSSGRENLYFQG)'
    expr_tag_match = re.search(expr_tag_regex, plasmid_aa_seq)
    if expr_tag_match != None:
        expr_tag_seq = plasmid_aa_seq[ slice(*expr_tag_match.span()) ]
        plasmid_insert_seq = plasmid_aa_seq[expr_tag_match.end() : ]
    else:
        expr_tag_seq = None
        plasmid_insert_seq = plasmid_aa_seq

    # Conduct alignment
    matrix = Bio.SubsMat.MatrixInfo.gonnet
    gap_open = -10
    gap_extend = -0.5
    aln = Bio.pairwise2.align.globalds(UniProt_seq, plasmid_insert_seq, matrix, gap_open, gap_extend)
    aln = [aln[0][0], aln[0][1]]

    # Add expression tag back into the alignmed plasmid seq
    if expr_tag_seq != None:
        plasmid_seq_aln_aa_start = re.search('[A-Za-z]', aln[1]).start()
        UniProt_aln_list = list(aln[0])
        plasmid_aln_list = list(aln[1])
        plasmid_seq_aln_expr_tag_start = plasmid_seq_aln_aa_start - len(expr_tag_seq)
        # where the expression tag extends beyond the aligned sequence, just add '-' for now
        if plasmid_seq_aln_expr_tag_start < 0:
            for a in range(plasmid_seq_aln_expr_tag_start, 0):
                UniProt_aln_list.insert(0, '-')
                plasmid_aln_list.insert(0, '-')
            plasmid_seq_aln_expr_tag_start = 0
        # now that the alignments are the correct length, add in the expression tag sequence
        for a in range(len(expr_tag_seq)):
            plasmid_aln_list[plasmid_seq_aln_expr_tag_start + a] = expr_tag_seq[a].lower()

        aln[0] = ''.join(UniProt_aln_list)
        aln[1] = ''.join(plasmid_aln_list)

    # Make mismatching residues in plasmid sequence lower case
    plasmid_aln_list = list(aln[1])
    for a in range(len(aln[0])):
        if aln[1][a] != aln[0][a]:
            plasmid_aln_list[a] = aln[1][a].lower()
    aln[1] = ''.join(plasmid_aln_list)

    # find UniProt domains and generate custom CSS assignments to highlight target domains of UniProt sequence in red ('c4')
    domains = DB_entry.findall('UniProt/domains/domain[@targetID]')
    aa_css_class_list = [None] * len(aln)
    aa_css_class_list[0] = ['bl'] * len(aln[0])
    for domain in domains:
        domain_seq = ''.join(domain.findtext('sequence').strip().split('\n'))
        domain_seq_regex = ''.join( [ aa + '-*' for aa in domain_seq ] )
        UniProt_aln_domain = re.search(domain_seq_regex, aln[0])
        aa_css_class_list[0][slice(*UniProt_aln_domain.span())] = ['c4'] * (UniProt_aln_domain.end() - UniProt_aln_domain.start())

    # Write to aligned sequences to text file
    ofile.write(aln[0] + '\n')
    ofile.write(aln[1] + '\n\n')

    # Generate html
    alnIDs = [UniProt_entry_name, cloneID]
    html_table = gen_html(aln, alnIDs, aa_css_class_list=aa_css_class_list)
    output_html_body.append(html_table)

ofile.close()

# write html
with open('aln.html', 'w') as htmlfile:
    htmlfile.write( etree.tostring(output_html_tree, pretty_print=True) )

