# Align PDB sequences against Harvard plasmid library sequences, and output pretty HTML file for viewing.
#
# Daniel L. Parton <partond@mskcc.org> - 3 Jan 2014
#

import sys, os, openpyxl, re
from lxml import etree
from lxml.builder import E
import TargetExplorer as clab
import Bio.Seq
import Bio.Alphabet

try:
    override_target = sys.argv[ sys.argv.index('-target') + 1 ]
    targets = [override_target]
except ValueError:
    targets = 'All'

# ===========
# Parse DB
# ===========

HOME_DIR = os.environ['HOME']
TargetExplorer_rootdir = os.path.join(HOME_DIR, 'dev/TargetExplorerDB')
DB_path = os.path.join(TargetExplorer_rootdir, 'database', 'database-stage.xml')
DB_root = etree.parse(DB_path).getroot()
nentries = len(DB_root)

# ===========
# Target selection
# ===========

if targets == 'All':
    targets = [ domain_node.get('target_id') for domain_node in DB_root.findall('entry/UniProt/domains/domain[@target_id]') ]

# ===========
# Get Harvard plasmid library sequences
# ===========

plasmid_library_dir = 'Harvard-DFHCC-HIP_human_kinase_collection-pJP1520'
Harvard_plasmid_library_filepath = os.path.join(plasmid_library_dir, 'Mehle_Kinase_VS1_pJP1520_new_plates.xlsx')

wb = openpyxl.load_workbook(Harvard_plasmid_library_filepath)
sheet_ranges = wb.get_sheet_by_name(name = 'Kinase_VS1_pJP1520_new_plates')
nrows = sheet_ranges.get_highest_row()
plasmid_aa_seqs = {}
for row in range(2,nrows):
    NCBI_GeneID = sheet_ranges.cell('H%d' % row).value    # type(int)
    Symbol = sheet_ranges.cell('I%d' % row).value
    DNA_seq = sheet_ranges.cell('L%d' % row).value 
    # Translate to DNA sequence (don't include stop codon)
    aa_seq = Bio.Seq.Seq(DNA_seq, Bio.Alphabet.generic_dna).translate(to_stop=True)
    plasmid_aa_seqs[NCBI_GeneID] = aa_seq

plasmid_NCBI_GeneIDs = plasmid_aa_seqs.keys()

# ===========
# Function definitions
# ===========

def generate_html_from_alignment(title, alignment, alignment_IDs, aa_css_class_list=None):
    html_body = E.body()
    html_body.append( E.h2(title) )
    html_table = E.table()
    html_body.append( html_table )

    max_alignment_ID_len = max( [ len(ID) for ID in alignment_IDs ] )
    for i in range(len(alignment)):
        #table = E.table( E.tr( E.td( E.div(alignment_IDs[i],CLASS='ali')), E.td( E.div(id='blank') ), E.td( E.div(id=alignment_IDs[i],CLASS='ali'),nowrap='' ) ) )
        row = E.tr( E.td( E.div(alignment_IDs[i],CLASS='ali')), E.td( E.div(id=alignment_IDs[i],CLASS='ali'), nowrap='') )
        seq_div = row.find('td/div[@id="%s"]' % alignment_IDs[i])
        seq_div.set('style','background-color:#dddddd;letter-spacing:-5px')
        if i == 0:
            prettyseq = clab.core.seq2pretty_html(alignment[i], aa_css_class_list=aa_css_class_list)
        else:
            prettyseq = clab.core.seq2pretty_html(alignment[i])
        for span in prettyseq:
            seq_div.append(span)

        html_table.append(row)

    return html_body

def write_output(ofilename, html_tree):
    ofile = open(ofilename, 'w')
    ofile.write( etree.tostring(html_tree, pretty_print=True) )
    ofile.close()

# ===========
# Iterate through targets
# ===========

for target in targets:
    print 'Working on target:', target
    # ===========
    # Generate html header
    # ===========
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
    # TODO fix
    css_path = os.path.join('..', '..', '..', 'pylib', 'choderalab', 'seqlib.cs')
    css_link.set('href',css_path)
    css_link.set('rel','stylesheet')

    # get DB entry
    DB_domain = DB_root.find('entry/UniProt/domains/domain[@target_id="%s"]' % target)
    DB_entry = DB_domain.getparent().getparent().getparent()

    # get IDs
    target_UniProt_entry_name = DB_entry.find('UniProt').get('entry_name')
    target_NCBI_GeneID = DB_entry.find('NCBI_Gene/entry').get('ID')
    if target_NCBI_GeneID == None:
        print 'Gene ID not found for target %s' % target_NCBI_GeneID
        continue
    target_NCBI_GeneID = int(target_NCBI_GeneID)

    if target_NCBI_GeneID not in plasmid_NCBI_GeneIDs:
        print 'Gene ID %s (%s) not found in Harvard plasmid library.' % (target_NCBI_GeneID, target_UniProt_entry_name)
        continue

    # get UniProt canonical isoform sequence
    UniProt_canonseq = clab.core.sequnwrap( DB_entry.findtext('UniProt/isoforms/canonical_isoform/sequence') )

    # get PDB sequences and store in dict e.g. { '3GKZ_B' : 'MGYL...' }
    PDB_seqs = { seq_node.getparent().getparent().getparent().get('ID') + '_' + seq_node.getparent().getparent().get('ID') : clab.core.sequnwrap(seq_node.text) for seq_node in DB_entry.findall('PDB/structure/chain/experimental_sequence/sequence') }

    # ===========
    # run MSA (using ClustalO)
    # ===========
    alignment_IDs = ['UniProt', 'plasmid_seq']
    alignment_seqs = [UniProt_canonseq, plasmid_aa_seqs[target_NCBI_GeneID]]
    for PDB_ID in PDB_seqs.keys():
        alignment_IDs.append(PDB_ID)
        alignment_seqs.append(PDB_seqs[PDB_ID])
    alignment = clab.align.run_clustalo(alignment_IDs, alignment_seqs)
    UniProt_canon_seq_aligned = alignment[0]

    # ===========
    # compare plasmid and PDB seqs with aligned UniProt canonical seq - put non-matching aas in lower case
    # ===========
    for seq_iter in range(len(alignment))[1:]:
        seq = list(alignment[seq_iter]) # seq needs to be mutable
        for aa_iter in range(len(seq)):
            if seq[aa_iter] != UniProt_canon_seq_aligned[aa_iter]:
                seq[aa_iter] = seq[aa_iter].lower()
        alignment[seq_iter] = ''.join(seq)

    # ===========
    # find the target domain within the aligned target domain seq
    # ===========
    target_domain_seq = clab.core.sequnwrap( DB_domain.findtext('sequence') )
    # to do this, we construct a regex which accounts for the possible presence of '-' chars within the target domain sequence.
    target_domain_seq_regex = ''.join([ aa + '-*' for aa in target_domain_seq ])[:-2] # ignore last '-*'
    aligned_UniProt_seq_target_domain_match = re.search(target_domain_seq_regex, UniProt_canon_seq_aligned)
    len_aligned_UniProt_seq_target_domain_match = aligned_UniProt_seq_target_domain_match.end() - aligned_UniProt_seq_target_domain_match.start()

    # ===========
    # calculate on the number of aas outside the target domain sequence
    # ===========
    ordered_alignment = alignment[0:2] # the UniProt seq and plasmid seq stay at the beginning
    # now iterate through the PDB construct seqs
    PDB_construct_seqs_aligned = alignment[2:]
    num_aas_outside_target_domain = [0] * len(PDB_construct_seqs_aligned)
    for i in range(len(PDB_construct_seqs_aligned)):
        # calculate number of aas present outside the target domain
        for aa_iter in range(len(alignment[i])):
            if PDB_construct_seqs_aligned[i][aa_iter] != '-' and ( i < aligned_UniProt_seq_target_domain_match.start() or i > aligned_UniProt_seq_target_domain_match.end() ):
                num_aas_outside_target_domain[i] += 1

    # ===========
    # score the alignment quality for the target domain sequence
    # ===========
    aln_scores = [0] * len(PDB_construct_seqs_aligned)
    for i in range(len(PDB_construct_seqs_aligned)):
        # extract target domain sequences
        UniProt_canon_seq_target_domain_seq = UniProt_canon_seq_aligned[ aligned_UniProt_seq_target_domain_match.start() : aligned_UniProt_seq_target_domain_match.end() ]
        PDB_seq_target_domain_seq = PDB_construct_seqs_aligned[i][ aligned_UniProt_seq_target_domain_match.start() : aligned_UniProt_seq_target_domain_match.end() ]
        # compare using PAM matrix
        aln_scores[i] = 0 - clab.align.score_aln(UniProt_canon_seq_target_domain_seq, PDB_seq_target_domain_seq) # subtract from 0 to reverse ordering

    # ===========
    # rank the PDB constructs based firstly on the number of aas outside the target domain sequence, and secondly on the alignment score 
    # ===========
    dict_for_sorting = { seq : (num_aas_outside_target_domain[i], aln_scores[i]) for i, seq in enumerate(PDB_construct_seqs_aligned) }
    PDB_construct_seqs_aligned= sorted( PDB_construct_seqs_aligned, key=dict_for_sorting.get )
    ordered_alignment += PDB_construct_seqs_aligned

    # ===========
    # generate custom set of css classes which will be used to highlight target domain of UniProt sequence in red ('c4')
    # ===========
    aa_css_class_list_UniProt_seq = ['bl'] * len(UniProt_canon_seq_aligned)
    aa_css_class_list_UniProt_seq[ aligned_UniProt_seq_target_domain_match.start() : aligned_UniProt_seq_target_domain_match.end() ] = ['c4'] * len_aligned_UniProt_seq_target_domain_match

    # ===========
    # generate html for the alignment
    # ===========
    alignment_html = generate_html_from_alignment(target, ordered_alignment, alignment_IDs, aa_css_class_list=aa_css_class_list_UniProt_seq)

    # add the alignment html to the main html tree
    for child in alignment_html.getchildren():
        output_html_body.append(child)

    # ===========
    # write to html file
    # ===========
    ofilename = os.path.join('PDB_construct_alignments', target + '.html')
    write_output(ofilename, output_html_tree)

