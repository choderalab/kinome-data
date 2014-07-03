import sys, os, urllib, urllib2
from lxml import etree
from StringIO import StringIO
import Bio.Seq, Bio.Alphabet
import TargetExplorer as clab
import pandas as pd
import yaml

parser = etree.HTMLParser()
maxreadlength = 10000000 # 10 MB

# read in manual exceptions file
manual_exceptions_filepath = 'manual_exceptions.yaml'
manual_exceptions = yaml.load( open(manual_exceptions_filepath, 'r').read() )

# =============
# definitions
# =============

def skip_manual_exceptions(plasmidID):
    plasmid_manual_exception = clab.core.parse_nested_dicts(manual_exceptions, [int(plasmidID), 'behavior'])
    if plasmid_manual_exception != None:
        if plasmid_manual_exception == 'ignore':
            manual_exception_comment = clab.core.parse_nested_dicts(manual_exceptions, [int(plasmidID), 'comment'])
            manual_exception_plasmid_name = clab.core.parse_nested_dicts(manual_exceptions, [int(plasmidID), 'name'])
            if manual_exception_comment != None and manual_exception_plasmid_name != None:
                print 'Skipping plasmid %s for reason: %s' % (manual_exception_plasmid_name, manual_exception_comment.strip())
            return True
    return False


# =============
# parse 6 plate pages and scrape plasmid names e.g. pDONR223-SNF1LK and URLs e.g. http://www.addgene.org/23838/
# =============

nplates = 6
plasmid_data_from_platepage = {
'plasmid_name' : [],
'plasmidID' : [],
}

print 'Scraping data from plasmid plate pages...'

for p in range(1,nplates+1):
    plate_url = 'http://www.addgene.org/human-kinase/Plate' + str(p)
    response = urllib2.urlopen(plate_url)
    page = response.read(maxreadlength)
    html_tree = etree.parse(StringIO(page), parser).getroot()
    plasmid_nodes = html_tree.findall('body/div/div/div/div[@id="body"]/table[@class="table-list std-hdr"]/tr/td/a')

    # iterate through plasmid nodes
    for plasmid_node in plasmid_nodes:
        plasmid_name = plasmid_node.text
        plasmidID = plasmid_node.get('href').replace('/', '')  # plasmid_node.get('href') returns e.g. '/23820/'
        plasmid_data_from_platepage['plasmidID'].append(plasmidID)
        plasmid_data_from_platepage['plasmid_name'].append(plasmid_name)

print 'Information obtained for %s plasmids.' % len(plasmid_data_from_platepage['plasmid_name'])


# =============
# iterate through plasmid pages
# =============

# plasmid output data
# data_fields = ['plasmidID', 'plasmid_name', 'NCBI_GeneID', 'NCBI_Gene_name', 'UniProtAC', 'UniProt_entry_name', 'insert_dna_seq', 'insert_aa_seq']
data_fields = ['plasmidID', 'plasmid_name', 'NCBI_GeneID', 'NCBI_Gene_name', 'UniProtAC', 'UniProt_entry_name', 'insert_dna_seq']
output_data = pd.DataFrame( [['None'] * len(data_fields)] * len(plasmid_data_from_platepage['plasmid_name']), columns=data_fields)

print 'Scraping data from plasmid pages...'

for p in output_data.index:
    plasmidID = plasmid_data_from_platepage['plasmidID'][p]
    plasmid_name = plasmid_data_from_platepage['plasmid_name'][p]
    output_data['plasmidID'][p] = plasmidID
    output_data['plasmid_name'][p] = plasmid_name

    # if skip_manual_exceptions(plasmidID):
    #     continue

    plasmid_seq_url = 'http://www.addgene.org/' + plasmid_data_from_platepage['plasmidID'][p]
    response = urllib2.urlopen(plasmid_seq_url)
    page = response.read(maxreadlength)
    html_tree = etree.parse(StringIO(page), parser).getroot()

    NCBI_Gene_node = html_tree.xpath('body/div/div/div/div/div/table/tr/td[@id="data"]/ul/li/label[text()="Entrez Gene:"]/../p/a')

    if len(NCBI_Gene_node) == 0:
        print 'WARNING: NCBI Gene node not found for plasmid name %s ID %s' % (output_data['plasmid_name'][p], output_data['plasmidID'][p])
        continue
    NCBI_Gene_href = NCBI_Gene_node[0].get('href')
    if 'http://www.ncbi.nlm.nih.gov/gene/' not in NCBI_Gene_href:
        raise Exception, 'Unexpected href found in plasmid NCBI Gene node for plasmid name %s ID %s: %s' % (output_data['plasmid_name'][p], output_data['plasmidID'][p], NCBI_Gene_href)
    NCBI_Gene_name = NCBI_Gene_node[0].text
    NCBI_GeneID = NCBI_Gene_href.replace('http://www.ncbi.nlm.nih.gov/gene/', '')

    output_data['NCBI_Gene_name'][p] = NCBI_Gene_name
    output_data['NCBI_GeneID'][p] = NCBI_GeneID

# =============
# iterate through plasmid sequence pages
# =============

print 'Scraping data from plasmid sequence pages...'

for p in output_data.index:
    plasmidID = output_data['plasmidID'][p]

    # if skip_manual_exceptions(plasmidID):
    #     continue

    plasmid_seq_url = 'http://www.addgene.org/' + plasmid_data_from_platepage['plasmidID'][p] + '/sequences'
    response = urllib2.urlopen(plasmid_seq_url)
    page = response.read(maxreadlength)
    html_tree = etree.parse(StringIO(page), parser).getroot()

    dna_seq_text_node = html_tree.find('body/div/div/div/div[@id="body"]/div[@class="sequence-div clear"]/textarea')
    if dna_seq_text_node == None:
        print 'Skipping as DNA sequence text not found for plasmid name %s ID %s.'  % (output_data['plasmid_name'][p], output_data['plasmidID'][p])
        continue

    dna_seq_text = dna_seq_text_node.text.strip().split('\n')

    if dna_seq_text[0][0:16] != '>Author sequence':
        print 'Skipping due to unexpected text found in DNA sequence box for plasmid name %s ID %s: %s' % (output_data['plasmid_name'][p], output_data['plasmidID'][p], dna_seq_text)
        continue

    dna_seq = ''.join(dna_seq_text[1:]).replace(' ', '')
   #  aa_seq = Bio.Seq.Seq(dna_seq, Bio.Alphabet.generic_dna).translate(to_stop=True)
    output_data['insert_dna_seq'][p] = dna_seq
   #  output_data['insert_aa_seq'][p] = aa_seq



# =============
# match plasmids to DB entries
# =============

print 'Matching plasmids to DB entries...'

DB_path = os.path.join('..', '..', '..', '..', 'database', 'database.xml')
DB_root = etree.parse(DB_path).getroot()

for p in output_data.index:
    plasmidID = output_data['plasmidID'][p]

    # if skip_manual_exceptions(plasmidID):
    #     continue

    # match using NCBI Gene ID
    plasmid_NCBI_GeneID = output_data['NCBI_GeneID'][p]

    DB_entry = DB_root.find('entry/NCBI_Gene/entry[@ID="%s"]/../..' % plasmid_NCBI_GeneID)
    if DB_entry == None:
        print 'Matching DB entry not found for Gene ID %s plasmid name %s' % (plasmid_NCBI_GeneID, output_data['plasmid_name'][p])
        continue

    DB_UniProt_node = DB_entry.find('UniProt')
    DB_domains = DB_UniProt_node.findall('domains/domain[@targetID]')
    UniProtAC = DB_UniProt_node.get('AC')
    UniProt_entry_name = DB_UniProt_node.get('entry_name')
    UniProt_canonseq = clab.core.sequnwrap( DB_UniProt_node.find('isoforms/canonical_isoform/sequence').text )

    output_data['UniProtAC'][p] = UniProtAC
    output_data['UniProt_entry_name'][p] = UniProt_entry_name


# =============
# output
# =============

output_data = pd.DataFrame(output_data)
output_data.to_csv('plasmid_insert_data.csv')

print 'Done.'


