import os, shutil, StringIO
import pandas as pd
import lxml.html
from lxml import etree
from lxml.html import builder as E
import jinja2

# ========
# Main html page
# ========

title = '96 kinases'
subtitle = 'Selected from SGC Oxford and HIP pJP1520 plasmid libraries'
html_output_cols = ['plasmid_lib', 'cloneID', 'targetID', 'nconflicts_target_domain_region', 'nextraneous_plasmid_residues', 'nPDBs', 'expr_tag', 'auth_score', 'DB_target_rank']

parser = etree.HTMLParser(remove_blank_text=True)

df = pd.DataFrame.from_csv('../PDB_construct_selection/96-kinases-sgc_and_hip.csv')

pd_html = df.to_html(columns=html_output_cols)
html = etree.parse(StringIO.StringIO(pd_html), parser).getroot()
table_children = html.find('body/table').getchildren()

# ========
# Strip whitespace added by pandas (seems to happen due to a bug in pandas)
# Also need to tell tablesorter to sort certain columns numerically (this is necessary if the column also contains non-numerical data, such as 'NaN' from pandas)
# ========

thead = table_children[0]
tbody = table_children[-1]
for tr in thead:
    for element in tr.getchildren():
        try:
            element.text = element.text.strip()
        except AttributeError:
            continue
        if element.text in ['nconflicts_target_domain_region', 'nextraneous_plasmid_residues', 'nPDBs', 'auth_score', 'DB_target_rank']:
            element.set("class", "{sorter: 'digit'}")
for tr in tbody:
    for element in tr.getchildren():
        element.text = element.text.strip()

# ========
# Add href links to alignments for HIP plasmids
# ========

for tr in tbody:
    if 'HIP pJP1520' in tr.findtext('td'):
        cloneID_td = tr.find('td[2]')
        targetID = tr.findtext('td[3]').strip()
        cloneID = cloneID_td.text
        cloneID_td.text = ''
        aelem = etree.SubElement(cloneID_td, 'a')
        aelem.set('href', 'plasmid-construct-alignments/HIP-pJP1520/' + targetID + '.html')
        aelem.text = cloneID
    if 'SGC Oxford' in tr.findtext('td'):
        cloneID_td = tr.find('td[2]')
        cloneID = cloneID_td.text
        cloneID_td.text = ''
        aelem = etree.SubElement(cloneID_td, 'a')
        aelem.set('href', 'plasmid-construct-alignments/SGC-Oxford/' + cloneID.strip() + '.html')
        aelem.text = cloneID

# ========
# Render html from template and output
# ========

table_content = '\n'.join([etree.tostring(element) for element in table_children])

env = jinja2.Environment(loader=jinja2.PackageLoader('app', 'templates'))
template = env.get_template('expression-construct-table.html')

with open('96-kinases-sgc_and_hip.html', 'w') as html_file:
    html_file.write( template.render(title=title, subtitle=subtitle, maintable=table_content) )

# ========
# Copy html alignment files from main branch
# ========

# copy CSS
shutil.copy('stylesheets/seqlib.css', 'plasmid-construct-alignments/HIP-pJP1520/seqlib.css')
shutil.copy('stylesheets/seqlib.css', 'plasmid-construct-alignments/SGC-Oxford/seqlib.css')

HIP_plasmid_targetIDs = list( df[df['plasmid_lib'] == 'HIP pJP1520']['targetID'] )
for targetID in HIP_plasmid_targetIDs:
    src_html_filepath = os.path.join('..', 'PDB_construct_selection', 'alignments', targetID + '.html') # XXX directory location will prob change soon
    dest_html_filepath = os.path.join('plasmid-construct-alignments', 'HIP-pJP1520', targetID + '.html')
    shutil.copy(src_html_filepath, dest_html_filepath)

SGC_plasmid_cloneIDs = list( df[df['plasmid_lib'] == 'SGC Oxford']['cloneID'] )
for cloneID in SGC_plasmid_cloneIDs:
    src_html_filepath = os.path.join('..', 'plasmids', 'SGC', 'Oxford_SGC_Clones', 'alignments', cloneID + '.html')
    dest_html_filepath = os.path.join('plasmid-construct-alignments', 'SGC-Oxford', cloneID + '.html')
    shutil.copy(src_html_filepath, dest_html_filepath)

