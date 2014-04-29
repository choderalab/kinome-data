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
html_output_cols = ['targetID', 'DB_target_rank', 'plasmid_source', 'plasmid_ID', 'plasmid_nconflicts_in_domain', 'plasmid_nextraneous_residues', 'nPDBs', 'top_pdb_ID', 'top_pdb_expr_tag', 'top_pdb_auth_score', 'top_pdb_nextraneous_residues', 'top_pdb_taxon']

parser = etree.HTMLParser(remove_blank_text=True)

df = pd.DataFrame.from_csv('../kinase-expression-constructs/96-kinases-sgc_and_hip-unfiltered.csv')

pd_html = df.to_html(columns=html_output_cols, index=True, na_rep='--')
html = etree.parse(StringIO.StringIO(pd_html), parser).getroot()

# ========
# Strip whitespace added by pandas (seems to happen due to a bug in pandas)
# Also need to tell tablesorter to sort certain columns numerically (this is necessary if the column also contains non-numerical data, such as 'NaN' from pandas)
# ========

table_children = html.find('body/table').getchildren()
thead = table_children[0]
tbody = table_children[-1]
for tr in thead:
    for element in tr.getchildren():
        try:
            element.text = element.text.strip()
        except AttributeError:
            continue
        # if element.text in ['plasmid_nconflicts_in_domain', 'plasmid_nextraneous_residues', 'nPDBs', 'top_pdb_auth_score', 'top_pdb_nextraneous_residues', 'DB_target_rank']:
        #     element.set("class", "{sorter: 'digit'}")
for tr in tbody:
    for element in tr.getchildren():
        element.text = element.text.strip()

# ========
# Add href links to alignments for HIP plasmids
# ========

for tr in tbody:
    for td in tr:
        if 'HIP pJP1520' in td.text:
            targetID_td = tr.find('td[1]')
            targetID = targetID_td.text
            targetID_td.text = ''
            aelem = etree.SubElement(targetID_td, 'a')
            aelem.set('href', 'plasmid-construct-alignments/' + targetID + '.html')
            aelem.text = targetID
            break
        elif 'SGC Oxford' in td.text:
            targetID_td = tr.find('td[1]')
            targetID = targetID_td.text
            targetID_td.text = ''
            aelem = etree.SubElement(targetID_td, 'a')
            aelem.set('href', 'plasmid-construct-alignments/' + targetID + '.html')
            aelem.text = targetID
            break

# ========
# Render html from template and output
# ========

table_content = '\n'.join([etree.tostring(element) for element in table_children])

env = jinja2.Environment(loader=jinja2.PackageLoader('app', 'templates'))
template = env.get_template('expression-construct-table.html')

with open('96-kinases-sgc_and_hip-unfiltered.html', 'w') as html_file:
    html_file.write( template.render(title=title, subtitle=subtitle, maintable=table_content) )

# ========
# Copy html alignment files from main branch
# ========

# copy CSS
# shutil.copy('stylesheets/seqlib.css', 'plasmid-construct-alignments/HIP-pJP1520/seqlib.css')
# shutil.copy('stylesheets/seqlib.css', 'plasmid-construct-alignments/SGC-Oxford/seqlib.css')
shutil.copy('stylesheets/seqlib.css', 'plasmid-construct-alignments/seqlib.css')

targetIDs = list(df['targetID'])
for targetID in targetIDs:
    src_html_filepath = os.path.join('..', 'kinase-expression-constructs', 'alignments', targetID + '.html') # XXX directory location will prob change soon
    dest_html_filepath = os.path.join('plasmid-construct-alignments', targetID + '.html')
    shutil.copy(src_html_filepath, dest_html_filepath)

# SGC_plasmid_cloneIDs = list( df[df['plasmid_source'] == 'SGC Oxford']['cloneID'] )
# for cloneID in SGC_plasmid_cloneIDs:
#     src_html_filepath = os.path.join('..', 'kinase-expression-constructs', 'alignments', targetID + '.html')
#     dest_html_filepath = os.path.join('plasmid-construct-alignments', 'SGC-Oxford', cloneID + '.html')
#     shutil.copy(src_html_filepath, dest_html_filepath)

