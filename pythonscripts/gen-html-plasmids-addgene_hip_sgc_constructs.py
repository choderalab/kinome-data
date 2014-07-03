import os, shutil, StringIO
import pandas as pd
from lxml import etree
import jinja2

# ========
# Main html page
# ========

title = 'kinase constructs'
subtitle = 'Selected from addgene kinase collection and HIP pJP1520 plasmid libraries'
html_output_cols = ['targetID', 'DB_target_rank', 'plasmid_ID', 'plasmid_source', 'plasmid_nconflicts', 'plasmid_nextraneous_residues', 'nPDBs', 'top_pdb_ID', 'top_pdb_expr_tag', 'top_pdb_auth_score', 'top_pdb_nextraneous_residues', 'family', 'top_pdb_taxname', 'selected_construct_source', 'selected_construct_nextraneous_residues']

parser = etree.HTMLParser(remove_blank_text=True)

df = pd.DataFrame.from_csv('../kinase-expression-constructs/addgene_hip_sgc_constructs/selected-kinases.csv')

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
            aelem.set('href', 'plasmid-construct-alignments/addgene_hip_sgc_constructs/' + targetID + '.html')
            aelem.text = targetID
            break
        elif 'addgene' in td.text:
            targetID_td = tr.find('td[1]')
            targetID = targetID_td.text
            targetID_td.text = ''
            aelem = etree.SubElement(targetID_td, 'a')
            aelem.set('href', 'plasmid-construct-alignments/addgene_hip_sgc_constructs/' + targetID + '.html')
            aelem.text = targetID
            break
        elif 'SGC' in td.text:
            targetID_td = tr.find('td[1]')
            targetID = targetID_td.text
            targetID_td.text = ''
            aelem = etree.SubElement(targetID_td, 'a')
            aelem.set('href', 'plasmid-construct-alignments/addgene_hip_sgc_constructs/' + targetID + '.html')
            aelem.text = targetID
            break

# ========
# Render html from template and output
# ========

table_content = '\n'.join([etree.tostring(element) for element in table_children])

env = jinja2.Environment(loader=jinja2.PackageLoader('app', 'templates'))
template = env.get_template('expr-cnstrct-custom-selection-infotext-addgene_hip_sgc_constructs.html')

with open('kinase_constructs-addgene_hip_sgc_constructs.html', 'w') as html_file:
    html_file.write( template.render(title=title, subtitle=subtitle, maintable=table_content, display_filters=False) )

# ========
# Copy html alignment files from main branch
# ========

# copy CSS
shutil.copy('stylesheets/seqlib.css', 'plasmid-construct-alignments/addgene_hip_sgc_constructs/seqlib.css')

targetIDs = list(df['targetID'])
for targetID in targetIDs:
    src_html_filepath = os.path.join('..', 'kinase-expression-constructs', 'addgene_hip_sgc_constructs', 'alignments', targetID + '.html')
    dest_html_filepath = os.path.join('plasmid-construct-alignments', 'addgene_hip_sgc_constructs', targetID + '.html')
    shutil.copy(src_html_filepath, dest_html_filepath)

