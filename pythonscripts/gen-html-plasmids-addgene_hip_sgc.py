import os, shutil, StringIO
import pandas as pd
from lxml import etree
import jinja2

include_results = True

# ========
# Main html page
# ========

title = 'kinase constructs'
subtitle = 'Selected from addgene kinase collection and HIP pJP1520 plasmid libraries'
html_output_cols = ['targetID', 'DB_target_rank', 'plasmid_ID', 'plasmid_source', 'plasmid_nconflicts', 'plasmid_nextraneous_residues', 'nPDBs', 'top_pdb_ID', 'top_pdb_expr_tag', 'top_pdb_auth_score', 'top_pdb_nextraneous_residues', 'family', 'top_pdb_taxname', 'selected_construct_source', 'selected_construct_nextraneous_residues', 'caliper(ng/ul)', 'yield(ug/mL_culture)']

parser = etree.HTMLParser(remove_blank_text=True)

df = pd.DataFrame.from_csv('../expression-constructs/addgene_hip_sgc/selected-kinases.csv')

if include_results:
    results_df = pd.read_pickle('../expression-constructs/addgene_hip_sgc/results-post_expression_testing/results.p')
    orig_colnames = ['target ID', 'Conc. (ng/ul)', 'expected mg / L culture']
    results_selected_cols = pd.DataFrame()
    results_selected_cols['targetID'] = results_df['target ID']
    # Fix conversion.
    results_df['expected mg / L culture'] = results_df['Conc. (ng/ul)'] * (120.0 / 900.0)
    #
    results_selected_cols['caliper(ng/ul)'] = results_df['Conc. (ng/ul)'].map('{:.0f}'.format)
    results_selected_cols['yield(ug/mL_culture)'] = results_df['expected mg / L culture'].map('{:.1f}'.format)
    df = df.merge(results_selected_cols, on='targetID', how='left')



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
# Add href links to alignments
# ========

for tr in tbody:
    for td in tr:
        if 'HIP pJP1520' in td.text:
            targetID_td = tr.find('td[1]')
            targetID = targetID_td.text
            targetID_td.text = ''
            aelem = etree.SubElement(targetID_td, 'a')
            aelem.set('href', 'plasmid-construct-alignments/addgene_hip_sgc/' + targetID + '.html')
            aelem.text = targetID
            break
        elif 'addgene' in td.text:
            targetID_td = tr.find('td[1]')
            targetID = targetID_td.text
            targetID_td.text = ''
            aelem = etree.SubElement(targetID_td, 'a')
            aelem.set('href', 'plasmid-construct-alignments/addgene_hip_sgc/' + targetID + '.html')
            aelem.text = targetID
            break
        elif 'SGC' in td.text:
            targetID_td = tr.find('td[1]')
            targetID = targetID_td.text
            targetID_td.text = ''
            aelem = etree.SubElement(targetID_td, 'a')
            aelem.set('href', 'plasmid-construct-alignments/addgene_hip_sgc/' + targetID + '.html')
            aelem.text = targetID
            break

# ========
# Render html from template and output
# ========

table_content = '\n'.join([etree.tostring(element) for element in table_children])

env = jinja2.Environment(loader=jinja2.PackageLoader('app', 'templates'))
template = env.get_template('expr-cnstrct-custom-selection-infotext-addgene_hip_sgc.html')

with open('kinase_constructs-addgene_hip_sgc.html', 'w') as html_file:
    html_file.write( template.render(title=title, subtitle=subtitle, maintable=table_content, display_filters=False) )

# ========
# Copy html alignment files from main branch
# ========

# copy CSS
shutil.copy('stylesheets/seqlib.css', 'plasmid-construct-alignments/addgene_hip_sgc/seqlib.css')

targetIDs = list(df['targetID'])
for targetID in targetIDs:
    src_html_filepath = os.path.join('..', 'expression-constructs', 'addgene_hip_sgc', 'alignments', targetID + '.html')
    dest_html_filepath = os.path.join('plasmid-construct-alignments', 'addgene_hip_sgc', targetID + '.html')
    shutil.copy(src_html_filepath, dest_html_filepath)

