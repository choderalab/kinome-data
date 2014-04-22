import pandas as pd
import lxml.html
from lxml import etree
from lxml.html import builder as E
import StringIO
import jinja2

title = '96 kinases'
subtitle = 'Selected from SGC and HIP pJP1520 plasmid libraries'
html_output_cols = ['cloneID', 'UniProt_entry_name', 'nconflicts_target_domain_region', 'nextraneous_plasmid_residues', 'nPDBs', 'expr_tag', 'auth_score', 'DB_target_rank']

parser = etree.HTMLParser(remove_blank_text=True)

df = pd.DataFrame.from_csv('../PDB_construct_selection/96-kinases-sgc_and_hip.csv')

pd_html = df.to_html(columns=html_output_cols)
html = etree.parse(StringIO.StringIO(pd_html), parser).getroot()
table_children = html.find('body/table').getchildren()
table_content = '\n'.join([etree.tostring(element) for element in table_children])

env = jinja2.Environment(loader=jinja2.PackageLoader('app', 'templates'))
template = env.get_template('expression-construct-table.html')

with open('96-kinases-sgc_and_hip.html', 'w') as html_file:
    html_file.write( template.render(title=title, subtitle=subtitle, maintable=table_content) )

