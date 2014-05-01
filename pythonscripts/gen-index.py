import os, StringIO
import jinja2

# ========
# Render html from template and output
# ========

env = jinja2.Environment(loader=jinja2.PackageLoader('app', 'templates'))
template = env.get_template('index.html')

with open('index.html', 'w') as html_file:
    html_file.write( template.render(head='', content='') )

