# Kinase expression panel results table

Pages are included for a few different sets of results.
The expression panel which was actually tested is named `addgene_hip_sgc`. The relevant HTML page is at:
`kinase_constructs-addgene_hip_sgc.html`

This file is generated by running:
`pythonscripts/gen-html-plasmids-addgene_hip_sgc.py`

Note that input data filepaths are currently hardcoded (lines 18, 21, and 104). These may need to be changed.

`plasmid-construct-alignments` contains HTML sequence alignments, which can be found in the repo [kinase-ecoli-expression-panel](https://github.com/choderalab/kinase-ecoli-expression-panel)

The Python script uses jinja templating. Templates are stored in `app/templates`

The table feature is implemented using DataTables js package.
DataTables configuration can be found in this template file:
`app/templates/expression-construct-table.html`
