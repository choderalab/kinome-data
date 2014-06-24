A panel of kinase constructs for expression testing.

Manifest
--------

* scripts/select-kinase-constructs.py
  * takes plasmid and PDB construct data (from ../plasmids and ../PDB-constructs), and conducts a custom ranking and filtering protocol to select a suitable panel of expression constructs
  * outputs:
    * selected-kinases-sgc\_and\_hip.csv - output data for the custom selection of kinase constructs
    * selected-kinases-sgc\_and\_hip.txt - ascii table showing the selected constructs with some basic data
* scripts/mk\_spreadsheet.py
  * generates a spreadsheet with the data required to generate primers and start an expression project
  * outputs:
    * selected-kinases-sgc\_and\_hip.xlsx - spreadsheet
    * selected-kinases-sgc\_and\_hip.fa - for each kinase, contains an alignment of the UniProt seq, the selected plasmid, and the selected constructs (if different from the plasmid)
* scripts/mk\_alignments.py
  * generates an html alignment (in alignments/) for each kinase, containing the UniProt seq, all matching plasmids (sorted), and all PDB constructs (sorted)
* Spreadsheets and other material for exploring possible construct synthesis by gen9
