Analysis of SGC Oxford expression clones
========================================

Manifest
--------

* oxford sgc clones for sbs\_020513.xlsx - plasmid data spreadsheet downloaded from Source Bioscience website:
  * http://www.lifesciences.sourcebioscience.com/clone-products/mammalian/genomic-clones-others/structural-genomics-consortium-expression-clones.aspx
* extract-kinase-plasmids.py - matches spreadsheet plasmid entries against a TargetExplorer database and outputs data to 'plasmid-data.txt' (pretty-formatted) and 'plasmid-data.csv'
* aln-against-UniProt-seq.py - each of the plasmids in 'plasmid-data.csv' is aligned against the matching UniProt sequence for comparison. Output is to 'aln.html' (pretty-formatted) and 'aln.txt' (not-so-pretty-formatted).
