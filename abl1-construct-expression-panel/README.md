Three "control" constructs were chosen, which had known expression data in the
literature. Five replicates of each were included in the plate.

Then 9 construct start points and 9 end points were chosen and combined to
produce 81 constructs. These were chosen to be relatively close to the
literature start/end points. The aim was to be able to study how expression
levels can be affected by single-residue start/end shifts, or by whether
start/end points are located at structured or unstructured regions.


* construct-variants.txt - notes on the selected construct variants

* python scripts/align-pdb-constructs.py --target ABL1\_HUMAN\_D0'
  * produced the HTML alignment 'pdb\_seqs.html'
  * pdb\_seqs.html - Abl1 PDB constructs aligned against the UniProt canonical isoform kinase domain
* abl1-macrolab-dna\_seq.txt - DNA sequence for the Harvard HIP Abl1 plasmid in stock at the MacroLab
* python scripts/gen-construct-variants.py
  * produced abl1-constructs.xlsx and abl1-constructs-data.txt
  * abl1-constructs.xlsx - Excel spreadsheet containing dna and aa sequences, sent to the MacroLab
  * abl1-constructs-data.csv - same data in csv format
* view-2fo0.pml - PyMOL script
* results\_from\_MacroLab/ - results of the expression tests carried out by the MacroLab
