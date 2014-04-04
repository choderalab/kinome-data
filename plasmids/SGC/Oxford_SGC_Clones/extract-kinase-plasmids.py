import os
import argparse
import openpyxl
import pandas as pd
from lxml import etree

# ========
# Command-line args
# ========

argparser = argparse.ArgumentParser()
argparser.add_argument('--database_path', type=str, help='Path to a TargetExplorer database XML file', required=True)
args = argparser.parse_args()

# ========
# Read in plasmid spreadsheet
# ========

plasmid_lib_filepath = os.path.join('oxford sgc clones for sbs_020513.xlsx')
wb = openpyxl.load_workbook(plasmid_lib_filepath)
ws = wb.get_sheet_by_name(name='Sheet1')
nrows = ws.get_highest_row()

# ========
# DataFrame for plasmid data
# ========

plasmid_df = {
'cloneID':[],
'HGNCSymbol':[],
'dna_seq':[],
'aa_seq':[],
'Expression system':[],
'Expression cell line':[],
'Availability comments':[],
'Special expression comments':[],
'Vector name':[],
'Antibiotic resistance':[],
'N-tag':[],
'C-tag':[],
'Vector comments':[],
'Mutations/sequence comments':[],
}

# ========
# Read in database
# ========

DB_root = etree.parse(args.database_path)

# ========
# Iterate through plasmid spreadsheet
# ========

for row in range(2, nrows-1):
    cloneID = ws.cell('C%d' % row).value
    if cloneID == None:
        continue
    HGNCSymbol = ws.cell('G%d' % row).value

    # Search for Hugo Symbol in database
    matching_DB_entry = DB_root.find('entry/HGNC/entry[@Approved_Symbol="%s"]/../..' % HGNCSymbol)
    if matching_DB_entry == None:
        continue

    plasmid_df['cloneID'].append(cloneID)
    plasmid_df['HGNCSymbol'].append(HGNCSymbol)
    plasmid_df['aa_seq'].append( ws.cell('T%d' % row).value )
    plasmid_df['dna_seq'].append( ws.cell('U%d' % row).value )
    plasmid_df['Expression system'].append( ws.cell('J%d' % row).value )
    plasmid_df['Expression cell line'].append( ws.cell('K%d' % row).value )
    plasmid_df['Availability comments'].append( ws.cell('L%d' % row).value )
    plasmid_df['Special expression comments'].append( ws.cell('M%d' % row).value )
    plasmid_df['Vector name'].append( ws.cell('N%d' % row).value )
    plasmid_df['Antibiotic resistance'].append( ws.cell('O%d' % row).value )
    plasmid_df['N-tag'].append( ws.cell('P%d' % row).value )
    plasmid_df['C-tag'].append( ws.cell('Q%d' % row).value )
    plasmid_df['Vector comments'].append( ws.cell('R%d' % row).value )
    plasmid_df['Mutations/sequence comments'].append( ws.cell('S%d' % row).value )

plasmid_df = pd.DataFrame(plasmid_df)

plasmid_df.to_csv('plasmid-data.csv')
with open('plasmid-data.txt', 'w') as plasmid_data_txt_file:
    plasmid_data_txt_file.write(plasmid_df.to_string())
