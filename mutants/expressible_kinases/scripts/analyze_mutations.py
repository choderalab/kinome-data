import os
from openpyxl import load_workbook
import pandas as pd
from targetexplorer.flaskapp import models, db

script_dir = os.path.dirname(__file__)
expression_panel_kinases_workbook_filename = os.path.join(script_dir, '..', '..', '..', 'expression-constructs', 'addgene_hip_sgc', 'selected-kinases.xlsx')
expressible_kinases_list_filename = os.path.join(script_dir, '..', 'expressible_kinases.csv')
out_csv_filepath = os.path.join(script_dir, '..', 'mutation_data.csv')
out_txt_filepath = os.path.join(script_dir, '..', 'mutation_data.txt')

expression_panel_kinases_worksheet = load_workbook(expression_panel_kinases_workbook_filename)['kinase constructs']
# 1-based
construct_aa_starts = {
    expression_panel_kinases_worksheet.range('A{0}'.format(r+2)).value: expression_panel_kinases_worksheet.range('H{0}'.format(r+2)).value for r in range(96)
}
construct_aa_ends = {
    expression_panel_kinases_worksheet.range('A{0}'.format(r+2)).value: expression_panel_kinases_worksheet.range('I{0}'.format(r+2)).value for r in range(96)
}

expressible_kinases_df = pd.read_csv(expressible_kinases_list_filename)
expressible_kinases_list_colnames = list(expressible_kinases_df.columns)
expressible_targetids = [t for t in expressible_kinases_df['targetid']]

expressible_uniprot_domain_rows = [row for row in models.UniProtDomain.query.all() if row.targetid in expressible_targetids]

cbioportal_mutation_table = models.CbioportalMutation.__table__
cbioportal_mutation_colnames = [c.name for c in cbioportal_mutation_table.columns]

data_dict = {colname: [] for colname in expressible_kinases_list_colnames + cbioportal_mutation_colnames + ['construct_aa_start', 'construct_aa_end']}

for uniprot_domain in expressible_uniprot_domain_rows:
    targetid = uniprot_domain.targetid
    mutation_rows = uniprot_domain.cbioportal_mutations.all()
    for mutation in mutation_rows:
        for colname in expressible_kinases_list_colnames:
            data_dict[colname].append(
                expressible_kinases_df[expressible_kinases_df['targetid'] == targetid][colname].values[0]
            )
        for colname in cbioportal_mutation_colnames:
            data_dict[colname].append(
                getattr(mutation, colname)
            )
        data_dict['construct_aa_start'].append(construct_aa_starts[targetid])
        data_dict['construct_aa_end'].append(construct_aa_ends[targetid])

df = pd.DataFrame(data_dict)

df.drop([
    'id',
    'dbentry_id',
    'cbioportal_case_id',
    'uniprot_domain_id',
    'crawl_number',
    'cbioportal_aa_change_string',
    'reference_dna_allele',
    'variant_dna_allele',
    'chromosome_index',
    'chromosome_startpos',
    'chromosome_endpos',
], axis=1, inplace=True)

out_columns = [
    'targetid',
    'conc_(ng/ul)',
    'expected_conc(mg/l)',
    'construct_aa_start',
    'construct_aa_end',
    'type',
    'oncotator_aa_pos',
    'oncotator_reference_aa',
    'oncotator_variant_aa',
    'validation_status',
    'functional_impact_score',
    'mutation_origin',
]

df.to_csv(out_csv_filepath, columns=out_columns)

with open(out_txt_filepath, 'w') as out_txt_file:
    out_txt_file.write(df.to_string(columns=out_columns))
