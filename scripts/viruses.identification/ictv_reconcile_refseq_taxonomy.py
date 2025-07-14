#!/usr/bin/env python

'''
ictv_reconcile_refseq_taxonomy.py

Propogate upper ranks of viralRefSeq taxonomy assignments based on ICTV taxonomy table.

Note: This was written to enable correcting the taxonomy predictions associated with vConTACT results (the version of which used an older RefSeq version with outdated taxonomy). Script written based on viralRefSeq v223 viralRefSeq_metadata.csv, and ICTV_Master_Species_List_2022_MSL38.v3.xlsx, and may break if there are formatting changes in other versions of these databases.

MIT License

Copyright (c) 2025 Michael Hoggard

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import re

parser = ArgumentParser()
parser.add_argument("-i", "--ictv_species_list", dest="ictv_species_list",
                    help="ICTV_Master_Species_List.xlsx downloaded from ICTV",
                    metavar="ICTV_Master_Species_List.xlsx", required=True)
parser.add_argument("-r", "--refseq_metadata", dest="refseq_metadata",
                    help="viralRefSeq_metadata.csv downloaded from RefSeq",
                    metavar="viralRefSeq_metadata.csv", required=True)
parser.add_argument("-o", "--output", dest="outfile",
                    help="(Optional) Output name for reconciled taxonomy tsv file. Default = viralRefSeq_ICTV_reconciled_taxonomy.tsv",
                    metavar='output_filename.tsv', default='viralRefSeq_ICTV_reconciled_taxonomy.tsv')
args = parser.parse_args()

def reconcile_taxonomy():
    ## ictv taxonomy sheet
    ictv_df = pd.read_excel(args.ictv_species_list, sheet_name='MSL', index_col=None)
    ictv_df = ictv_df[[x for x in ictv_df.columns if x not in ['Sort', 'Last Change', 'MSL of Last Change', 'Proposal for Last Change ', 'Taxon History URL']]].rename(columns={'Genome Composition': 'Genome_Composition'})
    ictv_df = ictv_df.add_prefix('ICTV_')
    ## refseq_update_df
    refseq_update_df = pd.read_csv(args.refseq_metadata)
    # Add "Unclassified" and/or propagate taxonomy for empty ranks
    refseq_update_df['Family'] = refseq_update_df['Family'].fillna('f__Unclassified')
    refseq_update_df['Genus'] = refseq_update_df['Genus'].fillna(refseq_update_df['Family']+';g__Unclassified')
    refseq_update_df['Species'] = refseq_update_df['Species'].fillna(refseq_update_df['Genus']+';s__Unclassified')
    refseq_update_df = refseq_update_df.add_prefix('refseq_')
    ## join ictv_df to refseq to add higher ranks.
    # with genus classification (join by genus)
    refseq_update_df_sub1 = refseq_update_df[~refseq_update_df['refseq_Genus'].str.contains('Unclassified')]
    ictv_genus = ictv_df[['ICTV_Realm', 'ICTV_Subrealm', 'ICTV_Kingdom', 'ICTV_Subkingdom', 'ICTV_Phylum', 'ICTV_Subphylum', 'ICTV_Class', 'ICTV_Subclass', 'ICTV_Order', 'ICTV_Suborder', 'ICTV_Family', 'ICTV_Subfamily', 'ICTV_Genus']]
    ictv_genus = ictv_genus.drop_duplicates()
    ictv_genus = ictv_genus[ictv_genus['ICTV_Genus'].notnull()]
    tax_update_sub1 = pd.merge(refseq_update_df_sub1,
                            ictv_genus,
                            how="left", left_on='refseq_Genus', right_on='ICTV_Genus')
    # rest (join by family instead)
    refseq_update_df_sub2 = refseq_update_df[refseq_update_df['refseq_Genus'].str.contains('Unclassified')]
    # ictv unique taxonomy down to family
    ictv_family = ictv_df[['ICTV_Realm', 'ICTV_Subrealm', 'ICTV_Kingdom', 'ICTV_Subkingdom', 'ICTV_Phylum', 'ICTV_Subphylum', 'ICTV_Class', 'ICTV_Subclass', 'ICTV_Order', 'ICTV_Suborder', 'ICTV_Family']]
    ictv_family = ictv_family.drop_duplicates()
    ictv_family = ictv_family[ictv_family['ICTV_Family'].notnull()]
    tax_update_sub2 = pd.merge(refseq_update_df_sub2,
                            ictv_family,
                            how="left", left_on='refseq_Family', right_on='ICTV_Family')
    tax_update_sub2['ICTV_Subfamily'] = 'Unclassified'
    tax_update_sub2['ICTV_Genus'] = 'Unclassified'
    # concatenate together
    tax_update = pd.concat([tax_update_sub1,tax_update_sub2],ignore_index=True)
    # update refseq_Organism_Name to replace spaces with '~' (to match vConTACT2 formatting)
    tax_update['refseq_Organism_Name'] = tax_update['refseq_Organism_Name'].str.replace(' ', '~').copy()
    # write out viralRefSeq_ICTV_reconciled_taxonomy.tsv
    tax_update.to_csv(args.outfile, index=False, sep='\t')

def main():
    print("\n--------------------\r\n")
    print("Running ictv_reconcile_refseq_taxonomy.py\r\n")
    reconcile_taxonomy()
    print("\nCompleted ictv_reconcile_refseq_taxonomy.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

