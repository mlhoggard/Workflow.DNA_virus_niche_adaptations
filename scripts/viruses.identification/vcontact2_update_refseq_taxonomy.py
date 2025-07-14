#!/usr/bin/env python

'''
vcontact2_update_refseq_taxonomy.py

Update the taxonomy assignments in vConTACT2's genome_by_genome_overview.csv file based on reconcilation with ICTV taxonomy. Takes genome_by_genome_overview.csv and viralRefSeq_ICTV_reconciled_taxonomy.tsv (generated via the script ictv_reconcile_refseq_taxonomy.py) as inputs.

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
parser.add_argument("-i", "--input", dest="genome_by_genome",
                    help="genome_by_genome_overview.csv file generated via vConTACT2",
                    metavar="genome_by_genome_overview.csv", required=True)
parser.add_argument("-t", "--taxonomy_update", dest="reconciled_taxonomy",
                    help="File of taxonomy reconcilations for viralRefSeq references based on latest ICTV updates (e.g. viralRefSeq_ICTV_reconciled_taxonomy.tsv generated via the script ictv_reconcile_refseq_taxonomy.py)",
                    metavar="viralRefSeq_ICTV_reconciled_taxonomy.tsv", required=True)
parser.add_argument("-o", "--output", dest="outfile",
                    help="(Optional) Output name for genome_by_genome_overview.csv with updated taxonomy. Default = genome_by_genome_overview.tax_update.csv",
                    metavar='genome_by_genome_overview.tax_update.csv', default='genome_by_genome_overview.tax_update.csv')
args = parser.parse_args()

def update_taxonomy():
    # genome_by_genome file
    df = pd.read_csv(args.genome_by_genome,
                    dtype={'Genome': str, 'Order': str, 'Family': str, 'Genus': str, 'VC': str,
                            'VC Status': str, 'Size': float, 'VC Subcluster': str, 'VC Subcluster Size': float,
                            'Quality': float, 'Adj P-value': float, 'Topology Confidence Score': float,
                            'Genera in VC': float, 'Families in VC': float, 'Orders in VC': float,
                            'Genus Confidence Score': float})
    # updated taxonomy
    tax_ref = pd.read_csv(args.reconciled_taxonomy, sep='\t', dtype=str, keep_default_na=False)
    tax_ref = tax_ref[['refseq_Organism_Name', 'refseq_Species', 'refseq_Genus', 'refseq_Family', 'refseq_Molecule_type',
                        'ICTV_Realm', 'ICTV_Subrealm', 'ICTV_Kingdom', 'ICTV_Subkingdom', 'ICTV_Phylum', 'ICTV_Subphylum',
                        'ICTV_Class', 'ICTV_Subclass', 'ICTV_Order', 'ICTV_Suborder', 'ICTV_Family', 'ICTV_Subfamily', 'ICTV_Genus']].drop_duplicates()
    # merge
    df_update = pd.merge(df, tax_ref, how="left", left_on='Genome', right_on='refseq_Organism_Name')
    # select taxonomy ranks to keep (family, genus, species from refseq update; rest from ranks filled in from ICTV)
    for tax_variable in ['Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', 'Subfamily']:
        df_update[tax_variable] = df_update['ICTV_'+tax_variable]
    for tax_variable in ['Family', 'Genus', 'Species']:
        df_update[tax_variable] = df_update['refseq_'+tax_variable]
    # drop extra taxonomy columns and reorder
    df_update = df_update[['Genome', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 'Genus', 'Species',
                            'preVC', 'VC Status', 'VC', 'VC Size', 'Quality', 'Adjusted P-value', 'VC Avg Distance', 'Topology Confidence Score', 'Genus Confidence Score',
                            'VC Kingdoms', 'VC Phyla', 'VC Classes', 'VC Orders', 'VC Families', 'VC Genera']]
    # write up updated genome_by_genome file
    df_update.to_csv(args.outfile, index=False)

def main():
    print("\n--------------------\r\n")
    print("Running vcontact2_update_refseq_taxonomy.py\r\n")
    update_taxonomy()
    print("\nCompleted vcontact2_update_refseq_taxonomy.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

