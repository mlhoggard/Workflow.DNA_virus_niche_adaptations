#!/usr/bin/env python

'''
dramv_compile_summary_table.py

Compile DRAMv annotations, trna, rrna, and AMG (distill) outputs, together with gene coords and contig (genome) lengths into one table for downstream use.

MIT License

Copyright (c) 2025 Michael Hoggard

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from argparse import ArgumentParser
import pandas as pd
import numpy as np
from glob import glob
import os
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = ArgumentParser()
parser.add_argument("-a", "--dramv_annotations", dest="dramv_annotations",
                    help="annotations.tsv file generated via DRAMv annotate",
                    metavar="annotations.tsv", required=True)
parser.add_argument("-d", "--dramv_distill_amg_summary", dest="dramv_distill_amg_summary",
                    help="amg_summary.tsv generated via DRAMv distill",
                    metavar="amg_summary.tsv", required=True)
parser.add_argument("-t", "--dramv_trnas", dest="dramv_trnas",
                    help="trnas.tsv file generated via DRAMv annotate",
                    metavar="trnas.tsv", required=True)
parser.add_argument("-r", "--dramv_rrnas", dest="dramv_rrnas",
                    help="rrnas.tsv file generated via DRAMv annotate",
                    metavar="rrnas.tsv", required=True)
parser.add_argument("-f", "--viruses_fna", dest="viruses_fna",
                    help="viruses fna file used as input into DRAMv",
                    metavar="viruses.fna", required=True)
parser.add_argument("-o", "--output", dest="outfile",
                    help="(Optional) Output file name for DRAMv results summary table. Default = DRAMv.summary_table.tsv",
                    metavar='output_filename.tsv', default='DRAMv.summary_table.tsv')
args = parser.parse_args()

def compile_summary_table():
    ## dramv annotations
    annot_df = pd.read_csv(args.dramv_annotations, sep='\t', usecols=lambda x: x not in ['fasta', 'Unnamed: 0'],
                        dtype={
                        'scaffold': str, 'gene_position': int,
                            'start_position': int, 'end_position': int, 'strandedness': int, 'rank': str,
                            'kegg_genes_id': str, 'ko_id': str, 'kegg_hit': str, 'kegg_RBH': str, 'kegg_identity': float, 'kegg_bitScore': float, 'kegg_eVal': float,
                            'viral_id': str, 'viral_hit': str, 'viral_RBH': str, 'viral_identity': float, 'viral_bitScore': float, 'viral_eVal': float,
                            'peptidase_id': str, 'peptidase_family': str, 'peptidase_hit': str, 'peptidase_RBH': str, 'peptidase_identity': float, 'peptidase_bitScore': float, 'peptidase_eVal': float,
                            'pfam_hits': str, 'cazy_id': str, 'cazy_hits': str, 'cazy_subfam_ec': str, 'cazy_best_hit': str,
                            'vogdb_id': str, 'vogdb_hits': str, 'vogdb_categories': str,
                            'heme_regulatory_motif_count': int, 'virsorter_category': float, 'auxiliary_score': float, 'is_transposon': str, 'amg_flags': str}).rename(
        {'scaffold': 'genomeID', 'auxiliary_score': 'virsorter_auxiliary_score', 'rank': 'virsorter_rank'}, axis=1)
    annot_df["genomeID"] = annot_df["genomeID"].str.replace(r"-cat_.*", "", regex=True)
    annot_df['geneID'] = annot_df['genomeID'] + '_' + annot_df['gene_position'].astype(str)
    ## AMGs
    amg_df = pd.read_csv(args.dramv_distill_amg_summary, sep='\t', usecols=lambda x: x not in ['auxiliary_score', 'amg_flags']).rename({'scaffold': 'genomeID'}, axis=1)
    amg_df['gene_position'] = amg_df['gene'].str.replace(r".*-cat_\d_", "", regex=True)
    amg_df["genomeID"] = amg_df["genomeID"].str.replace(r"-cat_.*", "", regex=True)
    amg_df['geneID'] = amg_df['genomeID'] + '_' + amg_df['gene_position'].astype(str)
    amg_df.drop(['gene_position', 'gene'], axis=1, inplace=True)
    amg_df['AMG'] = 'AMG'
    amg_df.columns = [f'AMG_{i}' if i not in ['geneID', 'genomeID', 'AMG'] else f'{i}' for i in amg_df.columns]
    annot_df = pd.merge(annot_df, amg_df, how='outer', on=['geneID', 'genomeID'])
    ## tRNAs
    trna_df = pd.read_csv(args.dramv_trnas, sep='\t', usecols=lambda x: x not in ['fasta']).rename(
        {'Name': 'genomeID', 'tRNA #': 'tRNA_count', 'Begin': 'start', 'End': 'end', 'Type': 'tRNA_species', 'Codon': 'tRNA_codon', 'Score': 'tRNA_score', 'Note': 'tRNA_Notes'}, axis=1)
    trna_df["genomeID"] = trna_df["genomeID"].str.replace(r"-cat_.*", "", regex=True)
    trna_df['geneID'] = trna_df['genomeID'] + '_tRNA_' + trna_df['tRNA_count'].astype(str)
    trna_df['strandedness'] = np.select([(trna_df['start'] < trna_df['end']),
                                (trna_df['start'] > trna_df['end'])
                                ], ['+', '-'], default=np.nan)
    trna_df['start_position'] = np.select([(trna_df['strandedness'] == '+'), (trna_df['strandedness'] == '-')], [trna_df['start'].astype(str), trna_df['end'].astype(str)], default=np.nan)
    trna_df['end_position'] = np.select([(trna_df['strandedness'] == '+'), (trna_df['strandedness'] == '-')], [trna_df['end'].astype(str), trna_df['start'].astype(str)], default=np.nan)
    trna_df.drop(['start', 'end', 'tRNA_count'], axis=1, inplace=True)
    ## rRNAs
    rrna_df = pd.read_csv(args.dramv_rrnas, sep='\t', usecols=lambda x: x not in ['fasta']).rename(
        {'scaffold': 'genomeID', 'begin': 'start_position', 'end': 'end_position', 'strand': 'strandedness', 'type': 'rRNA_type', 'e-value': 'rRNA_e_value', 'note': 'rRNA_Notes'}, axis=1)
    rrna_df["genomeID"] = rrna_df["genomeID"].str.replace(r"-cat_.*", "", regex=True)
    rrna_df['rRNA_count']=rrna_df.groupby('genomeID').cumcount()+1
    rrna_df['geneID'] = rrna_df['genomeID'] + '_rRNA_' + rrna_df['rRNA_count'].astype(str)
    rrna_df.drop('rRNA_count', axis=1, inplace=True)
    ## genome lengths
    fasta_contigs = {}
    with open(args.viruses_fna, 'r') as read_fasta:
        for name, seq in SimpleFastaParser(read_fasta):
            fasta_contigs[name] = len(seq)
    genome_lengths_df = pd.DataFrame(fasta_contigs.items(), columns=['genomeID', 'genome_sequence_length'])
    ## Merge annotations together, Modify strandedness to be + or -, sort by contigID and start_position, add contig lengths, reorder columns, and write out
    # merge, modify strandedness to be + or -, and sort
    combined_df = pd.concat([annot_df, trna_df, rrna_df], sort=False, ignore_index=True).merge(genome_lengths_df, on='genomeID', how='left')
    combined_df['strandedness'] = combined_df['strandedness'].astype(str).replace({'-1': '-'}).replace({'1': '+'})
    combined_df["start_position"] = pd.to_numeric(combined_df["start_position"], errors='coerce')
    combined_df["end_position"] = pd.to_numeric(combined_df["end_position"], errors='coerce')
    combined_df = combined_df.sort_values(['genomeID', 'start_position'], ascending=[True, True]).reset_index(drop=True)
    # Re-order columns
    cols_first = ['genomeID', 'genome_sequence_length', 'gene_position', 'geneID', 'start_position', 'end_position', 'strandedness', 'AMG']
    cols_kegg = [col for col in combined_df.columns if 'ko_' in col or 'kegg_' in col]
    cols_pfam = [col for col in combined_df.columns if 'pfam' in col]
    cols_peptidase = [col for col in combined_df.columns if 'peptidase_' in col]
    cols_cazy = [col for col in combined_df.columns if 'cazy_' in col]
    cols_viral = [col for col in combined_df.columns if 'viral_' in col]
    cols_vog = [col for col in combined_df.columns if 'vog' in col]
    cols_trna = [col for col in combined_df.columns if 'tRNA_' in col]
    cols_rrna = [col for col in combined_df.columns if 'rRNA_' in col]
    cols_amg_other = [col for col in combined_df.columns if 'is_transposon' in col or 'virsorter' in col]
    cols_amg = [col for col in combined_df.columns if 'AMG_' in col or 'amg_' in col]
    cols_rest = [col for col in combined_df.columns if col not in cols_first+cols_kegg+cols_pfam+cols_peptidase+cols_cazy+cols_viral+cols_vog+cols_trna+cols_rrna+cols_amg+cols_amg_other]
    combined_df = combined_df[cols_first+cols_kegg+cols_pfam+cols_peptidase+cols_cazy+cols_viral+cols_vog+cols_trna+cols_rrna+cols_rest+cols_amg_other+cols_amg]
    # Write out
    combined_df.to_csv(args.outfile, sep='\t', index=False)

def main():
    print("\n--------------------\r\n")
    print("Running dramv_compile_summary_table.py\r\n")
    compile_summary_table()
    print("\nCompleted dramv_compile_summary_table.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

