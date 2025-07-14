#!/usr/bin/env python

'''
concatenate_protein_alignments.py

Concantenate protein alignments for all "core genes" for all genomes into a single file

Filtering broadly based on criteria from Low et al. (2019) (https://doi.org/10.1038/s41564-019-0448-z)

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
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = ArgumentParser()
parser.add_argument("-r", "--reference_core_genes", dest="core_genes_infile",
                    help="core_genes.annotations.tsv file generated via identify_core_genes.py",
                    metavar="core_genes.annotations.tsv", required=True)
parser.add_argument("-a", "--alignment_files_directory", dest="alignment_files_directory",
                    help="directory containing aln.VOGID.faa files",
                    metavar="alignment_files_directory", required=True)
parser.add_argument("-o", "--output_filename", dest="outfile",
                    help="Output file name for concatenated alignment file. Default = concatenated_alignment.faa",
                    metavar='concatenated_alignment.faa', default='concatenated_alignment.faa')
args = parser.parse_args()

def concatenate_protein_alignments():
    # vog IDs
    marker_genes_df = pd.read_csv(args.core_genes_infile, sep='\t').drop('FunctionalCategory', axis=1)
    marker_genes_dict = dict(marker_genes_df.values)
    # establish df
    concat_alignments_df = pd.DataFrame(columns=['genomeID'])
    # loop throgh each alignment, create as df, join with running full df
    print('Processing individual protein alignments\r\n')
    for vog_id in marker_genes_dict.keys():
        print('Processing aln.'+vog_id+'.faa')
        aa_seqs_dict = {}
        with open(str(args.alignment_files_directory).rstrip("/") + '/aln.'+vog_id+'.faa', 'r') as read_fasta:
            for name, seq in SimpleFastaParser(read_fasta):
                aa_seqs_dict[name] = seq
        aa_seqs_df = pd.DataFrame(aa_seqs_dict.items(), columns=['geneID', 'seq_alignment'])
        aa_seqs_df['genomeID'] = aa_seqs_df['geneID'].str.replace('_\d+_VOG.*', '', regex=True)
        # if more than one gene for this vog in a genome, take the first alignment
        aa_seqs_df = aa_seqs_df.drop_duplicates(subset='genomeID', keep="first")
        # split alignment into columns
        aa_seqs_df = pd.concat([aa_seqs_df['genomeID'], aa_seqs_df['seq_alignment'].apply(lambda x: pd.Series(list(x)))], axis=1)
        aa_seqs_df = aa_seqs_df.rename(columns={c: str(vog_id)+'_pos_'+str(c) for c in aa_seqs_df.columns if c not in ['genomeID']})
        # replace all * with '-'
        aa_seqs_df = aa_seqs_df.replace('*', '-')
        ## Filter to remove columns (amino acid positions) present in less than 50% of genomes
        print('Filtering protein alignment to remove amino acid positions present in less than 50% of genomes')
        aa_seqs_df = aa_seqs_df.replace('-', np.nan)
        aa_seqs_df = aa_seqs_df.dropna(thresh=int(len(aa_seqs_df)*0.5), axis=1)
        print('Concatenating alignments')
        concat_alignments_df = concat_alignments_df.merge(aa_seqs_df, how="outer", on='genomeID')
    print('\nProcessing concatenated protein alignment\r\n')
    # sort by genomeID
    concat_alignments_df = concat_alignments_df.sort_values('genomeID').reset_index(drop=True)
    # remove genomes with <5% amino acid representation of the total alignment length
    genomes_count = len(concat_alignments_df)
    concat_alignments_df = concat_alignments_df.dropna(subset=[col for col in concat_alignments_df.columns if col != 'genomeID'], thresh=int((len(concat_alignments_df.columns)-1)*0.05), axis=0).reset_index(drop=True)
    genomes_filt_count = len(concat_alignments_df)
    print(str(genomes_count - genomes_filt_count) + ' genomes removed due to having <5% amino acid representation of the total alignment length')
    print(str(genomes_filt_count) + ' genomes remaining')
    # Calculate % missing data
    missing_data = 100 * round(concat_alignments_df[[col for col in concat_alignments_df.columns if col != 'genomeID']].isna().mean().mean(), 4)
    print('Concatenated alignment missing data = ' + str(missing_data) + '%')
    # replace all nan with '-'
    concat_alignments_df = concat_alignments_df.replace(np.nan, '-')
    # Write out concatenated alignment file
    concat_alignments_df['seq_alignment'] = concat_alignments_df[[col for col in concat_alignments_df.columns if col != 'genomeID']].astype(str).agg(''.join, axis=1)
    concat_alignments_df = concat_alignments_df[['genomeID', 'seq_alignment']]
    concat_alignments_dict = dict(concat_alignments_df.values)
    with open(args.outfile, 'w') as write_faa:
        for key,value in concat_alignments_dict.items():
            write_faa.write('>' + str(key) + '\n' + str(value) + '\n')


def main():
    print("\n--------------------\r\n")
    print("Running concatenate_protein_alignments.py\r\n")
    concatenate_protein_alignments()
    print("\nCompleted concatenate_protein_alignments.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

