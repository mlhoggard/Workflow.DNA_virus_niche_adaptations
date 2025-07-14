#!/usr/bin/env python

'''
featurecounts_make_feature_table.py

Generate an SAF formatted feature table from DRAM (prokaryotes) and/or DRAMv (viruses) annotations.tsv, rrna.tsv, and trna.tsv data to assign mapped read counts to features (e.g. genes) via featureCounts.

Note: this script also trims "-cat_n" from virus contig/genome IDs (introduced by VirSorter2 during the prep for DRAMv step) to keep IDs consistent with other analyses.

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
parser.add_argument("-a", "--dram_annotations", dest="dram_annotations",
                    help="annotations.tsv file generated via DRAM annotate",
                    metavar="dram_annotations.tsv", default=None)
parser.add_argument("-t", "--dram_trnas", dest="dram_trnas",
                    help="trnas.tsv file generated via DRAM annotate",
                    metavar="dram_trnas.tsv", default=None)
parser.add_argument("-r", "--dram_rrnas", dest="dram_rrnas",
                    help="rrnas.tsv file generated via DRAM annotate",
                    metavar="dram_rrnas.tsv", default=None)
parser.add_argument("-x", "--dramv_annotations", dest="dramv_annotations",
                    help="annotations.tsv file generated via DRAMv annotate",
                    metavar="dramv_annotations.tsv", default=None)
parser.add_argument("-y", "--dramv_trnas", dest="dramv_trnas",
                    help="trnas.tsv file generated via DRAMv annotate",
                    metavar="dramv_trnas.tsv", default=None)
parser.add_argument("-z", "--dramv_rrnas", dest="dramv_rrnas",
                    help="rrnas.tsv file generated via DRAMv annotate",
                    metavar="dramv_rrnas.tsv", default=None)
parser.add_argument("-o", "--output", dest="outfile",
                    help="(Optional) Output file name for SAF feature table. Default = gene_coords.SAF",
                    metavar='output_filename.SAF', default='gene_coords.SAF')
args = parser.parse_args()

def read_annotations(data):
    annot_df = pd.read_csv(data, sep='\t')[['Unnamed: 0', 'scaffold', 'start_position', 'end_position', 'strandedness']].rename(columns = {'Unnamed: 0': 'GeneID', 'scaffold': 'Chr', 'start_position': 'Start', 'end_position': 'End', 'strandedness': 'Strand'})
    # remove "-cat_n" and any spaces from IDs
    annot_df['Chr'] = annot_df['Chr'].str.replace(r'-cat_\d', r'', regex=True).str.replace(r' ', r'')
    annot_df['GeneID'] = annot_df['GeneID'].str.replace(r'-cat_\d', r'', regex=True).str.replace(r' ', r'')
    # update strand to +/-
    annot_df['Strand'] = annot_df['Strand'].astype(str).str.replace(r'-1', r'-')
    annot_df['Strand'] = annot_df['Strand'].astype(str).str.replace(r'1', r'+')
    annot_df = annot_df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
    return annot_df

def read_rrna(data):
    rrna_df = pd.read_csv(data, sep='\t')[['scaffold', 'begin', 'end', 'strand', 'type']].rename(columns = {'scaffold': 'Chr', 'begin': 'Start', 'end': 'End', 'strand': 'Strand'})
    # remove "-cat_n" and any spaces from IDs
    rrna_df['Chr'] = rrna_df['Chr'].str.replace(r'-cat_\d', r'', regex=True).str.replace(r' ', r'')
    # add 'GeneID'
    rrna_df['GeneID'] = rrna_df['Chr'] + '_' + rrna_df['type'].str.replace(r' ', r'_')
    rrna_df = rrna_df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
    return rrna_df

def read_trna(data):
    trna_df = pd.read_csv(data, sep='\t')[['Name', 'tRNA #', 'Begin', 'End', 'Type', 'Codon']].rename(columns = {'Name': 'Chr', 'Begin': 'start', 'End': 'end'})
    # remove "-cat_n" and any spaces from IDs
    trna_df['Chr'] = trna_df['Chr'].str.replace(r'-cat_\d', r'', regex=True).str.replace(r' ', r'')
    # add 'GeneID'
    trna_df['GeneID'] = trna_df['Chr'] + '_tRNA_' + trna_df['tRNA #'].astype(str) + '_' + trna_df['Type'] + '_' + trna_df['Codon']
    # add strand column based on start and end coords and correcting order of coords for -ve strand entries
    trna_df['Strand'] = np.select([(trna_df['start'] < trna_df['end']),
                                (trna_df['start'] > trna_df['end'])
                                ], ['+', '-'], default=np.nan)
    trna_df['Start'] = np.select([(trna_df['Strand'] == '+'), (trna_df['Strand'] == '-')], [trna_df['start'].astype(str), trna_df['end'].astype(str)], default=np.nan)
    trna_df['End'] = np.select([(trna_df['Strand'] == '+'), (trna_df['Strand'] == '-')], [trna_df['end'].astype(str), trna_df['start'].astype(str)], default=np.nan)
    trna_df = trna_df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
    return trna_df

def main():
    print("\n--------------------\r\n")
    print("Running featurecounts_make_feature_table.py\r\n")
    res = []
    if args.dram_annotations:
        print('Reading file '+args.dram_annotations+'\n')
        res.append(read_annotations(args.dram_annotations))
    if args.dram_rrnas:
        print('Reading file '+args.dram_rrnas+'\n')
        res.append(read_rrna(args.dram_rrnas))
    if args.dram_trnas:
        print('Reading file '+args.dram_trnas+'\n')
        res.append(read_trna(args.dram_trnas))
    if args.dramv_annotations:
        print('Reading file '+args.dramv_annotations+'\n')
        res.append(read_annotations(args.dramv_annotations))
    if args.dramv_rrnas:
        print('Reading file '+args.dramv_rrnas+'\n')
        res.append(read_rrna(args.dramv_rrnas))
    if args.dramv_trnas:
        print('Reading file '+args.dramv_trnas+'\n')
        res.append(read_trna(args.dramv_trnas))
    if len(res) == 0:
        print('Error: no valid input files provided. Exiting.')
        exit
    # concatenate results
    df = pd.concat(res, ignore_index=True, sort=False).reset_index(drop=True).sort_values(['GeneID', 'Chr', 'Start'], ascending=[True, True, True]).reset_index(drop=True)
    # overwrite any instances of Start < 0 (overwrite as 0)
    df['Start'][df['Start'].astype(int) < 0] = '0'
    # Write out SAF
    df.to_csv(args.outfile, sep='\t', index=False)
    print("\nCompleted featurecounts_make_feature_table.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

