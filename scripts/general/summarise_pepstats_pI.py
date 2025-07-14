#!/usr/bin/env python

'''
summarise_pepstats_pI.py

Extract and summarise pI results generated via pepstats.

MIT License

Copyright (c) 2025 Michael Hoggard

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from argparse import ArgumentParser
import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re
import os
from Bio import SeqIO
from itertools import groupby, chain

parser = ArgumentParser()
parser.add_argument("-i", "--input_file", dest="in_filename",
                    help="Input file generated via pepstats.",
                    metavar="results.pepstats", required=True)
parser.add_argument("-s", "--sample_id", dest="sample_id",
                    help="ID string (optional). Adds sample ID column to summary table (useful for compiling results from multiple samples downstream)",
                    metavar="sampleID", default=None)
parser.add_argument("-o", "--output_filename", dest="out_filename",
                    help="Filename for results summary table (tsv format). Default = pepstats.summaryTable.tsv.",
                    metavar='output_filename.tsv', default='pepstats.summaryTable.tsv')
args = parser.parse_args()

def get_gene_pepstats(file):
    with open(file) as f:
        grps = groupby(f, key=lambda x: x.lstrip().startswith("PEPSTATS of"))
        for k, v in grps:
            if k:
                yield chain([next(v)], (next(grps)[1]))

def summarise_results():
    iep_dict = {}
    for gene_pepstats in get_gene_pepstats(args.in_filename):
        gene_pepstats_tmp = (list(gene_pepstats))
        geneID = [x for x in gene_pepstats_tmp if 'PEPSTATS of' in x][0].split(' ')[2]
        iso_point = [x for x in gene_pepstats_tmp if 'Isoelectric Point' in x][0].split(' ')[3].rstrip()
        iep_dict[geneID] = iso_point
    # compile df
    iep_df = pd.DataFrame(iep_dict.items(), columns=['geneID', 'isoelectric_point'])
    if args.sample_id:
        iep_df['sampleID'] = args.sample_id
    # write out
    iep_df.to_csv(args.out_filename, sep='\t', index=False)


def main():
    print("\n--------------------\r\n")
    print("Running summarise_pepstats_pI.py\r\n")
    summarise_results()
    print("\nCompleted summarise_pepstats_pI.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

