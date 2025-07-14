#!/usr/bin/env python

'''
Filter virus genomes fna file to retain viruses >= completeness_threshold based on CheckV results.

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

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="fna_in",
                    help="viral genomes fna file to be filtered by completeness threshold",
                    metavar="viruses.fna", required=True)
parser.add_argument("-c", "--checkv_quality_summary", dest="checkv_quality_summary",
                    help="quality_summary.tsv file generated via CheckV",
                    metavar="checkv.quality_summary.tsv", required=True)
parser.add_argument("-t", "--completeness_threshold", dest="completeness_threshold",
                    help="retain genomes >= this completeness percentage threshold (based on CheckV completeness prediction). Default = 50",
                    metavar="completeness_threshold", default=50)
parser.add_argument("-o", "--output", dest="outfile",
                    help="Output file name for filtered virus genomes fna file. Default = <input_filename>.completeness_<completeness_threshold>.fna",
                    metavar='concatenated_alignment.faa', default=None)
args = parser.parse_args()

def filter_viruses_by_completeness():
    # Read in checkV stats
    checkv_df = pd.read_csv(args.checkv_quality_summary, sep='\t', usecols=['contig_id', 'contig_length', 'completeness', 'completeness_method'])
    checkv_df['completeness'] = checkv_df['completeness'].fillna(0.0).astype(float)
    # subset by completeness threshold
    vOTUs_subset_IDs = checkv_df[checkv_df['completeness'] >= float(args.completeness_threshold)].sort_values(by=['contig_length'], ascending=[False])['contig_id'].values.tolist()
    # generate output filename
    if args.outfile:
        out_filename = args.outfile
    else:
        out_filename = str(args.fna_in).rsplit(".", 1)[0] + '.completeness_' + str(args.completeness_threshold) + '.fna'
    # Write new multifasta files of genomes
    with open(args.fna_in, 'r') as read_fna:
        with open(out_filename, 'w') as write_fna:
            for name, seq in SimpleFastaParser(read_fna):
                if name in vOTUs_subset_IDs:
                    write_fna.write(">" + str(name) + "\n" + str(seq) + "\n")

def main():
    print("\n--------------------\r\n")
    print("Running filter_viruses_by_completeness.py\r\n")
    filter_viruses_by_completeness()
    print("\nCompleted filter_viruses_by_completeness.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

