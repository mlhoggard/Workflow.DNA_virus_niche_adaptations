#!/usr/bin/env python

'''
dvfpred_extract_fasta.py

Extract sequences for contigs identified as putatively 'viral' by DeepVirFinder

MIT License

Copyright (c) 2025 Michael Hoggard

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from argparse import ArgumentParser
import os
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = ArgumentParser()
parser.add_argument("-d", "--deepvirfinder_results", dest="infile_dvf",
                    help="Results file from DeepVirFinder. e.g. SampleX.fasta_gt1000bp_dvfpred.txt",
                    metavar='...dvfpred.txt', required=True)
parser.add_argument("-a", "--assembly_fasta", dest="infile_fasta",
                    help="Assembly fasta file used as input for DeepVirFinder. e.g. assembly/SampleX.fasta",
                    metavar='assembly_file.fasta', required=True)
parser.add_argument('-o', '--output', dest="outfile_fasta",
                    help="Name for output fasta file of sequences identified by DeepVirFinder", required=True)
args = parser.parse_args()

def main():
    print("\n--------------------\r\n")
    print("Running dvfpred_extract_fasta.py on "+args.infile_dvf+"\r\n")
    # Import DeepVirFinder results
    DVF_df = pd.read_csv(args.infile_dvf, sep='\t')
    print("Extracting sequences from assembly fasta...\r\n")   
    ## Read original assmebly fasta, extract those in in DeepVirFinder results (DVF_df) and write to new fasta file.
    with open(args.infile_fasta, 'r') as read_fasta:
        with open(args.outfile_fasta, 'w') as write_fasta:
            for name, seq in SimpleFastaParser(read_fasta):
                if name in list(set(DVF_df['name'])):
                    write_fasta.write(">" + str(name) + "\n" + str(seq) + "\n")
    print("Output:\r\n")
    print("<sample>_dvfpred.fasta:\nOutput fasta file of putative viral contigs initially identified by dvfind\r\n")
    print("Completed dvfpred_extract_fasta.py on "+args.infile_dvf+"\r\n")
    print("--------------------\r\n")

if __name__ == '__main__':
    main()
