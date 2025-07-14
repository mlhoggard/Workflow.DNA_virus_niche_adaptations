#!/usr/bin/env python

'''
filter_refseq_by_taxonomy.py

Filter viralRefSeq genomic.gbff based on taxonomy string and write out genomic.fna file.

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
import glob
import os
from Bio import SeqIO

parser = ArgumentParser()
parser.add_argument("-i", "--input_file", dest="infile",
                    help="viralRefSeq genomic GenBank file (viral.genomic.gbff)",
                    metavar="viral.genomic.gbff", required=True)
parser.add_argument("-t", "--taxonomy_string", dest="taxonomy_string",
                    help="Taxonomy search string to subset sequences based on GenBank seq_record.annotations['taxonomy'] entries. E.g. 'Caudoviricetes'.",
                    required=True)
parser.add_argument("-o", "--output_directory", dest="outdir",
                    help="Directory to write output files to (tax_subset.genomic.fna; tax_subset.genomic.metadata.tsv). Default = current directory.",
                    metavar='output_directory/', default='./')
args = parser.parse_args()

def subset_by_taxonomy_string():
    ### Genomic files
    metadata_dfs = []
    with open(args.infile, 'r') as genbank_infile:
        with open(str(args.outdir).rstrip("/") + '/' + str(args.taxonomy_string) + '.genomic.fna', 'w') as write_fasta:
            for seq_record in SeqIO.parse(genbank_infile, "genbank") :
                if str(args.taxonomy_string) in seq_record.annotations['taxonomy']:
                    write_fasta.write(">%s\n%s\n" % (
                        seq_record.id,
                        str(seq_record.seq)))
                    try:
                        host = seq_record.features[0].qualifiers['host'][0]
                    except Exception as e:
                        host = np.nan
                    metadata_dfs.append(pd.DataFrame({'genomeID': seq_record.id, 'full_taxonomy': ';'.join(seq_record.annotations['taxonomy']), 'host': host}, index=[0]))
    # compile metadata df
    df = pd.concat(metadata_dfs).reset_index(drop=True)
    # separate taxonomy into ranks based on suffix extraction
    df['domain'] = np.where(df['full_taxonomy'].astype(str).str.contains('Viruses.*'), df['full_taxonomy'].astype(str).str.replace(pat='(Viruses).*', repl='\\1', regex=True), 'Unclassified')
    tax_rank = ['realm', 'subrealm', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'suborder', 'family', 'subfamily', 'genus']
    suffix_str = ['viria', 'vira', 'virae', 'virites', 'viricota', 'viricotina', 'viricetes', 'viricetidae', 'virales', 'virineae', 'viridae', 'virinae', 'virus']
    for i in range(len(tax_rank)):
        df[tax_rank[i]] = np.where(df['full_taxonomy'].astype(str).str.contains('.*;.*'+suffix_str[i]+'.*'), df['full_taxonomy'].astype(str).str.replace(pat='.*;(.*'+suffix_str[i]+').*', repl='\\1', regex=True), 'Unclassified')
    df['species'] = np.where(df['full_taxonomy'].astype(str).str.contains('.*;.*virus .*'), df['full_taxonomy'].astype(str).str.replace(pat='.*;(.*virus .*)', repl='\\1', regex=True), 'Unclassified')
    # write out metadata file
    df.to_csv(args.outdir.rstrip("/") + '/' + args.taxonomy_string + '.genomic.metadata.tsv', sep='\t', index=False)


def main():
    print("\n--------------------\r\n")
    print("Running filter_refseq_by_taxonomy.py\r\n")
    subset_by_taxonomy_string()
    print("\nCompleted filter_refseq_by_taxonomy.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

