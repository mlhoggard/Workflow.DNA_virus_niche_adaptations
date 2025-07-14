#!/usr/bin/env python

'''
antiphage_padloc_summary_table.py

Compile summary table of PADLOC results for multiple prokaryote genomes

Note: the script incorporates a 'genomeID' to the summary table based on the input filename (as derived by PADLOC from the input fna filename)

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

parser = ArgumentParser()
parser.add_argument("-p", "--padloc_results_path", dest="padloc_results_path", 
                    help='Path to PADLOC results directory (i.e. location of *_padloc.csv results files)',
                    metavar="path/to/PADLOC_results_directory", required=True)
parser.add_argument("-o", "--output", dest="outfile", 
                    help="(Optional) Output name for PADLOC antiphage defence results summary table. Default = padloc.summary_table.csv",
                    metavar="output_filename.csv", default='padloc.summary_table.csv')
args = parser.parse_args()

def compile_padloc_results():
    infiles = glob.glob(args.padloc_results_path.rstrip('/')+'/*.csv')
    dfs_list = []
    for file in infiles:
        df = pd.read_csv(file)
        # add genome id
        genomeID = re.sub(r'(.*)\..*_padloc.csv', '\\1', re.sub(r'.*/', '', file))
        df['genomeID'] = genomeID
        dfs_list.append(df)
    df = pd.concat(dfs_list, ignore_index=True)
    df = df[['genomeID'] + [col for col in df.columns if col != 'genomeID']]
    df.to_csv(args.outfile, index=False)

def main(): 
    print("\n--------------------\r\n")
    print("Running antiphage_padloc_summary_table.py\r\n")
    compile_padloc_results()
    print("\nCompleted antiphage_padloc_summary_table.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

