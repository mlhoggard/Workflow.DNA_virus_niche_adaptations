#!/usr/bin/env python

'''
compile_dram_annotations.py

Compile annotations.tsv, trnas.tsv, and rrnas.tsv output files from DRAM subset runs (e.g. from manually parallelised run)

MIT License

Copyright (c) 2025 Michael Hoggard

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import os.path
from glob import glob

parser = ArgumentParser()
parser.add_argument('-i', '--input', dest="inpath",
                    help="Path to directory containing output directories for each DRAM(-v) subset", required=True)
parser.add_argument('-o', '--output', dest="outprefix",
                    help="Optional prefix for output filenames. Can include directory path (e.g. 'path/to/output/collated_dram').", default='')
args = parser.parse_args()

def main(): 
    print("\n--------------------\n")
    print("Running compile_dram_annotations.py\n")
    # Loop through the three file types (annotations, rrnas, trnas)
    for file in ['annotations', 'rrnas', 'trnas']:   
        # List of dram subset output directories
        dataframe_paths = sorted(glob(args.inpath+'/*/'+file+'.tsv'))
        # Concatenate subsets (if no objects to concatenated, output empty file)
        try:
            merged_subsets = pd.concat([pd.read_csv(i, sep='\t', index_col=0) for i in dataframe_paths])
        except ValueError:
            merged_subsets = pd.DataFrame()
        # Write out concatenated dataframe
        merged_subsets.to_csv(args.outprefix+file+'.tsv', sep='\t', index=True)
    print("Completed compile_dram_annotations.py\n")
    print("--------------------\n")


if __name__ == '__main__':
    main()

