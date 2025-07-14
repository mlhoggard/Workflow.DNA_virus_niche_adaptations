#!/usr/bin/env python

'''
identify_core_genes.py

Assess Caudoviricetes viralRefSeq VOGdb gene hits for putative single copy core genes.

Method:

- References filtered for >= x% completeness (predicted via CheckV) (default = 95% completeness threshold)
- Marker genes selected based on the following criteria:
  (1) present in >= 10% of referece virus genomes
  (2) average gene copy number <= 1.2
  (3) average predicted protein length > 100 amino acid residues.

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
import os
from glob import glob
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = ArgumentParser()
parser.add_argument("-v", "--vogdb_hmmsearch_tblout", dest="vogdb_hmmsearch_tblout",
                    help="tblout results file from hmmsearch seaches of all proteins against the VOG database  (e.g. viralRefSeq.all_genes.vogdb)",
                    metavar="viralRefSeq.all_genes.vogdb", required=True)
parser.add_argument("-a", "--vogdb_annotations_file", dest="vodgb_annotations",
                    help="VOG database annotations file (vog.annotations.tsv)",
                    metavar="vog.annotations.tsv", required=True)
parser.add_argument("-p", "--proteins_faa", dest="proteins_faa",
                    help="Protein sequences generated via prodigal or DRAM  (genes.faa)",
                    metavar="genes.faa", required=True)
parser.add_argument("-c", "--checkv_quality_summary", dest="checkv_quality_summary",
                    help="quality_summary.tsv file generated via checkV.",
                    required=True)
parser.add_argument("-t", "--completeness_threshold", dest="completeness_threshold",
                    help="retain genomes >= this completeness percentage threshold (based on CheckV completeness prediction). Default = 95",
                    metavar="completeness_threshold", default=95)
parser.add_argument("-o", "--output_directory", dest="outdir",
                    help="Directory to write output files to. Default = current directory.",
                    metavar='output_directory/', default='./')
args = parser.parse_args()

def identify_core_genes():
    df = pd.read_csv(args.vogdb_hmmsearch_tblout, comment='#', header=None, delimiter=r"\s+", usecols=[*range(0,18)]).reset_index(drop=True)
    df.columns = ['target_name','accession','query_name','accession','fullSeq_eval','fullSeq_score','fullSeq_bias','bestDomain_eval','bestDomain_score','bestDomain_bias','exp','reg','clu','ov','env','dom','rep','inc']
    df['genomeID'] = df['target_name'].str.replace('_\\d+$', '', regex=True)
    # Only keep match with lowest full seq eval for each geneID
    df = df.sort_values(by=["target_name", "fullSeq_eval"], ascending=[True, True])
    df = df.groupby("target_name", as_index=False).first()
    # filter for high-quality and/or complete genomes
    checkv_df = pd.read_csv(args.checkv_quality_summary, delimiter="\t")
    checkv_df = checkv_df[checkv_df['completeness'] >= float(args.completeness_threshold)].reset_index(drop=True)
    df = df[df['genomeID'].isin(checkv_df['contig_id'].values)]
    # Write out vog annotations summary
    df.to_csv(str(args.outdir).rstrip("/") + '/vogdb_hmmsearch.compiled_results.tsv', sep='\t', index=False)
    # VOGdb hits per genome
    vog_counts_df = df.groupby('genomeID')["query_name"].apply(lambda x: x.groupby(x).size()).reset_index().pivot(index='genomeID', columns='level_1', values='query_name').reset_index()
    # Select vog hits based on desired criteria
    ## present in more than 10% of genomes
    vog_counts_df = vog_counts_df.dropna(thresh=int(len(vog_counts_df)*0.1), axis=1)
    ## average copy number ≤1.2
    cols_mean = vog_counts_df.mean(axis=0, numeric_only=True)
    drop_cols = cols_mean[cols_mean >= 1.2].index
    vog_counts_df.drop(columns=drop_cols, inplace=True)
    vogIDs_to_keep = [col for col in vog_counts_df.columns if col != 'genomeID']
    ## average protein length >100 amino acid residues
    seq_dict = {}
    with open(args.proteins_faa, 'r') as read_fasta:
        for name, seq in SimpleFastaParser(read_fasta):
            name_trim = re.sub(r' .*', r'', name)
            seq_dict[name_trim] = len(seq)
    genes_df = df.merge(pd.DataFrame(seq_dict.items(), columns=['geneID','aa_seq_length']), how='left', left_on='target_name', right_on='geneID')
    genes_df = genes_df[genes_df['query_name'].isin(vogIDs_to_keep)]
    mean_seq_len_df = genes_df.groupby('query_name', as_index=False)['aa_seq_length'].mean()
    # final filtered list of vogIDs
    final_vogIDs = mean_seq_len_df[mean_seq_len_df['aa_seq_length'] >= 100]['query_name'].values.tolist()
    # counts table of vogs per genome
    vog_counts_final_df = vog_counts_df[['genomeID']+[col for col in vog_counts_df.columns if col in final_vogIDs]]
    # write out
    vog_counts_final_df.to_csv(str(args.outdir).rstrip("/") + '/core_genes.countTable.tsv', sep='\t', index=False)
    # compile table of final vogIDs w/annotations
    vog_annot = pd.read_csv(args.vodgb_annotations, delimiter="\t")
    vog_annot = vog_annot[vog_annot['#GroupName'].isin(final_vogIDs)].reset_index(drop=True).rename(columns={'#GroupName': 'vogdb_id'})
    vog_annot = vog_annot[['vogdb_id', 'FunctionalCategory', 'ConsensusFunctionalDescription']]
    # write out
    vog_annot.to_csv(str(args.outdir).rstrip("/") + '/core_genes.annotations.tsv', sep='\t', index=False)


def main():
    print("\n--------------------\r\n")
    print("Running identify_core_genes.py\r\n")
    identify_core_genes()
    print("\nCompleted identify_core_genes.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

