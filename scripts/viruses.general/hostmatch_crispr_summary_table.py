#!/usr/bin/env python

'''
hostmatch_crispr_summary_table.py

Compile summary table of crispr spacer blast results for virus genomes and putative prokaryote hosts.

Note: blast results must be in the following format: "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

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
parser.add_argument("-p", "--prokaryote_host_crispr_blast_results", dest="host_crispr_blast_results",
                    help='CRISPR spacer blastn results for (prokaryote) hosts (Must be in the following format: "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs")',
                    metavar="blastn_crisprSpacers.hosts.txt", required=True)
parser.add_argument("-v", "--viruses_crispr_blast_results", dest="viruses_crispr_blast_results",
                    help='CRISPR spacer blastn results for viruses (Must be in the following format: "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs")',
                    metavar="blastn_crisprSpacers.viruses.txt", required=True)
parser.add_argument("-b", "--blast_hits_threshold", dest="blast_hits_threshold",
                    help="Number of top blast hits to keep. Default = 100",
                    metavar='100', default=100)
parser.add_argument("-m", "--mismatch_threshold", dest="mismatch_threshold",
                    help="Threshold of allowable mismatches from blast search. Default = 0",
                    metavar='0', default=0)
parser.add_argument("-t", "--min_virhost_matches_threshold", dest="min_virhost_matches_threshold",
                    help="(Optional) Minimum number of distinct spacer matches between individual viruses and hosts to retain results for each virus-host pair  (e.g. for -t 2, virus-host pairs are excluded if they share <2 unique spacer sequences between them)",
                    metavar='1', default=None)
parser.add_argument("-o", "--output", dest="outfile",
                    help="(Optional) Output name for crispr spacer host matching summary table. Default = virhost.crispr_spacers.summary_table.tsv",
                    metavar='', default='virhost.crispr_spacers.summary_table.tsv')
args = parser.parse_args()

def host_results():
    df = pd.read_csv(args.host_crispr_blast_results, sep='\t', header=None)
    df.columns = ['host_contig_ID', 'query_length', 'spacer_ID', 'spacer_len', 'pident', 'match_length', 'mismatch', 'gapopen', 'query_start', 'query_end', 'spacer_start', 'spacer_end', 'evalue', 'bitscore', 'query_covs']
    # Filter to only keep matches with: <= n mismatch over the *full* length of the spacer sequence, & zero gaps
    df = df[(df['spacer_len'] == df['match_length']) & (df['mismatch'] == int(args.mismatch_threshold)) & (df['gapopen'] == 0)]
    return df

def viruses_results():
    df = pd.read_csv(args.viruses_crispr_blast_results, sep='\t', header=None)
    df.columns = ['virus_ID', 'query_length', 'spacer_ID', 'spacer_len', 'pident', 'match_length', 'mismatch', 'gapopen', 'query_start', 'query_end', 'spacer_start', 'spacer_end', 'evalue', 'bitscore', 'query_covs']
    # Filter to only keep matches with: <= n mismatch over the *full* length of the spacer sequence, & zero gaps
    df = df[(df['spacer_len'] == df['match_length']) & (df['mismatch'] == int(args.mismatch_threshold)) & (df['gapopen'] == 0)]
    return df

def merge_results(virus_df, host_df):
    # merge virus and host results (by spacer_ID) and filter to keep only rows that have hits for both a viral contig and a binned contig
    df = pd.merge(virus_df, host_df, how="left", on="spacer_ID", suffixes=("_virus", "_host")).reset_index(drop=True)
    df = df[df['host_contig_ID'].notnull()].sort_values(by=['virus_ID']).reset_index(drop=True)
    # Filter to keep columns of interest
    df = df[['virus_ID', 'host_contig_ID', 'spacer_ID', 'pident_virus', 'evalue_virus', 'bitscore_virus', 'pident_host', 'evalue_host', 'bitscore_host']]
    df.columns = ['virus_ID', 'host_contig_ID', 'crispr_blast_spacer_ID', 'crispr_blast_pident_virus', 'crispr_blast_evalue_virus', 'crispr_blast_bitscore_virus', 'crispr_blast_pident_host', 'crispr_blast_evalue_host', 'crispr_blast_bitscore_host']
    return(df)

def filter_by_minimum_matches(df):
    # Note, filters based on minimum number of virus-host matches per host *contig* rather than per genome (i.e. assumes the whole CRISPR spacer-repeat region will be contained on a single contig).
    df['VirHost_paired_ID'] = df['virus_ID'] + '__' + df['host_contig_ID']
    # identify virus-mag matches with >=n distinct spacer match
    df_counts_tmp = pd.DataFrame(df[['virus_ID','host_contig_ID', 'crispr_blast_spacer_ID']].sort_values(by=['virus_ID', 'host_contig_ID']).drop_duplicates().reset_index(drop=True)[['virus_ID','host_contig_ID']].groupby(['virus_ID']).value_counts()).reset_index()
    df_counts_tmp.columns = ['virus_ID', 'host_contig_ID', 'counts']
    df_counts_tmp = df_counts_tmp[df_counts_tmp['counts'] >= int(args.min_virhost_matches_threshold)].reset_index(drop=True)
    df_counts_tmp['VirHost_paired_ID'] = df_counts_tmp['virus_ID'] + '__' + df_counts_tmp['host_contig_ID']
    df = df[df['VirHost_paired_ID'].isin(df_counts_tmp['VirHost_paired_ID'])]
    df = df.drop(columns=['VirHost_paired_ID'])
    return(df)

def filter_results(df):
    # ERROR handling: If n_hits_threshold greater than or equal to max counts, need to modify n_hits_threshold
    MAX_VALUE_COUNTS = df.groupby('virus_ID')['virus_ID'].value_counts().max()
    if int(args.blast_hits_threshold) >= MAX_VALUE_COUNTS:
        n_hits_threshold = MAX_VALUE_COUNTS-1
    else:
        n_hits_threshold = int(args.blast_hits_threshold)
    df = df[df.index.isin(df.groupby('virus_ID')['crispr_blast_evalue_virus'].nsmallest(n_hits_threshold).index.get_level_values(1))].sort_values(by=['virus_ID', 'crispr_blast_evalue_virus'])
    ## pivot wider and apply suffix to multiple hits
    df['idx'] = '_'+(df.groupby(['virus_ID']).cumcount() + 1).astype(str)
    df = (
        df.pivot_table(
            index=['virus_ID'],
            columns=['idx'],
            values=['host_contig_ID', 'crispr_blast_spacer_ID', 'crispr_blast_pident_virus', 'crispr_blast_evalue_virus', 'crispr_blast_bitscore_virus', 'crispr_blast_pident_host', 'crispr_blast_evalue_host', 'crispr_blast_bitscore_host'],
            aggfunc='first'
        )
    )
    df.columns = [''.join(col) for col in df.columns]
    df = df.reset_index()
    return df

def main():
    print("\n--------------------\r\n")
    print("Running hostmatch_crispr_summary_table.py\r\n")
    host_df = host_results()
    virus_df = viruses_results()
    summary_df = merge_results(virus_df, host_df)
    if args.min_virhost_matches_threshold is not None:
        summary_df = filter_by_minimum_matches(summary_df)
    filt_summary_df = filter_results(summary_df)
    filt_summary_df.to_csv(args.outfile, sep='\t', index=False)
    print("\nCompleted hostmatch_crispr_summary_table.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()
