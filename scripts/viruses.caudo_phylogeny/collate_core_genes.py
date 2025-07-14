#!/usr/bin/env python

'''
collate_core_genes.py

Collate putative single copy core gene protein sequences from viral genomic datasets.

Method:

- filter datasets to retain only genomes with predicted completeness >= X% (e.g. >= 85%)
- extract gene matches based on vogdb_id core genes identified previously
- for each gene: write out faa file of sequences to align (and then concatenate after alignment)

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
parser.add_argument("-v", "--vogdb_hmmsearch_tblout", dest="vogdb_hmmsearch_tblout",
                    help="tblout results file from hmmsearch seaches of all proteins against the VOG database  (e.g. viralRefSeq.all_genes.vogdb)",
                    metavar="viralRefSeq.all_genes.vogdb", required=True)
parser.add_argument("-p", "--proteins_faa", dest="proteins_faa",
                    help="Protein sequences generated via prodigal or DRAM  (genes.faa)",
                    metavar="genes.faa", required=True)
parser.add_argument("-c", "--checkv_quality_summary", dest="checkv_quality_summary",
                    help="quality_summary.tsv file generated via checkV.",
                    metavar="checkv.quality_summary.tsv", required=True)
parser.add_argument("-t", "--completeness_threshold", dest="completeness_threshold",
                    help="retain genomes >= this completeness percentage threshold (based on CheckV completeness prediction). Default = 85",
                    metavar="completeness_threshold", default=85)
parser.add_argument("-o", "--output_directory", dest="outdir",
                    help="Directory to write output files to. Default = current directory.",
                    metavar='output_directory/', default='./')
args = parser.parse_args()

def collate_core_genes():
    ## marker gene IDs
    marker_genes_df = pd.read_csv(args.core_genes_infile, sep='\t').drop('FunctionalCategory', axis=1)
    marker_genes_dict = dict(marker_genes_df.values)
    ## VOG search results for dataset
    vog_df = pd.read_csv(args.vogdb_hmmsearch_tblout, comment='#', header=None, delimiter=r"\s+", usecols=[*range(0,18)]).reset_index(drop=True)
    vog_df.columns = ['target_name','accession','query_name','accession','fullSeq_eval','fullSeq_score','fullSeq_bias','bestDomain_eval','bestDomain_score','bestDomain_bias','exp','reg','clu','ov','env','dom','rep','inc']     
    vog_df['genomeID'] = vog_df['target_name'].str.replace('_\\d+$', '', regex=True)
    # Only keep match with lowest full seq eval for each geneID
    vog_df = vog_df.sort_values(by=["target_name", "fullSeq_eval"], ascending=[True, True])
    vog_df = vog_df.groupby("target_name", as_index=False).first()
    # filter based on completeness
    checkv_df = pd.read_csv(args.checkv_quality_summary, delimiter="\t")
    checkv_df = checkv_df[checkv_df['completeness'] >= float(args.completeness_threshold)].reset_index(drop=True)
    df = vog_df[vog_df['genomeID'].isin(checkv_df['contig_id'].values)]
    ## Extract marker genes, write sequences to new faa file
    for marker_gene in marker_genes_dict.keys():
        marker_gene_geneIDs = df[df['query_name'] == marker_gene].reset_index(drop=True)['target_name'].values.tolist()
        # generate faa files of amino acid sequences
        with open(args.proteins_faa, 'r') as read_fasta:
            with open(str(args.outdir).rstrip("/") + '/'+marker_gene+'.faa', 'w') as write_faa:
                for name, seq in SimpleFastaParser(read_fasta):
                    gene_id = re.sub(r' .*', r'', name)
                    if gene_id in marker_gene_geneIDs:
                        write_faa.write('>' + str(gene_id) + '_' + marker_gene + '\n' + str(seq) + '\n')


def main():
    print("\n--------------------\r\n")
    print("Running collate_core_genes.py\r\n")
    collate_core_genes()
    print("\nCompleted collate_core_genes.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

