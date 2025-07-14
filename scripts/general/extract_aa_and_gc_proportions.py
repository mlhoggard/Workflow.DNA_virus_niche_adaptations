#!/usr/bin/env python

'''
extract_aa_and_gc_proportions.py

Calculate proportions of amino acid sequences in predicted proteins and extract gene GC content from files generated via DRAM, prodigal, or prodigal-gv.

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

parser = ArgumentParser()
parser.add_argument("-f", "--input_format", dest="format",
                    help="Input files format. Options = ['DRAM', 'prodigal']. If -f DRAM, --protein_sequences, --dram_genbank, and --dram_annotations inputs required. If -f prodigal, only --protein_sequences input required.",
                    metavar="[DRAM|prodigal]", required=True)
parser.add_argument("-p", "--protein_sequences", dest="genes_faa",
                    help="Predicted protein sequences generated via DRAM or prodigal (genes.faa)",
                    metavar='genes.faa', required=True)
parser.add_argument("-g", "--dram_genbank", dest="genbank",
                    help="GenBank formatted file generated via DRAM (genes.gbk)",
                    metavar='genes.gbk', default=None)
parser.add_argument("-a", "--dram_annotations", dest="annotations",
                    help="Gene annotations file generated via DRAM (annotations.tsv)",
                    metavar='annotations.tsv', default=None)
parser.add_argument("-s", "--sample_id", dest="sample_id",
                    help="ID string (optional). Adds sample ID column to summary table (useful for compiling results from multiple samples downstream)",
                   metavar='sampleID', default=None)
parser.add_argument("-o", "--output_filename", dest="out_filename",
                    help="Filename for results summary table (tsv format). Default = AA_and_GC.summary_table.tsv.",
                    metavar='output_filename.tsv', default='AA_and_GC.summary_table.tsv')
args = parser.parse_args()

def from_DRAM():
    # Establish empty variable lists to append results
    contigID_list = []
    geneID_list = []
    gene_count_list = []
    seq_aa_list = []
    aa_seq_length_list = []
    glu_list = []
    asp_list = []
    arg_list = []
    lys_list = []
    hist_list = []
    glycine_list = []
    serine_list = []
    tyrosine_list = []
    cysteine_list = []
    glutamine_list = []
    asparagine_list = []
    threonine_list = []
    phenylalanine_list = []
    leucine_list = []
    tryptophan_list = []
    proline_list = []
    isoleucine_list = []
    methionine_list = []
    valine_list = []
    alanine_list = []
    aa_acidic_list = []
    aa_basic_list = []
    aa_polar_list = []
    aa_nonpolar_list = []
    aa_charged_list = []
    # Loop through each sample genes.faa file to extract values for each variable and append to lists
    # Also calculate propotions for amino  acids and summaries of aa types.
    with open(args.genes_faa, 'r') as read_fasta:
        for name, seq in SimpleFastaParser(read_fasta):
            headers = name.split(' ')
            contigID_list.append(re.sub('_\d+$', '', headers[0].strip()))
            geneID_list.append(headers[0].strip())
            gene_count_list.append(re.sub(r'.*_(\d+)$', r'\1', headers[0].strip()))
            seq_aa_list.append(seq)
            aa_seq_length = len(seq)
            aa_seq_length_list.append(aa_seq_length)
            glu_list.append(seq.count('E')/aa_seq_length)
            asp_list.append(seq.count('D')/aa_seq_length)
            arg_list.append(seq.count('R')/aa_seq_length)
            lys_list.append(seq.count('K')/aa_seq_length)
            hist_list.append(seq.count('H')/aa_seq_length)
            glycine_list.append(seq.count('G')/aa_seq_length)
            serine_list.append(seq.count('S')/aa_seq_length)
            tyrosine_list.append(seq.count('Y')/aa_seq_length)
            cysteine_list.append(seq.count('C')/aa_seq_length)
            glutamine_list.append(seq.count('Q')/aa_seq_length)
            asparagine_list.append(seq.count('N')/aa_seq_length)
            threonine_list.append(seq.count('T')/aa_seq_length)
            phenylalanine_list.append(seq.count('F')/aa_seq_length)
            leucine_list.append(seq.count('L')/aa_seq_length)
            tryptophan_list.append(seq.count('W')/aa_seq_length)
            proline_list.append(seq.count('P')/aa_seq_length)
            isoleucine_list.append(seq.count('I')/aa_seq_length)
            methionine_list.append(seq.count('M')/aa_seq_length)
            valine_list.append(seq.count('V')/aa_seq_length)
            alanine_list.append(seq.count('A')/aa_seq_length)
            aa_acidic_list.append((seq.count('E')/aa_seq_length)+(seq.count('D')/aa_seq_length))
            aa_basic_list.append((seq.count('R')/aa_seq_length)+(seq.count('K')/aa_seq_length)+(seq.count('H')/aa_seq_length))
            aa_polar_list.append((seq.count('D')/aa_seq_length)+(seq.count('E')/aa_seq_length)+(seq.count('H')/aa_seq_length)+(seq.count('K')/aa_seq_length)+(seq.count('N')/aa_seq_length)+(seq.count('Q')/aa_seq_length)+(seq.count('R')/aa_seq_length)+(seq.count('S')/aa_seq_length)+(seq.count('T')/aa_seq_length)+(seq.count('Z')/aa_seq_length))
            aa_nonpolar_list.append((seq.count('A')/aa_seq_length)+(seq.count('C')/aa_seq_length)+(seq.count('F')/aa_seq_length)+(seq.count('G')/aa_seq_length)+(seq.count('I')/aa_seq_length)+(seq.count('L')/aa_seq_length)+(seq.count('M')/aa_seq_length)+(seq.count('P')/aa_seq_length)+(seq.count('V')/aa_seq_length)+(seq.count('W')/aa_seq_length)+(seq.count('Y')/aa_seq_length))
            aa_charged_list.append((seq.count('B')/aa_seq_length)+(seq.count('D')/aa_seq_length)+(seq.count('E')/aa_seq_length)+(seq.count('H')/aa_seq_length)+(seq.count('K')/aa_seq_length)+(seq.count('R')/aa_seq_length)+(seq.count('Z')/aa_seq_length))
    # collate df
    df = pd.DataFrame({'contigID': contigID_list,
                    'geneID': geneID_list,
                    'gene_count': gene_count_list,
                    'aa_seq': seq_aa_list,
                    'aa_seq_length': aa_seq_length_list,
                    'aa_prop_Glutamic_acid': glu_list,
                    'aa_prop_Aspartic_acid': asp_list,
                    'aa_prop_Arginine': arg_list,
                    'aa_prop_Lysine': lys_list,
                    'aa_prop_Histidine': hist_list,
                    'aa_prop_Glycine': glycine_list,
                    'aa_prop_Serine': serine_list,
                    'aa_prop_Tyrosine': tyrosine_list,
                    'aa_prop_Cysteine': cysteine_list,
                    'aa_prop_Glutamine': glutamine_list,
                    'aa_prop_Asparagine': asparagine_list,
                    'aa_prop_Threonine': threonine_list,
                    'aa_prop_Phenylalanine': phenylalanine_list,
                    'aa_prop_Leucine': leucine_list,
                    'aa_prop_Tryptophan': tryptophan_list,
                    'aa_prop_Proline': proline_list,
                    'aa_prop_Isoleucine': isoleucine_list,
                    'aa_prop_Methionine': methionine_list,
                    'aa_prop_Valine': valine_list,
                    'aa_prop_Alanine': alanine_list,
                    'aa_prop_Acidic': aa_acidic_list,
                    'aa_prop_Basic': aa_basic_list,
                    'aa_prop_Polar': aa_polar_list,
                    'aa_prop_Nonpolar': aa_nonpolar_list,
                    'aa_prop_charged': aa_charged_list
                    })
    # Pull in annotations and gbk files to make df of other fields of interest (e.g. start, end, seq length, gc content)
    gbk_seq_ids = []
    gbk_gene_ids = []
    gbk_gene_gc = []
    with open(args.genbank, 'r') as genbank_infile:
        for seq_record in SeqIO.parse(genbank_infile, "genbank"):
            for gene in seq_record.features:
                if gene.type == 'CDS':
                    gbk_seq_ids.append(seq_record.id)
                    gbk_gene_ids.append(gene.qualifiers['gene'][0])
                    gbk_gene_gc.append(gene.qualifiers['gc_cont'][0])
    annot_df = pd.read_csv(args.annotations, sep='\t').rename(columns={'Unnamed: 0':'geneID', 'start_position': 'start','end_position': 'end'})
    annot_df = annot_df[['geneID', 'start', 'end']+[col for col in annot_df if 'vogdb' in col]]
    annot_df = pd.merge(
        pd.DataFrame({'contigID': gbk_seq_ids, 'geneID': gbk_gene_ids, 'gene_gc': gbk_gene_gc}),
        annot_df,
        how = 'outer',
        on = 'geneID')
    # collate with AA results
    full_df = pd.merge(df, annot_df, how='left', on=['contigID', 'geneID'])
    full_df = full_df.rename(columns={'start': 'nt_start', 'end': 'nt_end', 'gene_gc': 'gc_content'})
    full_df['nt_seq_length'] = full_df['nt_end'] - full_df['nt_start']
    # write out summary table
    if args.sample_id:
        full_df['sampleID'] = args.sample_id
    full_df.to_csv(args.out_filename, sep='\t', index=False)

def from_prodigal():
    # Establish empty variable lists to append results
    contigID_list = []
    geneID_list = []
    gene_count_list = []
    start_list = []
    end_list = []
    nt_seq_length_list = []
    gc_content_list = []
    seq_aa_list = []
    aa_seq_length_list = []
    glu_list = []
    asp_list = []
    arg_list = []
    lys_list = []
    hist_list = []
    glycine_list = []
    serine_list = []
    tyrosine_list = []
    cysteine_list = []
    glutamine_list = []
    asparagine_list = []
    threonine_list = []
    phenylalanine_list = []
    leucine_list = []
    tryptophan_list = []
    proline_list = []
    isoleucine_list = []
    methionine_list = []
    valine_list = []
    alanine_list = []
    aa_acidic_list = []
    aa_basic_list = []
    aa_polar_list = []
    aa_nonpolar_list = []
    aa_charged_list = []
    # prodigal.faa file to extract values for each variable and append to lists
    # Also calculate propotions for amino  acids
    with open(args.genes_faa, 'r') as read_fasta:
        for name, seq in SimpleFastaParser(read_fasta):
            headers = name.split('#')
            contigID_list.append(re.sub('_\d+$', '', headers[0].strip()))
            geneID_list.append(headers[0].strip())
            gene_count_list.append(re.sub(r'.*_(\d+)$', r'\1', headers[0].strip()))
            start_list.append(float(headers[1].strip()))
            end_list.append(float(headers[2].strip()))
            nt_seq_length_list.append(float(headers[2].strip())-float(headers[1].strip()))
            gc_content_list.append(float(headers[-1].split('gc_cont=')[-1].strip()))
            seq_aa_list.append(seq)
            aa_seq_length = len(seq)
            aa_seq_length_list.append(aa_seq_length)
            glu_list.append(seq.count('E')/aa_seq_length)
            asp_list.append(seq.count('D')/aa_seq_length)
            arg_list.append(seq.count('R')/aa_seq_length)
            lys_list.append(seq.count('K')/aa_seq_length)
            hist_list.append(seq.count('H')/aa_seq_length)
            glycine_list.append(seq.count('G')/aa_seq_length)
            serine_list.append(seq.count('S')/aa_seq_length)
            tyrosine_list.append(seq.count('Y')/aa_seq_length)
            cysteine_list.append(seq.count('C')/aa_seq_length)
            glutamine_list.append(seq.count('Q')/aa_seq_length)
            asparagine_list.append(seq.count('N')/aa_seq_length)
            threonine_list.append(seq.count('T')/aa_seq_length)
            phenylalanine_list.append(seq.count('F')/aa_seq_length)
            leucine_list.append(seq.count('L')/aa_seq_length)
            tryptophan_list.append(seq.count('W')/aa_seq_length)
            proline_list.append(seq.count('P')/aa_seq_length)
            isoleucine_list.append(seq.count('I')/aa_seq_length)
            methionine_list.append(seq.count('M')/aa_seq_length)
            valine_list.append(seq.count('V')/aa_seq_length)
            alanine_list.append(seq.count('A')/aa_seq_length)
            aa_acidic_list.append((seq.count('E')/aa_seq_length)+(seq.count('D')/aa_seq_length))
            aa_basic_list.append((seq.count('R')/aa_seq_length)+(seq.count('K')/aa_seq_length)+(seq.count('H')/aa_seq_length))
            aa_polar_list.append((seq.count('D')/aa_seq_length)+(seq.count('E')/aa_seq_length)+(seq.count('H')/aa_seq_length)+(seq.count('K')/aa_seq_length)+(seq.count('N')/aa_seq_length)+(seq.count('Q')/aa_seq_length)+(seq.count('R')/aa_seq_length)+(seq.count('S')/aa_seq_length)+(seq.count('T')/aa_seq_length)+(seq.count('Z')/aa_seq_length))
            aa_nonpolar_list.append((seq.count('A')/aa_seq_length)+(seq.count('C')/aa_seq_length)+(seq.count('F')/aa_seq_length)+(seq.count('G')/aa_seq_length)+(seq.count('I')/aa_seq_length)+(seq.count('L')/aa_seq_length)+(seq.count('M')/aa_seq_length)+(seq.count('P')/aa_seq_length)+(seq.count('V')/aa_seq_length)+(seq.count('W')/aa_seq_length)+(seq.count('Y')/aa_seq_length))
            aa_charged_list.append((seq.count('B')/aa_seq_length)+(seq.count('D')/aa_seq_length)+(seq.count('E')/aa_seq_length)+(seq.count('H')/aa_seq_length)+(seq.count('K')/aa_seq_length)+(seq.count('R')/aa_seq_length)+(seq.count('Z')/aa_seq_length))
    # collate df
    df = pd.DataFrame({'contigID': contigID_list,
                    'geneID': geneID_list,
                    'aa_prop_Lysine': lys_list,
                    'aa_prop_Histidine': hist_list,
                    'aa_prop_Glycine': glycine_list,
                    'aa_prop_Serine': serine_list,
                    'aa_prop_Tyrosine': tyrosine_list,
                    'aa_prop_Cysteine': cysteine_list,
                    'aa_prop_Glutamine': glutamine_list,
                    'aa_prop_Asparagine': asparagine_list,
                    'aa_prop_Threonine': threonine_list,
                    'aa_prop_Phenylalanine': phenylalanine_list,
                    'aa_prop_Leucine': leucine_list,
                    'aa_prop_Tryptophan': tryptophan_list,
                    'aa_prop_Proline': proline_list,
                    'aa_prop_Isoleucine': isoleucine_list,
                    'aa_prop_Methionine': methionine_list,
                    'aa_prop_Valine': valine_list,
                    'aa_prop_Alanine': alanine_list,
                    'aa_prop_Acidic': aa_acidic_list,
                    'gene_count': gene_count_list,
                    'nt_start': start_list,
                    'nt_end': end_list,
                    'nt_seq_length': nt_seq_length_list,
                    'gc_content': gc_content_list,
                    'aa_seq': seq_aa_list,
                    'aa_seq_length': aa_seq_length_list,
                    'aa_prop_Glutamic_acid': glu_list,
                    'aa_prop_Aspartic_acid': asp_list,
                    'aa_prop_Arginine': arg_list,
                    'aa_prop_Basic': aa_basic_list,
                    'aa_prop_Polar': aa_polar_list,
                    'aa_prop_Nonpolar': aa_nonpolar_list,
                    'aa_prop_charged': aa_charged_list
                    })
    # write out summary table
    if args.sample_id:
        df['sampleID'] = args.sample_id
        df_filt = df[['sampleID', 'contigID','geneID','gc_content']+[col for col in df.columns if 'aa_prop' in col]]
    else:
        df_filt = df[['contigID','geneID','gc_content']+[col for col in df.columns if 'aa_prop' in col]]
    df_filt.to_csv(args.out_filename, sep='\t', index=False)


def main():
    print("\n--------------------\r\n")
    print("Running extract_aa_and_gc_proportions.py\r\n")
    if args.format.lower() == 'dram'.lower():
        from_DRAM()
    elif args.format.lower() == 'prodigal'.lower():
        from_prodigal()
    else:
        print('Error: Invalid input file format (-f) specificed. Must be one of: ["DRAM", "prodigal"]. Exiting.')
        exit()
    print("\nCompleted extract_aa_and_gc_proportions.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

