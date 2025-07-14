#!/usr/bin/env python

'''
virus_niche_summarise_vConTACT2_clusters.py

Summarise vConTACT2 clusters by virus environment type for vOTUs and IMG/VR sequences.

In it's current form this script generates summary tables of environment types associated with clusters for: 1. all clusters that contain a vOTU; and 2. all clusters. For the latter, it also generates Venn diagrams of cluster and environment associations based on: 1. non-saline, saline, terrestrial, and sediment; and 2. non-saline, saline, and terrestrial.

Note: These steps were originally done manually via python code. The code has been put into this script for simplicity, but includes hard-coded portions that are study-specific (e.g. IMG/VR environment types/categorisation of interest, vOTU estuary zone classifications, etc.). For broader use cases this script will need to be redeveloped.

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
import venn
import pylab as plt

parser = ArgumentParser()
parser.add_argument("-v", "--vcontact_results", dest="vcontact_results",
                    help='genome_by_genome_overview.csv generated via vConTACT2',
                    metavar="genome_by_genome_overview.csv", required=True)
parser.add_argument("-i", "--imgvr_sequence_information", dest="imgvr_sequence_information",
                    help='IMG/VR database metadata file IMGVR_all_Sequence_information-high_confidence.tsv',
                    metavar="IMGVR_all_Sequence_information-high_confidence.tsv", required=True)
parser.add_argument("-w", "--waiwera_votus_ecosystems", dest="waiwera_votus_ecosystems",
                    help='IMG/VR database metadata file IMGVR_all_Sequence_information-high_confidence.tsv',
                    metavar="IMGVR_all_Sequence_information-high_confidence.tsv", required=True)
parser.add_argument("-o", "--out_directory", dest="out_directory",
                    help="(Optional) Output directory. Default = current directory",
                    metavar="ouput_directory/", default='./')
args = parser.parse_args()

def read_vcontact():
    # vcontact2 data
    col_names = pd.read_csv(args.vcontact_results, nrows=0).columns
    types_dict = {'VC Size': float, 'Quality': float, 'Adjusted P-value': float, 'VC Avg Distance': float,
                'Topology Confidence Score': float, 'Genus Confidence Score': float, 'VC Kingdoms': float,
                'VC Phyla': float, 'VC Classes': float, 'VC Orders': float, 'VC Families': float,
                'VC Genera': float, 'Genus Confidence Score': float}
    types_dict.update({col: str for col in col_names if col not in types_dict})
    vcontact_df = pd.read_csv(args.vcontact_results, dtype=types_dict)
    # filter for only vOTUs and IMG-VR (exclude refseq)
    vcontact_df = vcontact_df[vcontact_df['Genome'].str.contains('UViG') | vcontact_df['Genome'].str.contains('vOTU')].copy().reset_index(drop=True)
    vcontact_df.columns = vcontact_df.columns.str.replace(' ', '_')
    # Add subcluster column, and strip subcluster string ('_n') from VC column
    vcontact_df['VC_Subcluster'] = vcontact_df['VC']
    vcontact_df['VC'] = vcontact_df['VC'].str.replace(r'(VC_.*)_.*', r'\1', regex=True)
    vcontact_df['UVIG'] = vcontact_df['Genome'].str.extract('(IMGVR[^\|]+)\|')
    return vcontact_df

def read_imgvr():
    # IMG-VR data
    col_names = pd.read_csv(args.imgvr_sequence_information, sep='\t', nrows=0).columns
    types_dict = {'Length': float, 'geNomad score': float, 'Estimated completeness': float, 'Estimated contamination': float}
    types_dict.update({col: str for col in col_names if col not in types_dict})
    img_df = pd.read_csv(args.imgvr_sequence_information, sep='\t', dtype=types_dict)
    return img_df

def read_votu_env():
    votus_env_df = pd.read_csv(args.waiwera_votus_ecosystems, sep='\t').add_prefix('vOTU_')
    return votus_env_df

def add_environment_types(vcontact_df, img_df, votus_env_df):
    # imgvr
    merged_df = pd.merge(
        vcontact_df,
        img_df[['UVIG', 'Ecosystem classification']].copy().reset_index(drop=True),
        how="left", on='UVIG')
    # vOTUs
    merged_df = pd.merge(
        merged_df,
        votus_env_df,
        how="left", left_on='Genome', right_on='vOTU_genomeID').drop(columns=['vOTU_genomeID'])
    return merged_df

def analyse_clusters_votus_only(merged_df):
    clusters_df = merged_df.groupby("VC").filter(lambda x: (x["Genome"].str.contains("UViG")).any() & (x["Genome"].str.contains("vOTU")).any()).reset_index(drop=True)
    # sort by VC and by vOTU>IMGVR>RefSeq
    conds = {}
    conds['IMG_VR'] = clusters_df['Genome'].str.contains('IMGVR', na=False)
    conds['vOTU'] = clusters_df['Genome'].str.contains('vOTU', na=False)
    conds['RefSeq'] = ~clusters_df['Genome'].str.contains('|'.join(['vOTU', 'IMGVR']), na=False)
    clusters_df['Dataset'] = np.select(condlist=conds.values(), choicelist=conds.keys(), default='other')
    conditions = list(map(clusters_df['Dataset'].str.contains, ['vOTU', 'IMG_VR', 'RefSeq']))
    clusters_df['dataset_sort'] = np.select(conditions, [1, 2, 3], 0)
    clusters_df = clusters_df.sort_values(by=['VC', 'dataset_sort', 'Genome'], ascending=[True, True, True]).reset_index(drop=True).drop(columns=['dataset_sort'])
    ## Add 'cluster_environment_summary' (similar to virus taxonomy summaries above) listing all environment types included in each cluster
    # Split Ecosystem classification into separate columns
    clusters_df[['Ecosystem_A','Ecosystem_B','Ecosystem_C','Ecosystem_D']] = clusters_df['Ecosystem classification'].str.split(';',expand=True)
    # create ecosystem_summary column
    conditions = [(clusters_df['vOTU_sample_type'] == 'water_column') & (clusters_df['vOTU_salinity'] == 'fresh'),
                (clusters_df['vOTU_sample_type'] == 'water_column') & (clusters_df['vOTU_salinity'].isin(['brackish', 'marine'])),
                (clusters_df['vOTU_sample_type'] == 'sediment'),
                clusters_df['Ecosystem_C'] == 'Soil',
                clusters_df['Ecosystem_C'] == 'Roots',
                clusters_df['Ecosystem_C'] == 'Grass',
                clusters_df['Ecosystem_C'] == 'Plant litter',
                clusters_df['Ecosystem_C'] == 'Agricultural field',
                clusters_df['Ecosystem_C'] == 'Peat moss',
                clusters_df['Ecosystem_C'] == 'Soil microcosm',
                clusters_df['Ecosystem_C'] == 'Marine',
                clusters_df['Ecosystem_C'] == 'Freshwater',
                clusters_df['Ecosystem_C'] == 'Sediment',
                clusters_df['Ecosystem_D'] == 'Groundwater'
                ]
    choices1 = ['aquatic_freshwater', 'aquatic_saline', 'sediment', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'aquatic_saline', 'aquatic_freshwater', 'sediment', 'aquatic_freshwater']
    clusters_df['ecosystem_summary'] = np.select(conditions, choices1, default='other')
    # filter out 'other' and 'sediment' and 'groundwater'
    clusters_df = clusters_df[clusters_df['ecosystem_summary'] != 'other']
    clusters_df = clusters_df[clusters_df['ecosystem_summary'] != 'sediment']
    # remove any "clusters" of 1 (now that 'other' has been filtered out)
    clusters_df = clusters_df.groupby('VC').filter(lambda x : len(x)>1)
    # Summarise ecosystem for each VC
    clusters_df = clusters_df.join(clusters_df.groupby(by='VC')['ecosystem_summary'].apply(lambda s: list({x for x in s if x != "Unclassified"})), how='right', on='VC', rsuffix='_cluster')
    # add counts dictionary for ecosystem_summary
    clusters_df = clusters_df.join(clusters_df.groupby(by='VC')['ecosystem_summary'].apply(lambda s: str(dict(s.value_counts(dropna=False)))), how='right', on='VC', rsuffix='_cluster_counts')
    # Add proportion for each ecosystem summary
    prop_summary_df = clusters_df.groupby(by='VC')['ecosystem_summary'].apply(lambda s: s.value_counts(normalize=True)).reset_index().pivot_table(index='VC', columns='level_1', values='ecosystem_summary').reset_index().rename(columns={'aquatic_freshwater': 'cluster_prop_freshwater', 'aquatic_saline': 'cluster_prop_saline', 'terrestrial': 'cluster_prop_terrestrial'}).fillna(0)
    clusters_df = clusters_df.merge(prop_summary_df, how='right', on='VC')
    # print VC counts
    print("Total viral clusters containing vOTUs = "+str(len(clusters_df['VC'].unique())))
    print("Viral clusters with more than one associated ecosystem type = "+str(len(clusters_df[clusters_df['ecosystem_summary_cluster_counts'].astype(str).str.contains(',')]['VC'].unique()))+"\r\n")
    # Remove separate ecosystem columns
    keep_columns = [x for x in clusters_df.columns if x not in ['Ecosystem_A','Ecosystem_B','Ecosystem_C','Ecosystem_D']]
    clusters_out_df = clusters_df[keep_columns].copy().reset_index(drop=True)
    # Reorder columns
    clusters_out_df = clusters_out_df[['Dataset', 'Genome',
                                    'VC_Status', 'preVC', 'VC', 'VC_Subcluster', 'VC_Size',
                                    'Ecosystem classification', 'ecosystem_summary', 'ecosystem_summary_cluster', 'ecosystem_summary_cluster_counts',
                                    'cluster_prop_freshwater', 'cluster_prop_saline']]
    # write out
    clusters_out_df.to_csv(args.out_directory.rstrip('/')+'/vcontact2_clusters.ecosystem_summary.vOTUs_only.FreshSaline.tsv', sep='\t', index=False)

def analyse_clusters_votus_and_imgvr_prep(merged_df):
    clusters_df = merged_df.groupby("VC").filter(lambda x: (x["Genome"].str.contains("UViG")).any() | (x["Genome"].str.contains("vOTU")).any()).reset_index(drop=True)
    # sort by VC and by vOTU>IMGVR>RefSeq
    conds = {}
    conds['IMG_VR'] = clusters_df['Genome'].str.contains('IMGVR', na=False)
    conds['vOTU'] = clusters_df['Genome'].str.contains('vOTU', na=False)
    conds['RefSeq'] = ~clusters_df['Genome'].str.contains('|'.join(['vOTU', 'IMGVR']), na=False)
    clusters_df['Dataset'] = np.select(condlist=conds.values(), choicelist=conds.keys(), default='other')
    conditions = list(map(clusters_df['Dataset'].str.contains, ['vOTU', 'IMG_VR', 'RefSeq']))
    clusters_df['dataset_sort'] = np.select(conditions, [1, 2, 3], 0)
    clusters_df = clusters_df.sort_values(by=['VC', 'dataset_sort', 'Genome'], ascending=[True, True, True]).reset_index(drop=True).drop(columns=['dataset_sort'])
    return clusters_df

def analyse_clusters_votus_and_imgvr_four_sets(clusters_df):
    ### Summarise by 4 environment types: non-saline, saline, terrestrial, sediment
    clusters_df[['Ecosystem_A','Ecosystem_B','Ecosystem_C','Ecosystem_D']] = clusters_df['Ecosystem classification'].str.split(';',expand=True)
    ## create ecosystem_summary column
    conditions = [(clusters_df['vOTU_sample_type'] == 'water_column') & (clusters_df['vOTU_salinity'] == 'fresh'),
                (clusters_df['vOTU_sample_type'] == 'water_column') & (clusters_df['vOTU_salinity'].isin(['brackish', 'marine'])),
                (clusters_df['vOTU_sample_type'] == 'sediment'),
                clusters_df['Ecosystem_C'] == 'Soil',
                clusters_df['Ecosystem_C'] == 'Roots',
                clusters_df['Ecosystem_C'] == 'Grass',
                clusters_df['Ecosystem_C'] == 'Plant litter',
                clusters_df['Ecosystem_C'] == 'Agricultural field',
                clusters_df['Ecosystem_C'] == 'Peat moss',
                clusters_df['Ecosystem_C'] == 'Soil microcosm',
                clusters_df['Ecosystem_C'] == 'Marine',
                clusters_df['Ecosystem_C'] == 'Freshwater',
                clusters_df['Ecosystem_C'] == 'Sediment',
                clusters_df['Ecosystem_D'] == 'Groundwater'
                ]
    choices1 = ['freshwater', 'saline', 'sediment', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'saline', 'freshwater', 'sediment', 'groundwater']
    clusters_df['ecosystem_summary'] = np.select(conditions, choices1, default='other')
    # filter out 'other' and and 'groundwater'
    clusters_df = clusters_df[clusters_df['ecosystem_summary'] != 'other']
    clusters_df = clusters_df[clusters_df['ecosystem_summary'] != 'groundwater']
    # remove any "clusters" of 1 (now that 'other' has been filtered out)
    clusters_df = clusters_df.groupby('VC').filter(lambda x : len(x)>1)
    ## Venn diagram
    # make dictionary of sets ({ecosystem: {VCs}}) for Venn
    clusters_venn_dict = {}
    for ecosystem in clusters_df['ecosystem_summary'].unique():
        clusters_venn_dict[ecosystem] = set(clusters_df[clusters_df['ecosystem_summary'] == ecosystem]['VC'])
    # fix freshwater to 'non-saline'
    clusters_venn_dict = { k.replace('freshwater', 'non-saline'): v for k, v in clusters_venn_dict.items() }
    # plot
    col_palette = ["#8ed7b3","#205872",'#960000','#e69e5c']
    vd = venn.venn(clusters_venn_dict, cmap=col_palette, alpha=0.4, fontsize=5, figsize=(4, 4), ax=None)
    plt.savefig(args.out_directory.rstrip('/')+'/vcontact2_clusters.ecosystem_summary.vOTUs_and_IMGVR.FreshSalineTerrestrialSediment.venn.pdf',format="pdf")

def analyse_clusters_votus_and_imgvr_three_sets(clusters_df):
    ### Summarise by 3 environment types: non-saline, saline, terrestrial
    clusters_df[['Ecosystem_A','Ecosystem_B','Ecosystem_C','Ecosystem_D']] = clusters_df['Ecosystem classification'].str.split(';',expand=True)
    # create ecosystem_summary column
    conditions = [(clusters_df['vOTU_sample_type'] == 'water_column') & (clusters_df['vOTU_salinity'] == 'fresh'),
                (clusters_df['vOTU_sample_type'] == 'water_column') & (clusters_df['vOTU_salinity'].isin(['brackish', 'marine'])),
                (clusters_df['vOTU_sample_type'] == 'sediment'),
                clusters_df['Ecosystem_C'] == 'Soil',
                clusters_df['Ecosystem_C'] == 'Roots',
                clusters_df['Ecosystem_C'] == 'Grass',
                clusters_df['Ecosystem_C'] == 'Plant litter',
                clusters_df['Ecosystem_C'] == 'Agricultural field',
                clusters_df['Ecosystem_C'] == 'Peat moss',
                clusters_df['Ecosystem_C'] == 'Soil microcosm',
                clusters_df['Ecosystem_C'] == 'Marine',
                clusters_df['Ecosystem_C'] == 'Freshwater',
                clusters_df['Ecosystem_C'] == 'Sediment',
                clusters_df['Ecosystem_D'] == 'Groundwater'
                ]
    choices1 = ['freshwater', 'saline', 'sediment', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'terrestrial', 'saline', 'freshwater', 'sediment', 'groundwater']
    clusters_df['ecosystem_summary'] = np.select(conditions, choices1, default='other')
    # filter out 'other' and 'sediment' and 'groundwater'
    clusters_df = clusters_df[clusters_df['ecosystem_summary'] != 'other']
    clusters_df = clusters_df[clusters_df['ecosystem_summary'] != 'sediment']
    clusters_df = clusters_df[clusters_df['ecosystem_summary'] != 'groundwater']
    # remove any "clusters" of 1 (now that 'other' has been filtered out)
    clusters_df = clusters_df.groupby('VC').filter(lambda x : len(x)>1)
    ## Venn diagram
    # make dictionary of sets ({ecosystem: {VCs}}) for venn
    clusters_venn_dict = {}
    for ecosystem in clusters_df['ecosystem_summary'].unique():
        clusters_venn_dict[ecosystem] = set(clusters_df[clusters_df['ecosystem_summary'] == ecosystem]['VC'])
    # fix freshwater to 'non-saline'
    clusters_venn_dict = { k.replace('freshwater', 'non-saline'): v for k, v in clusters_venn_dict.items() }
    # plot
    col_palette = ["#8ed7b3","#205872",'#960000']
    vd = venn.venn(clusters_venn_dict, cmap=col_palette, alpha=0.4, fontsize=5, figsize=(4, 4), ax=None)
    plt.savefig(args.out_directory.rstrip('/')+'/vcontact2_clusters.ecosystem_summary.vOTUs_and_IMGVR.FreshSalineTerrestrial.venn.pdf',format="pdf")
    ## Environment type summaries per cluster
    # Summarise ecosystem for each VC
    clusters_df = clusters_df.join(clusters_df.groupby(by='VC')['ecosystem_summary'].apply(lambda s: list({x for x in s if x != "Unclassified"})), how='right', on='VC', rsuffix='_cluster')
    # add counts dictionary for ecosystem_summary
    clusters_df = clusters_df.join(clusters_df.groupby(by='VC')['ecosystem_summary'].apply(lambda s: str(dict(s.value_counts(dropna=False)))), how='right', on='VC', rsuffix='_cluster_counts')
    # Add proportions for ecosystem types
    prop_summary_df = clusters_df.groupby(by='VC')['ecosystem_summary'].apply(lambda s: s.value_counts(normalize=True)).reset_index().pivot_table(index='VC', columns='level_1', values='ecosystem_summary').reset_index().rename(columns={'freshwater': 'cluster_prop_freshwater', 'saline': 'cluster_prop_saline', 'terrestrial': 'cluster_prop_terrestrial'}).fillna(0)
    prop_summary_df.to_csv(args.out_directory.rstrip('/')+'/vcontact2_clusters.ecosystem_summary.vOTUs_and_IMGVR.FreshSalineTerrestrial.VC_summary.tsv', sep='\t', index=False)
    clusters_df = clusters_df.merge(prop_summary_df, how='right', on='VC')
    # print VC counts
    print("Total viral clusters containing vOTUs or IMG/VR sequences = "+str(len(clusters_df['VC'].unique())))
    print("Viral clusters with more than one associated ecosystem type = "+str(len(clusters_df[clusters_df['ecosystem_summary_cluster_counts'].astype(str).str.contains(',')]['VC'].unique())))
    print("Percentage of viral clusters with only one associated ecosystem type = "+str(round(((1-((len(clusters_df[clusters_df['ecosystem_summary_cluster_counts'].astype(str).str.contains(',')]['VC'].unique())) / (len(clusters_df['VC'].unique()))))*100), 2))+"%\r\n")
    ## Proportions of each ecosystem type
    results_tmp = []
    for ecosystem in ['freshwater', 'saline', 'terrestrial']:
        results_tmp.append(clusters_df.groupby('VC').first()[['cluster_prop_'+ecosystem]].value_counts().rename_axis('unique_values').reset_index(name='cluster_prop_'+ecosystem).set_index('unique_values'))
    counts_table = pd.concat(results_tmp, axis=1, ignore_index=False).reset_index()
    # proportion where 100% one ecosystem type
    counts_table_pct = counts_table[counts_table['unique_values'] == 1].copy().reset_index(drop = True)
    counts_table_pct[[col for col in counts_table_pct if 'cluster_prop' in col]] /= len(clusters_df['VC'].unique())
    counts_table_pct = pd.melt(counts_table_pct, id_vars=["unique_values"], var_name="ecosystem_type").drop(columns='unique_values')
    print("Summary table of propotions of clusters that contain exclusively one environment type (freshwater, saline, terrestrial)")
    print(counts_table_pct)
    print("Percentage of clusters that contain exclusively one environment type (freshwater, saline, terrestrial) = "+str(round((100*counts_table_pct['value'].sum()), 2))+"\r\n")
    # proportion where >=90% one ecosystem type
    counts_table_pct = counts_table[counts_table['unique_values'] >= 0.9].copy().reset_index(drop = True).drop(columns='unique_values').sum(numeric_only=True, axis=0) / len(clusters_df['VC'].unique())
    print("Summary table of propotions of clusters that contain >=90% one environment type (freshwater, saline, terrestrial)")
    print(counts_table_pct)
    print("Percentage of clusters that contain >=90% environment type (freshwater, saline, terrestrial) = "+str(round((100*counts_table_pct.sum()), 2))+"\r\n")
    ## write out summary table
    # Remove separate ecosystem columns
    keep_columns = [x for x in clusters_df.columns if x not in ['Ecosystem_A','Ecosystem_B','Ecosystem_C','Ecosystem_D']]
    clusters_out_df = clusters_df[keep_columns].copy().reset_index(drop=True)
    # Reorder columns
    clusters_out_df = clusters_out_df[['Dataset', 'Genome',
                                    'VC_Status', 'preVC', 'VC', 'VC_Subcluster', 'VC_Size',
                                    'Ecosystem classification', 'ecosystem_summary', 'ecosystem_summary_cluster', 'ecosystem_summary_cluster_counts',
                                    'cluster_prop_freshwater', 'cluster_prop_saline']]
    # write out
    clusters_out_df.to_csv(args.out_directory.rstrip('/')+'/vcontact2_clusters.ecosystem_summary.vOTUs_and_IMGVR.FreshSalineTerrestrial.tsv', sep='\t', index=False)

def main():
    print("\n--------------------\r\n")
    print("Running virus_niche_summarise_vConTACT2_clusters.py\r\n")
    print("Reading in files\r\n")
    vcontact_df = read_vcontact()
    img_df = read_imgvr()
    votus_env_df = read_votu_env()
    merged_df = add_environment_types(vcontact_df, img_df, votus_env_df)
    print("Analysing clusters that contain at least one Waiwera vOTU\r\n")
    analyse_clusters_votus_only(merged_df)
    print("Analysing all clusters that contain Waiwera vOTUs or IMG/VR sequences\r\n")
    clusters_df = analyse_clusters_votus_and_imgvr_prep(merged_df)
    print("Processing clusters based on four environment types: freshwater, saline, terrestrial, and sediment\r\n")
    analyse_clusters_votus_and_imgvr_four_sets(clusters_df)
    print("Processing clusters based on three environment types: freshwater, saline, and terrestrial\r\n")
    analyse_clusters_votus_and_imgvr_three_sets(clusters_df)
    print("\nCompleted virus_niche_summarise_vConTACT2_clusters.py\n")
    print("--------------------\n")

if __name__ == '__main__':
    main()

