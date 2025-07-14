# Data processing overview for the study: 

#### DNA virus adaptations are influenced by environment and proliferation is constrained by environmental niche  

Michael Hoggard, Emilie Gios, Hwee Sze Tee, Jemma L. Geoghegan, Kim M. Handley

----

## Background

To investigate estuarine viral diversity, niche constraints, and genomic traits of environmental adaptation, we analysed metagenomic and metatranscriptomic data from across an estuarine salinity gradient, including water and sediment habitats. We then expanded our analysis to globally distributed viral genomes from the IMG/VR database. 

These docs present an overview of the data processing workflow used in this study.

Software (and versions), custom scripts, and additional databases used are listed below.

## Software used throughout this workflow

- Trimmomatic v0.38 / v0.39
- SortMeRNA v4.3.6 (database: smr_v4.3_fast_db)
- metaSPAdes v3.11.1
- seqmagick v0.7.0
- BBMap v37.93 / v39.01
- CONCOCT v0.4.1
- MetaBAT v2.12.1
- MaxBin v2.2.4
- DAS_Tool v1.1.1
- CheckM v1.2.1
- GTDB-TK v2.4.0 (database v214)
- dRep v1.4.3
- DRAM v1.3.5 / v1.4.6
- CRISPRDetect v2.4
- PADLOC v2.0.0 (database: PADLOC-DB v2.0.0)
- seqmagick v0.7.0
- VirSorter2 v2.2.3
- VIBRANT v1.2.1
- DeepVirFinder v1.0
- Kraken2 v2.0.9
- extract_kraken_reads.py (from KrakenTools)
- seqmagick v0.7.0
- CheckV v0.7.0
- Cluster_genomes_5.1.pl
- vConTACT2 v2.0.11.3
- Cytoscape v3.8.2
- prodigal-gv v2.9.0
- HMMER v3.3.2
- Clustal-Omega v1.2.4
- IQ-TREE v2.2.2.2
- crass v1.0.1
- blast v2.13.0
- pepstats: EMBOSS v6.6.0
- SAMtools v1.19
- featureCounts (subread-2.0.6-Linux-x86_64)

## Custom scripts used throughout this workflow (available in scripts/)

- viruses.identification/dvfind_add_fdr.R
- viruses.identification/dvfpred_extract_fasta.py
- viruses.identification/dvfpred_filter_euk.py
- viruses.identification/summarise_viral_contigs.py
- viruses.identification/virome_per_sample_derep.py
- viruses.identification/checkv_filter_contigs.py
- viruses.identification/filter_viruses_by_completeness.py
- viruses.identification/ictv_reconcile_refseq_taxonomy.py
- viruses.identification/vcontact2_update_refseq_taxonomy.py 
- viruses.identification/tax_predict_vConTACT2.0.11.3.tax_update_202309.py
- viruses.caudo_phylogeny/filter_refseq_by_taxonomy.py
- viruses.caudo_phylogeny/identify_core_genes.py
- viruses.caudo_phylogeny/collate_core_genes.py
- viruses.caudo_phylogeny/concatenate_protein_alignments.py
- viruses.general/virus_niche_summarise_vConTACT2_clusters.py
- viruses.general/hostmatch_crispr_summary_table.py
- general/compile_dram_annotations.py
- general/dramv_compile_summary_table.py
- general/antiphage_padloc_summary_table.py
- general/extract_aa_and_gc_proportions.py
- general/summarise_pepstats_pI.py
- general/featurecounts_make_feature_table.py
- general/summarise_counts.py
- general/summarise_counts.R

Note: The following MIT licence applies to all custom scripts

> MIT License
> 
> Copyright (c) 2025 Michael Hoggard
> 
> Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
> 
> The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
> 
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Additional databases used in this workflow

- IMG/VR v7.1
- viralRefSeq v211 (within vConTACT2)
- viralRefSeq v223 (for taxonomy updates)
- ICTV taxonomy MSL38.v3 (propagation of higher ranks for taxonomy updates to viralRefSeq)
- VOGdb v222

----