# Indian Gut Microbiome Project

In this project, we are going to explore the differences in the gut microbiomes of Indian populations living in India and Canada, compared to European controls. 

## Introduction 
The purpose of this study was to investigate the impact of immigration to Canada on the gut microbiome of Indians. Stool samples were collected from the following cohorts:
1. Indians living in India **(Indians)**
2. Indian migrants to Canada **(Indo-Immigr)**
3. Indians born in Canada **(Indo-Can)**
4. Euro-Canadians **(Euro-Can)**
5. Euro-Immigrants **(Euro-Immigr)**

Both 16S and shotgun sequencing were done on all samples to compare bacterial compositions and the functional potential of their microbiomes. Dietary data and demographic data was also collected and analyzed in Prism (analysis not shown here).

## Diagram 
PENDING: chart of all the analysis

## How to Use

For 16S sequence analysis, the following scripts are uploaded:

### 16S-taxonomic-analysis_QIIME2 ###
In this script, paired-end demultiplexed reads were imported into QIIME 2 (Version 2023.9). Quality control was completed with DADA2, which included filtering, chimera removal, dereplication, denoising and merging paired-end reads. For taxonomic classification, the q2-feature-classifer was trained using the GreenGenes2 database (10.28.22). Amplicon sequence variants (ASVs) were filtered for unclassified ASVs and mitochondrial/chloroplast DNA. Using q2-alignment, ASVs were aligned with mafft to construct a phylogeny with fastree via q2-phylogeny. Initial core metric results and taxa bar plots were then generated, then filed were exported for further analysis in R. 

### 16S-taxonomic-analysis_R.R ###
In this script, alpha and beta diversity metrics were plotted using the core libraries phyloseq, ggplot2, vegan and qiime2R. Additionally, LEFSE analysis is included, along with a sample size calculation based on the Bray Curtis distance matrix. 

### 16S-dbRDA-weighted-unifrac.R ### 
In this script, a distance-based redundancy analysis (dbDRA) was done in R using the Weighted UniFrac distance matrix generated from 16S amplicon data in QIIME 2. The purpose of this dbRDA analysis was to understand associations between the distinctions in beta diversity from the taxonomic data with the dietary patterns and baseline characteristics that were significantly different.


For shotgun sequence analysis, the following scripts are uploaded:

### Shotgun-taxonomic-analysis_QIIME2 ###
This script includes the generation of taxa bar plots and LEFSE analysis.

### Shotgun-taxonomic-analysis_R.R ###
In this script, beta diversity analysis was conducted and plotted using the core libraries of phyloseq, ggplot2, and vegan. _Prevotella copri_ clade abundances were also calculated to generate a heatmap. 

### Shotgun-functional-analysis.R ###
In this script, HUMAnN outputs were analyzed. Using LEfSe, differential analysis was conducted to compare abundances of MetaCyc pathways, CAZymes and KEGG modules across cohorts. Barplots were also generated on the stratified CAZy and KEGG data to identify taxonomic contributions for differentially abundant gene families. 


