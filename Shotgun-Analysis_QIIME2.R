################################################################################
################################################################################
###############################   QIIME 2   ####################################
################################################################################
################################################################################

##### SESSION INFO #####
# QIIME 2 Version: 2023.9 
# R Version: 4.2.2

###

## First export files in R into versions for QIIME 2
library(phyloseq)

# Assuming your phyloseq object names are otu_table and taxonomy
otu_matrix <- as.matrix(otu_table(ps))
taxonomy_table <- as.matrix(tax_table(ps))

write.table(otu_matrix, "otu_table.tsv", sep="\t")

write.table(taxonomy_table, "taxonomy.tsv", sep="\t")


conda activate /opt/homebrew/Caskroom/miniconda/base/envs/qiime2-amplicon-2023.9

###

# Convert otu table to biom
biom convert -i otu_table.tsv -o otu-table.biom --to-hdf5

# Import the otu table into QIIME 2
qiime tools import \
--input-path otu-table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path otu-table.qza

# Import the taxonomy table into QIIME 2
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format TSVTaxonomyFormat \
--input-path taxonomy.tsv \
--output-path taxonomy.qza

qiime taxa barplot \
--i-table otu-table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file sample-metadata_V5.txt \
--o-visualization taxa-bar-plots-meta.qzv

# Fix to just for genus filter
qiime taxa filter-table \
--i-table otu-table.qza \
--i-taxonomy taxonomy.qza \
--p-include g__ \
--o-filtered-table genus-only-table.qza

qiime taxa barplot \
--i-table genus-only-table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file sample-metadata_V5.txt \
--o-visualization taxa-bar-plots-meta-genus.qzv
