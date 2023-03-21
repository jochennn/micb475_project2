##### install libraries #####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

##### load libraries #####
library(readr)
library(phyloseq)
library(phyloseq)
library(ape)
library(tidyverse)

##### load in data needed #####
metafp <- "./qiime_export_for_R/binned_fish_metadata.txt"
metadata <- read_delim(metafp, delim="\t")

asvfp <- "C:/Users/17782/Desktop/micb475/micb475_project2/qiime_export_for_R/final_filt-table_EXPORT/feature-table.txt"
asvs_sample_data <- read_delim(file = asvfp, delim="\t", skip=1)

taxfp <- "C:/Users/17782/Desktop/micb475/micb475_project2/qiime_export_for_R/taxonomy_EXPORT/taxonomy.tsv"
taxonomy_data <- read_delim(taxfp, delim="\t")

phylotreefp <- "C:/Users/17782/Desktop/micb475/micb475_project2/qiime_export_for_R/rooted_tree_EXPORT/tree.nwk"
phylotree <- read.tree(phylotreefp)

##### make data phyloseq compatible #####

##### format ASV table #####
# save everything except first column (OTU ID) into a matrix
asv_mat <- as.matrix(asvs_sample_data[,-1])

# Make first column (#OTU ID) the rownames of the new matrix
rownames(asv_mat) <- asvs_sample_data$`#OTU ID`

# Use the "otu_table" function to make an OTU table
ASV_table <- otu_table(asv_mat, taxa_are_rows = TRUE) 
class(ASV_table)


#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(metadata[,-1])
# Make sampleids the rownames
rownames(samp_df)<- metadata$`#SampleID`
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- taxonomy_data %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- taxonomy_data$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
fish_phyloseq <- phyloseq(ASV_table, SAMP, TAX, phylotree)


