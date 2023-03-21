##### load libraries #####
library(readr)
library(phyloseq)
library(ape)
library(tidyverse)

##### load in data needed #####
metafp <- "./create_phyloseq_obj/inputs/binned_fish_metadata.txt"
metadata <- read_delim(metafp, delim="\t")

asvfp <- "./create_phyloseq_obj/inputs/feature-table.txt"
asvs_sample_data <- read_delim(file = asvfp, delim="\t", skip=1)

taxfp <- "./create_phyloseq_obj/inputs/taxonomy.tsv"
taxonomy_data <- read_delim(taxfp, delim="\t")

phylotreefp <- "./create_phyloseq_obj/inputs/tree.nwk"
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


