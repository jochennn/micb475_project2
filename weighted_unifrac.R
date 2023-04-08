##### load libraries #####
library(readr)
library(phyloseq)
library(ape)
library(tidyverse)

# load phyloseq object
load("./create_phyloseq_obj/outputs/fish_phyloseq.Rdata")

# make a distance matrix 
wunifrac_dm <- distance(fish_phyloseq, method="wunifrac")
pcoa_wu <- ordinate(fish_phyloseq, method="PCoA", distance = wunifrac_dm)

plot_ordination(fish_phyloseq, pcoa_wu, color = "anatomical_location", shape = "spatial_region") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Weighted Unifrac") 

plot_ordination(fish_phyloseq, pcoa_wu, color = "sample_type", shape = "habitat_depth_level3") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Weighted Unifrac") 

plot_ordination(fish_phyloseq, pcoa_wu, color = "anat_space_combine") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Weighted Unifrac") 

plot_ordination(fish_phyloseq, pcoa_wu, color = "sample_type") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Weighted Unifrac") 




  

  