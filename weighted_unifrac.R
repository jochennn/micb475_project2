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

plot_ordination(fish_phyloseq, pcoa_wu, color = "elevation", shape = "transectname") + 
  scale_color_gradient(low="darkblue", high ="lightblue") + 
  labs(pch = "Transect name", col ="Elevation")

  