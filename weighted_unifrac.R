##### load libraries #####
library(readr)
library(phyloseq)
library(ape)
library(tidyverse)

# load phyloseq object
load("./create_phyloseq_obj/outputs/fish_phyloseq.Rdata")

# make a distance matrix 
wunifrac_dm <- distance(fish_phyloseq, method="wunifrac")