## load packages
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(tidyverse)

## load phyloseq object
load("./create_phyloseq_obj/outputs/fish_phyloseq.Rdata")

## Convert to relative abundance
fish_phyloseq_RA <- 
  transform_sample_counts(fish_phyloseq, fun=function(x) x/sum(x))

metafp <- "./create_phyloseq_obj/inputs/binned_fish_metadata2.txt"
metadata <- read_delim(metafp, delim="\t")

data.frame(metadata$habitat_depth_level3, metadata$sample_type) %>%
  table()

## Filter dataset by "anat_hab"
#meso_benthopelagic
mg_phyloseq <- 
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="meso_benthopelagic_gill")

ms_phyloseq <- 
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="meso_benthopelagic_skin")

mm_phyloseq <- 
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="meso_benthopelagic_midgut")

mh_phyloseq <- 
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="meso_benthopelagic_hindgut")

#neritic
ng_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="neritic_gill")

ns_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="neritic_skin")

nm_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="neritic_midgut")

nh_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="neritic_hindgut")

#intertidal
ig_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="intertidal_gill")

is_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="intertidal_skin")

im_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="intertidal_midgut")

ih_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="intertidal_hindgut")

#abyssalpelagic
ag_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="abyssalpelagic_gill")

as_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="abyssalpelagic_skin")

am_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="abyssalpelagic_midgut")

#bathypelagic***
bg_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="bathypelagic_gill")

bs_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="bathypelagic_skin")

bm_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="bathypelagic_midgut")

bh_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_hab =="bathypelagic_hindgut")

## What ASVs are found in more than 50% of samples in each category?
#meso_benthopelagic
mg <- core_members(mg_phyloseq, detection=0.001, prevalence=0.5)
ms <- core_members(ms_phyloseq, detection=0.001, prevalence=0.5)
mm <- core_members(mm_phyloseq, detection=0.001, prevalence=0.5)
mh <- core_members(mh_phyloseq, detection=0.001, prevalence=0.5)

#neritic
ng <- core_members(ng_phyloseq, detection=0.001, prevalence=0.5)
ns <- core_members(ns_phyloseq, detection=0.001, prevalence=0.5)
nm <- core_members(nm_phyloseq, detection=0.001, prevalence=0.5)
nh <- core_members(nh_phyloseq, detection=0.001, prevalence=0.5)

#intertidal
ig <- core_members(ig_phyloseq, detection=0.001, prevalence=0.5)
is <- core_members(is_phyloseq, detection=0.001, prevalence=0.5)
im <- core_members(im_phyloseq, detection=0.001, prevalence=0.5)
ih <- core_members(ih_phyloseq, detection=0.001, prevalence=0.5)

#abyssalpelagic
ag <- core_members(ag_phyloseq, detection=0.01, prevalence=0.95)
as <- core_members(as_phyloseq, detection=0.01, prevalence=0.95)
am <- core_members(am_phyloseq, detection=0.01, prevalence=0.95)

#bathypelagic
bg <- core_members(bg_phyloseq, detection=0.001, prevalence=0.5)
bs <- core_members(bs_phyloseq, detection=0.001, prevalence=0.5)
bm <- core_members(bm_phyloseq, detection=0.001, prevalence=0.5)
bh <- core_members(bh_phyloseq, detection=0.001, prevalence=0.5)

## venn diagrams
ggVennDiagram(
  x=list(mg, ms, mm, mh),
  category.names=list("mg", "ms", "mm", "mh"),
  label="count"
)

ggVennDiagram(
  x=list(ng, ns, nm, nh),
  category.names=list("ng", "ns", "nm", "nh"),
  label="count"
)

ggVennDiagram(
  x=list(ig, is, im, ih),
  category.names=list("ig", "is", "im", "ih"),
  label="count"
)

ggVennDiagram(
  x=list(ag, as, am),
  category.names=list("ag", "as", "am"),
  label="count"
)

ggVennDiagram(
  x=list(bg, bs, bm, bh),
  category.names=list("bg", "bs", "bm", "bh"),
  label="count"
)

## What are these ASVs?

# found pseudomonas, staphylococcus, streptococcus core microbiome in abyssalpelagic midgut
prune_taxa(am,fish_phyloseq_RA) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_hab, scales="free")

# found lots of synechococcus in the gut and skin
prune_taxa(am,fish_phyloseq_RA) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_hab, scales="free")

prune_taxa(ms,fish_phyloseq) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_hab, scales="free")

prune_taxa(mm,fish_phyloseq) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_hab, scales="free")

prune_taxa(ns,fish_phyloseq) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_hab, scales="free")

prune_taxa(is,fish_phyloseq) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_hab, scales="free")

prune_taxa(ih,fish_phyloseq) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_hab, scales="free")

prune_taxa(bm,fish_phyloseq) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_hab, scales="free")

# visualise synechococcus as common core microbiome genus
ggVennDiagram(
  x=list(am, ms, mm, ns, is, ih, bm),
  category.names=list("A.M.", "M.S.", "M.M.", "N.S.", "I.S.", "I.H.", "B.M."),
  label="count"
)

