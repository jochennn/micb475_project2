library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggplot2)


## load phyloseq object
load("./create_phyloseq_obj/outputs/fish_phyloseq.Rdata")

fish_phyloseq_CIvOI <- 
  subset_samples(
    fish_phyloseq, 
    anat_space_combine == "open ocean/internal" | 
      anat_space_combine == "coastal/internal"
    )

#adding one to sample counts
fish_phyloseq_plus1 <- 
  transform_sample_counts(fish_phyloseq_CIvOI, function(x) x+1)

#make the deseq object
fish_phyloseq_deseq <- 
  phyloseq_to_deseq2(fish_phyloseq_plus1, ~anat_space_combine)

#run deseq
fish_phyloseq_deseq_output <- DESeq(fish_phyloseq_deseq)

#visualise results
res <- results(fish_phyloseq_deseq_output, tidy=TRUE)

#significang ASVs
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)

# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

#volacano plot
volcano_deseq <- ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

ggsave("./deseq_OIvCI/volcano_deseq.png",
       volcano_deseq
       ,height=4,width=4)

#bar plot
fish_phyloseq_DESeq <- prune_taxa(sigASVs_vec,fish_phyloseq_CIvOI)
sigASVs <- tax_table(fish_phyloseq_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

barplot_deseq <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))

ggsave("./deseq_OIvCI/CIvOI_barplot_deseq.png",
       barplot_deseq
       ,height=5,width=16)
