# install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2",force = TRUE)
BiocManager::install("Maaslin2")
BiocManager::install("lefser")

install.packages("ggpicrust2")
devtools::install_github("cafferychen777/ggpicrust2")


# load libraries
library(data.table)
library(phyloseq)
library(tidyverse)
library(ggpicrust2)
library(MicrobiomeStat)
library(dplyr)

# libraries for bar plot function
library(patchwork)
library(ggprism)


##### metadata loading and filtering #####
# read in metadata 
fish_binned_metadata <- read.delim("./binning_metadata/outputs/binned_fish_metadata.txt", sep = "\t")

# filter for internal only metadata for internal only comparison
internal_mdata <- filter(fish_binned_metadata, anatomical_location == "internal")

##### KO abundance load and filtering #####
# Convert KO abundance in picrust2 export files to KEGG pathway abundance - TAKES LONG TO RUN
kegg_abundance <- ko2kegg_abundance(file="./KO_metagenome_out/pred_metagenome_unstrat.tsv")

# filter metadata so that samples are matching (for external v internal comp)
kept_samples_ext_int <- colnames(kegg_abundance)
matching_metadata_ext_int <- subset(fish_binned_metadata, X.SampleID %in% kept_samples_ext_int)

# filter metadata so that samples are matching (for internal, spatial region comp)
kept_samples_int_only <- colnames(kegg_abundance)
matching_metadata_int_only <- subset(internal_mdata, X.SampleID %in% kept_samples_int_only)

# filter abundances so that samples match with metadata (for external v internal comp)
kegg_abundance_match_ext_int <- dplyr::select(kegg_abundance, matches(matching_metadata_ext_int$X.SampleID))

# filter abundances so that samples match with metadata (for internal, spatial region comp)
kegg_abundance_match_int_only <- dplyr::select(kegg_abundance, matches(matching_metadata_int_only$X.SampleID))

##### EXTERNAL V INTERNAL PATHWAY COMP BAR PLOT WORKFLOW #####

# differential abundance analysis using limma voom 
anat_daa_results_df <- pathway_daa(abundance = kegg_abundance_match_ext_int,
                                   metadata = matching_metadata_ext_int,
                                   daa_method = "limma voom",
                                   group = "anatomical_location" ,
                                   select = NULL,
                                   reference = NULL)



# filter to 30 lowest adj p-values
anat_daa_results_sig_df <- anat_daa_results_df[order(anat_daa_results_df$p_adjust),]
anat_daa_results_sig_df <- slice(anat_daa_results_sig_df,1:30)

# save table of pathways with 30 lowest p-adjusted values
write.csv(anat_daa_results_sig_df, file = "ext_int_kegg_pathway_table.csv", row.names = FALSE)


# pathway annotation
anat_daa_results_sig_df <- pathway_annotation(pathway = "KO", daa_results_df = anat_daa_results_sig_df, ko_to_kegg = TRUE)
# pathway annotation lead to different number of pathways for the different categories?? Not too sure why
# need to ensure total number of pathways annotated doesn't exceed 30


# plot anatomical location daa
anat_ebar_plot <- pathway_errorbar(abundance = kegg_abundance_match_ext_int,
                                   daa_results_df = anat_daa_results_sig_df,
                                   Group = matching_metadata_ext_int$anatomical_location,
                                   ko_to_kegg = TRUE,
                                   p_values_threshold = 0.05,
                                   order = "pathway_class",
                                   p_value_bar = TRUE,
                                   colors = NULL,
                                   x_lab = NULL)
print(anat_ebar_plot)

# save plot
ggsave(filename = "ext_int_mscript_plot.png"
       , anat_ebar_plot
       , height=10, width=17)






##### INTERNAL DIFFERENTIAL ANALYSIS BY SPATIAL REGION BAR PLOT WORKFLOW #####

# differential abundance analysis using limma voom 
int_daa_results_df <- pathway_daa(abundance = kegg_abundance_match_int_only,
                              metadata = matching_metadata_int_only,
                              daa_method = "limma voom",
                              group = "spatial_region" ,
                              select = NULL,
                              reference = NULL)



# filter to 30 lowest adj p-values
int_daa_results_sig_df <- int_daa_results_df[order(int_daa_results_df$p_adjust),]
int_daa_results_sig_df <- slice(int_daa_results_sig_df,1:30)

# save table of pathways with 30 lowest p-adjusted values
write.csv(int_daa_results_sig_df, file = "int_only_kegg_pathway_table.csv", row.names = FALSE)

# pathway annotation
int_daa_results_sig_df <- pathway_annotation(pathway = "KO", daa_results_df = int_daa_results_sig_df, ko_to_kegg = TRUE)
# pathway annotation lead to different number of pathways for the different categories?? Not too sure why
# need to ensure total number of pathways annotated doesn't exceed 30


# plot anatomical location daa
int_ebar_plot <- pathway_errorbar(abundance = kegg_abundance_match_int_only,
                  daa_results_df = int_daa_results_sig_df,
                  Group = matching_metadata_int_only$spatial_region,
                  ko_to_kegg = TRUE,
                  p_values_threshold = 0.05,
                  order = "pathway_class",
                  p_value_bar = TRUE,
                  colors = NULL,
                  x_lab = NULL) 
print(int_ebar_plot)

# save plot
ggsave(filename = "int_spatial_mscript_plot.png"
       , int_ebar_plot
       , height=7, width=18.5)



