##### Packages needed #####
library(tidyr)
library(dplyr)

##### clear workspace ######
rm(list = ls())

##### Read in data #####
fish_metadata <- read.delim("./binning_metadata/inputs/fish_metadata.txt", sep = "\t")

##### Drop rows with missing values and habitat_depth_level2 == "not applicable" #####
fish_mdata_filter<- fish_metadata%>%
  drop_na() %>%
  subset(habitat_depth_level2 != "not applicable")

##### Bin sample_type to internal/external ##### 

# make new column
fish_mdata_filter$anatomical_location <- NA

for (i in 1:nrow(fish_mdata_filter)) {
  if (fish_mdata_filter$sample_type[i] == "hindgut"|fish_mdata_filter$sample_type[i] == "midgut") {
    fish_mdata_filter$anatomical_location[i] <- "internal"
  } else {
    fish_mdata_filter$anatomical_location[i] <- "external"
  }
}

##### Bin habitat_depth_level2 to coastal and open ocean #####

# make new column
fish_mdata_filter$spatial_region <- NA

for (i in 1:nrow(fish_mdata_filter)) {
  if (fish_mdata_filter$habitat_depth_level2[i] == "neritic"|fish_mdata_filter$habitat_depth_level2[i] == "intertidal") {
    fish_mdata_filter$spatial_region[i] <- "coastal"
  } else {
    fish_mdata_filter$spatial_region[i] <- "open ocean"
  }
}


##### Make new column of spatial region & anatomical location combo #####

# make new column
fish_mdata_filter$anat_space_combine<- NA

for (i in 1:nrow(fish_mdata_filter)) {
  if (fish_mdata_filter$spatial_region[i] == "coastal"&fish_mdata_filter$anatomical_location[i] == "internal") {
    fish_mdata_filter$anat_space_combine[i] <- "coastal/internal"
  } else if (fish_mdata_filter$spatial_region[i] == "coastal"&fish_mdata_filter$anatomical_location[i] == "external") {
    fish_mdata_filter$anat_space_combine[i] <- "coastal/external"
  }
    else if (fish_mdata_filter$spatial_region[i] == "open ocean"&fish_mdata_filter$anatomical_location[i] == "internal") {
      fish_mdata_filter$anat_space_combine[i] <- "open ocean/internal"
    }
    else {
    fish_mdata_filter$anat_space_combine[i] <- "open ocean/external"
  }
}

fish_mdata_filter

colnames(fish_mdata_filter)[1]  <- "#SampleID" 
head(fish_mdata_filter)

##### save data as txt #####
write.table(fish_mdata_filter, file="./binning_metadata/outputs/binned_fish_metadata.txt", sep="\t",row.names = FALSE,quote = FALSE)
