#load packages
library(tidyverse)

#load metadata file
metadata <- 
  read.delim("./binning_metadata/outputs/binned_fish_metadata.txt", sep = "\t")

#creating modified metadata file
metadata2 <- select(metadata, 
                    "X.SampleID",
                    "habitat_depth_level2", 
                    "sample_type", 
                    "anatomical_location",
                    "spatial_region",
                    "anat_space_combine"
                    )

colnames(metadata2)[1]  <- "#SampleID"

# combining "mesopelagic" and "mesopelagic and benthopelagic" samples together
metadata2$habitat_depth_level3 <- metadata2$habitat_depth_level2

for (i in 1:nrow(metadata2)) {
  if (
    metadata2$habitat_depth_level3[i] == "mesopelagic"|
    metadata2$habitat_depth_level3[i] == "mesopelagic and benthopelagic"
    ) 
    {
    metadata2$habitat_depth_level3[i] <- "meso_benthopelagic"
  }
}

#creating column combining habitat location and anatomical location
metadata2 <- metadata2 %>%
  unite("anat_hab", c("habitat_depth_level3", "sample_type"), remove = F)

##### save data as txt #####
write.table(
  metadata2, 
  file="./binning_metadata/outputs/binned_fish_metadata2.txt", 
  sep="\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  metadata2, 
  file="./create_phyloseq_obj/inputs/binned_fish_metadata2.txt", 
  sep="\t",
  row.names = FALSE,
  quote = FALSE
)




