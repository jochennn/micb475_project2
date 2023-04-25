##### clear workspace #####
rm(list = ls())

##### packages needed #####
library(stringr)

##### load data #####

fish_manifest <- read.delim("fish_manifest.txt", sep = "\t")

binned_fish_metadata <- read.delim("binned_fish_metadata.txt", sep = "\t")

##### merge and filter for rows in common #####

filtered_fish_manifest <- merge(fish_manifest, binned_fish_metadata, by.x = "sample.id", by.y = "X.SampleID", all = FALSE)

##### select relevant columns ######
filtered_fish_manifest <- filtered_fish_manifest%>%
  select(c("sample.id","absolute.filepath"))

##### rename columns #####
colnames(filtered_fish_manifest)[1]  <- "sample-id" 
colnames(filtered_fish_manifest)[2]  <- "absolute-filepath" 


##### replace $pwd with absolute filepath #####
filtered_fish_manifest$`absolute-filepath` <- str_replace(filtered_fish_manifest$`absolute-filepath`,"[$]PWD","/root/data/project2/fish")

##### save data as txt #####
write.table(filtered_fish_manifest, file="filtered_fish_manifest.txt", sep="\t",row.names = FALSE,quote = FALSE)



