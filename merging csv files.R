#loading libraries

library(dplyr)
library(readr)

#merging csv files (Colon)

df <- list.files(path="/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Colon", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Merged files")

write.csv(df, "MergedColonFiles.csv")

#merging csv files (Liver)

df <- list.files(path="/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Liver",  full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Merged files")

write.csv(df, "MergedLiverfiles.csv")

#merging csv files (pan-tissue)

df <- list.files(path="/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Pan-tissue", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Merged files")

write.csv(df, "MergedPanTissueFiles.csv")

#merging csv files (Colon treatment only)

df <- list.files(path="/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Colon treatment only", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Merged files")

write.csv(df, "MergedColonTreatmentFiles.csv")

#merging csv files (Colon healthy only)

df <- list.files(path="/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Colon healthy only", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Merged files")

write.csv(df, "MergedColonHealthyFiles.csv")

#merging csv files (Liver treatment only)

df <- list.files(path="/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Liver treatment only", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Merged files")

write.csv(df, "MergedLiverTreatmentFiles.csv")

#merging csv files (Liver healthy only)

df <- list.files(path="/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Liver healthy only", full.names = TRUE) %>%
  lapply(read_csv) %>% 
  bind_rows 

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Merged files")

write.csv(df, "MergedLiverHealthyFiles.csv")

##merging pan tissue files including mutation_type column

df <- list.files(path="/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files 2", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming")

write.csv(df, "MergedPanTissueFilesWithMutType.csv")
write_tsv(df, "mut_type_data.txt")