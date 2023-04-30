#install and load packages

install.packages("vcfR")
install.packages("tidyverse")

library(vcfR)
library(tidyverse)

#reading vcf files

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Datasets/Kuijk et al/Extracted data")
vcf <- read.vcfR("5FU-PATIENT1-N-CLONE1.purple.somatic.vcf", verbose = FALSE )

# extract all the INFO and FORMAT fields into a list of tidy
# data frames: fix, gt, and meta. Here we don't coerce columns
# to integer or numeric types...
Z <- vcfR2tidy(vcf)
names(Z)

# here is the meta data in a table
Z$meta
# here is the fixed info
Z$fix

# here are the GT fields.  Note that ChromKey and POS are keys
# back to Z$fix
Z$gt

# here we put the data into a single, joined data frame (list component
# dat in the returned list) and the meta data.  Let's just pick out a 
# few fields:
df <- vcfR2tidy(vcf, single_frame = TRUE)
colnames(df$dat)

head(df)

#working out data formatting for signature attribution
test <- df$dat[c("CHROM", "POS", "REF", "ALT", "Indiv", "SGT")]
colnames(test)[5] <- "SAMPLE"
colnames(test)[6] <- "MUTATION TYPE"
test

#creating a loop that will process each vcf file into a dataframe at once

temp = list.files(pattern="*.vcf")

for (i in 1:length(temp)) {
  vcf_temp <- read.vcfR(temp[i], verbose = FALSE)
  sample_id <- str_replace(temp[i], ".purple.somatic.vcf", "")
  X <- vcfR2tidy(vcf_temp, single_frame = TRUE, verbose = FALSE)
  df_temp <- X$dat
  saved_df <- df_temp[c("CHROM", "POS", "REF", "ALT", "Indiv")]
  colnames(saved_df)[5] <- "SAMPLE"
  var_name <- paste("df", sample_id, sep = "_")
  assign(var_name, saved_df)
}

#creating a list of data frame names

for (i in 1:length(temp)) {
  sample_name <- str_replace(temp[i], ".purple.somatic.vcf", "")
  temp[i] <- paste("df", sample_name, sep = "_")
}

#now we are doing the same thing, but including mutation type data for signature attribution

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Datasets/Kuijk et al/Extracted data")

temp = list.files(pattern="*.vcf")

for (i in 1:length(temp)) {
  vcf_temp <- read.vcfR(temp[i], verbose = FALSE)
  sample_id <- str_replace(temp[i], ".purple.somatic.vcf", "")
  X <- vcfR2tidy(vcf_temp, single_frame = TRUE, verbose = FALSE)
  df_temp <- X$dat
  saved_df <- df_temp[c("CHROM", "POS", "REF", "ALT", "Indiv", "SGT")]
  colnames(saved_df) <- c("chr", "pos", "ref", "alt", "sample", "mutation_type")
  var_name <- paste("df2", sample_id, sep = "_")
  assign(var_name, saved_df)
}

#creating a list of data frame names

for (i in 1:length(temp)) {
  sample_name <- str_replace(temp[i], ".purple.somatic.vcf", "")
  temp[i] <- paste("df2", sample_name, sep = "_")
}

######

#saving data frames as csv files

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files")

for (i in 1:length(temp)) {
  temp_df = get(temp[i])
  write.csv(temp_df, paste(temp[i],".csv", sep=""), row.names = FALSE)
}

#creating one big csv file

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/One big VCF summary")

data_frame_list <- list()
for (i in 1:length(temp)) {
  data_frame_list[i] <- get(temp[i])
}

###saving csv files but for the second loop (signature assignment)

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files 2")

for (i in 1:length(temp)) {
  temp_df = get(temp[i])
  write.csv(temp_df, paste(temp[i],".csv", sep=""), row.names = FALSE)
}

###

big_df <- rbind(data_frame_list)
head(big_df)


