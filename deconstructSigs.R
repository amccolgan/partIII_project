#install and load libraries

install.packages("deconstructSigs")
library(deconstructSigs)
library(dplyr)
library(readr)
library(tidyverse)
library(RColorBrewer)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)

#read vcf files

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/Processed VCF files/Merged files")

input_data <- read.csv("MergedPanTissueFiles.csv")

input_data <- subset(input_data, select = -c(X)) #removing the numerical column added by merging

#formatting data frame for deconstructSigs

input_data <- input_data %>%
  mutate(CHR_col = "chr")  #adds column containing repeat chr string

input_data$CHR <- with(input_data, paste0(CHR_col, CHROM))  #combines columns to reformat CHR

#Convert to deconstructSigs input

sigs.input <- mut.to.sigs.input(mut.ref = input_data,
                                sample.id = "SAMPLE",
                                chr = "CHR",
                                pos = "POS",
                                ref = "REF",
                                alt = "ALT",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)

#counting number of mutations for each sample
total_mutations <- data.frame(rowSums(sigs.input))

#formatting signatures file

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/COSMIC data")
sig_df <- read_tsv("COSMIC_v3.2_SBS_GRCh37.tsv", 
                   col_types = list(Type = col_character(), SBS1 = col_number(), SBS2 = col_number(), SBS17a = col_number(), SBS17b = col_number(), SBS18 = col_number(), SBS31 = col_number(), SBS35 = col_number())) %>%
  select(Type, SBS1, SBS5, SBS17a, SBS17b, SBS18, SBS31, SBS35)
colnames(sig_df)[1] <- NA

cosmic_sigs_temp <- data.frame(t(sig_df[,-1]))
cosmic_names <- t(sig_df[1])

colnames(cosmic_sigs_temp) <- cosmic_names

#rename row names for dataframe recognition by whichSignatures? (didn't work)

#rownames(cosmic_sigs) <- c("Signature.1", "Signature.2", "Signature.3", "Signature.4", "Signature.5", "Signature.6", "Signature.7")

#creating a list of sample id names

sample_names = unique(input_data$SAMPLE)

#checking what signatures are present in samples (1 patient)

test_1 <- whichSignatures(tumor.ref = sigs.input, 
                sample.id = "5FU-PATIENT2-N-CLONE1",
                signatures.ref = cosmic_sigs_temp, 
                associated = c(),
                signatures.limit = NA, 
                signature.cutoff = 0.06, 
                contexts.needed = TRUE,
                tri.counts.method = "default")

#checking signatures for all patients at once

var_names <- list() #to list all the data frames we generate below

for (i in 1:length(sample_names)) {
  temp_variable <- whichSignatures(tumor.ref = sigs.input, 
                                   sample.id = sample_names[i],
                                   signatures.ref = cosmic_sigs_temp, 
                                   associated = c(),
                                   signatures.limit = NA, 
                                   signature.cutoff = 0.06, 
                                   contexts.needed = TRUE,
                                   tri.counts.method = "default")
  var_name <- paste(sample_names[i], "sigs", sep = "_")
  assign(var_name, temp_variable)
  var_names <- append(var_names, var_name)
}

#saving each data frame as a tsv file

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/DeconstructSigs output")

for(i in 1:length(var_names)) {
  var_df <- get(toString(var_names[i]))
  write.table(var_df, paste(var_names[i], ".tsv", sep=""), sep="\t", quote=F, row.names=F)
}

#checking weight differences

sig_diffs <- list()

for(i in 1:length(var_names)) {
  var_df <- get(toString(var_names[i]))
  var_diffs <- var_df$diff
  #name <- var_names[i]
  #assign(name, var_diffs)
  sig_diffs <- append(sig_diffs, var_diffs)
}

max(unlist(sig_diffs))
diffs <- unlist(sig_diffs)
ordered_diffs <- diffs[order(diffs,decreasing = TRUE)]
ordered_diffs_2 <- diffs[order(diffs, decreasing = FALSE)]

######plotting mutational signatures across the entire cohort

# initialize a list of the length of samples 
results <- vector("list", length(sample_names))
names(results) <- sample_names
# run the estimation of exposures for each sample and save the results in the list
for( i in 1:length(sample_names) ){
  sID <- sample_names[i]
  results[[i]] <- whichSignatures(tumor.ref = sigs.input, 
                                  sample.id = sample_names[i],
                                  signatures.ref = cosmic_sigs_temp, 
                                  associated = c(),
                                  signatures.limit = NA, 
                                  signature.cutoff = 0.06, 
                                  contexts.needed = TRUE,
                                  tri.counts.method = "default")
}

#extracting weights data from results

weights_df <- map_dfr(results, "weights")

#format weights_df for signature attribution

library(tibble)
formatted_weights <- tibble::rownames_to_column(weights_df, "sample")
setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming")
write_tsv(formatted_weights, "DeconstructSigsWeightsData.txt")

new_names <- c("Donor1_Clone1_Colon", "Donor1_Clone2_Colon", "Donor10_Clone1_Liver", "Donor10_Clone4_Liver", 
               "Donor11_Clone4_Liver", "Donor12_Clone4_Colon", "Donor13_Clone1_Colon", "Donor13_Clone2_Colon",
               "Donor14_Clone1_Liver", "Donor2_Clone1_Colon", "Donor2_Clone2_Colon", "Donor2_Clone3_Colon", 
               "Donor2_Clone6_Colon", "Donor3_Clone2_Colon", "Donor3_Clone4_Colon", "Donor4_Clone1_Colon",
               "Donor4_Clone2_Colon", "Donor5_Clone2_Colon", "Donor5_Clone3_Colon", "Donor6_Clone2_Liver",
               "Donor6_Clone5_Liver", "Donor7_Clone1_Liver", "Donor7_Clone4_Liver", "Donor8_Clone13_Liver", 
               "Donor8_Clone15_Liver", "Donor9_Clone3_Liver", "c080216D_Clone2_Liver", "c080216D_Clone3_Liver",
               "c110116D_Clone2_Liver", "c150216_Clone2_Liver",  "C30913D_Clone10_Liver", "C30913D_Clone12_Liver",
               "C30913D_Clone13_Liver", "EGAF00000827138_Colon", "EGAF00000827139_Colon", "EGAF00000827140_Colon",
               "EGAF00000827141_Colon", "EGAF00000827249_Colon", "EGAF00000860183_Colon", "EGAF00000860184_Colon",
               "EGAF00000860185_Colon", "EGAF00000860186_Colon",  "EGAF00000860187_Colon", "EGAF00001014392_Colon",
               "EGAF00001014393_Colon", "EGAF00001014394_Colon", "EGAF00001014395_Colon", "ste01166a_Colon",
               "ste01167f_Colon", "ste01205a_Colon", "ste01205f_Colon", "ste1207a_Colon")

new_name_order <- c("EGAF00000827138_Colon", "EGAF00000827139_Colon", "EGAF00000827140_Colon",
                    "EGAF00000827141_Colon", "EGAF00000827249_Colon", "EGAF00000860183_Colon", "EGAF00000860184_Colon",
                    "EGAF00000860185_Colon", "EGAF00000860186_Colon",  "EGAF00000860187_Colon", "EGAF00001014392_Colon",
                    "EGAF00001014393_Colon", "EGAF00001014394_Colon", "EGAF00001014395_Colon", "ste01166a_Colon",
                    "ste01167f_Colon", "ste01205a_Colon", "ste01205f_Colon", "ste1207a_Colon",
                    "Donor1_Clone1_Colon", "Donor1_Clone2_Colon","Donor2_Clone1_Colon", "Donor2_Clone2_Colon", "Donor2_Clone3_Colon", 
                    "Donor2_Clone6_Colon", "Donor3_Clone2_Colon", "Donor3_Clone4_Colon", "Donor4_Clone1_Colon",
                    "Donor4_Clone2_Colon", "Donor5_Clone2_Colon", "Donor5_Clone3_Colon",
                    "Donor12_Clone4_Colon", "Donor13_Clone1_Colon", "Donor13_Clone2_Colon",
                    "c080216D_Clone2_Liver", "c080216D_Clone3_Liver",
                    "c110116D_Clone2_Liver", "c150216_Clone2_Liver",  "C30913D_Clone10_Liver", "C30913D_Clone12_Liver",
                    "C30913D_Clone13_Liver", "Donor6_Clone2_Liver",
                    "Donor6_Clone5_Liver", "Donor7_Clone1_Liver", "Donor7_Clone4_Liver", "Donor8_Clone13_Liver", 
                    "Donor8_Clone15_Liver", "Donor9_Clone3_Liver",
                    "Donor10_Clone1_Liver", "Donor10_Clone4_Liver", "Donor11_Clone4_Liver", "Donor14_Clone1_Liver")

#generating total mutations histogram
rownames(total_mutations) <- new_names
ordered_mutations <- total_mutations[match(new_name_order, row.names(total_mutations)),]
mut_df <- data.frame(new_name_order, ordered_mutations)
rownames(mut_df) <- new_name_order
mut_df$new_name_order <- NULL

colon_1 <- rep(list("untreated colon"), 19)
colon_2 <- rep(list("treated colon"), 15)
liver_1 <- rep(list("untreated liver"), 7)
liver_2 <- rep(list("treated liver"), 11)
sample_type <- c(colon_1, colon_2, liver_1, liver_2)
mut_df$sample_type <- unlist(sample_type)

palette <- RColorBrewer::brewer.pal(length(unique(mut_df$sample_type)),name = 'Set1')
mut_df$color <- palette[as.factor(mut_df$sample_type)]   

par(mar=c(9, 2, 3, 10), mgp = c(6,0,-0.5))
barplot(mut_df$ordered_mutations, main = "Total SBS mutations per sample", xlab = "Sample", names.arg = rownames(mut_df), cex.names = 0.5, las=2,
        col = mut_df$color)

#generating the plot
# convert the exposures for each sample into a sample x signatures matrix
expo <- do.call("rbind", sapply(results, "[", 1))
# add the unknown value to the matrix such that the contributions add up to 1 per sample
Signature.unknown <- unlist(sapply(results, "[", 5))
expo <- cbind(expo, Signature.unknown)
#change names of expo variables using new_names
row.names(expo) <- new_names

#reorder expo rows
new_expo <- expo[match(new_name_order, row.names(expo)),]

o <- row.names(new_expo)
# reorder samples by similarity in their signature profiles
#o <- row.names(expo)[hclust(dist(expo))$order]
# trick base graphics into putting the legend outside of the plot
#coul <- brewer.pal(8, "Set3") 
par(mar=c(9, 2, 3, 9) +0.5)
barplot(t(as.matrix(new_expo))[,o], las=2, col=coul, cex.names = 0.5, main = "Signature Weights",
        ylab = "Signature weights (contribution per sample)")
legend("topright", inset=c(-0.15,0), legend=colnames(new_expo), fill=coul, title="Signature", cex = 0.7, ncol=1, xpd=TRUE)

######

#saving unknown signature weights

setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming/DeconstructSigs output")
expo_IDS <- expo
expo_IDS$sample_ID <- rownames(expo_IDS)
unknown_sigs <- data.frame(expo_IDS$sample_ID, expo_IDS$Signature.unknown)
colnames(unknown_sigs) <- c("Sample ID", "Unknown Signature Weight")
write.csv(unknown_sigs, "unknown_signatures_summary.csv")

#discard samples without SBS31/35 activity

SBS31_35samples <- subset(expo_IDS, SBS31 != 0 | SBS35 != 0)
#if both use & instead, but only gives 4 samples

#plot results of whichSignatures

plotSignatures(test_1, sub = "")

#saving cosmic sigs info for signature attrition
setwd("/Users/ainem/Documents/Cambridge/Part III SysBio/Project/Programming")
cosmic_sigs_temp[1]
rownames(cosmic_sigs_temp) <- c("Type", "SBS1", "SBS5", "SBS17a", "SBS17b", "SBS18", "SBS31", "SBS35")
write_tsv(cosmic_sigs_temp, file="cosmic_sigs_modified2.tsv")

cos_sigs <- read_tsv("cosmic_sigs_modified.tsv")
write_tsv(cos_sigs, "cos_sigs2.txt")