#Note following code can be done for both COI/18S datasets

# Load the installed packages
suppressPackageStartupMessages({
  library(ggpubr)
  library(phyloseq)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
})

# Display package versions
cat("Package Versions:\n")
cat("ggplot2:", packageVersion("ggplot2"), "\n")
cat("phyloseq:", packageVersion("phyloseq"), "\n")

# Define the path to the directories
script_directory <- "C:/Users/miche/OneDrive - Cawthron/Google_Drive/PhD Project/Chapter 3_dispersal of eDNA_eRNA/Data Analyze/Metabarcoding_data_2/Pest Alert tool"
input_rds_directory <- "C:/Users/miche/OneDrive - Cawthron/Google_Drive/PhD Project/Chapter 3_dispersal of eDNA_eRNA/Data Analyze/Metabarcoding_data_2/input_ps"
output_rds_directory <- "C:/Users/miche/OneDrive - Cawthron/Google_Drive/PhD Project/Chapter 3_dispersal of eDNA_eRNA/Data Analyze/Metabarcoding_data_2/input_ps"
fig_output_directory <- "C:/Users/miche/OneDrive - Cawthron/Google_Drive/PhD Project/Chapter 3_dispersal of eDNA_eRNA/Data Analyze/Metabarcoding_data_2/Figures"

# Load phyloseq object data processed by DADA2, with contamination from controls removed by subtracting the maximum sequence count across controls from ASVs
Opua18S.ps <- readRDS(file.path(input_rds_directory, "Opua18S.subtract.ps.rds"))
OpuaCOI.ps <- readRDS(file.path(input_rds_directory, "OpuaCOI.subtract.ps.rds"))

# Remove singletons from the phyloseq objects
Opua18S.ps <- phyloseq::filter_taxa(Opua18S.ps, function(x) sum(x) > 0, TRUE) 
Opua18S.ps <- prune_samples(sample_sums(Opua18S.ps) > 0, Opua18S.ps)

# Note: The following code can be done for both COI and 18S datasets.

# Instructions for the PAT tool
# Go to PAT tool via this link: https://pestalert.builtbyupshift.dev/
# Choose your 18S/COI files which have been processed after DADA2 (i.e., after primer removal, merging, chemical removal, etc.)
# Download the filtered NIS as a text file and format it like this:
# bb7b18bae54454c87557df8217b5422b, Polydora_cornuta
# 06881212cb8d13d0b6ea902a9b331b28, Amphibalanus_amphitrite
# 43efaa6e1de0e23ac2ad8cba7ff086a9, Oratosquilla_oratoria

# Subset the marine NIS from the COI/18S phyloseq objects
imp_ASVs <- read.csv("18S_PAT.csv", header = FALSE)
colnames(imp_ASVs) <- c("ASVs", "PATSpecies")

newTaxa <- imp_ASVs$ASVs
allTaxa <- taxa_names(Opua18S.ps)
allTaxa <- allTaxa[!(allTaxa %in% newTaxa)]

phyloseq_subset <- prune_taxa(allTaxa, Opua18S.ps)
phyloseq_NIS_18S <- prune_taxa(colnames(otu_table(Opua18S.ps)) %in% imp_ASVs$ASVs, Opua18S.ps)

# Summarize the phyloseq object
summarize_phyloseq(phyloseq_NIS_18S)
ntaxa(phyloseq_NIS_18S)

# Export taxonomic data
NISr_taxa <- as(tax_table(phyloseq_NIS_18S), "matrix")
write.csv(NISr_taxa, "NIS_taxa_18S.csv")

# Read taxonomic data and check row names
tax_NISr <- read.csv("NIS_taxa_18S1.csv", header = TRUE, row.names = 1) %>% as.matrix()
stopifnot(all(row.names(tax_NISr) == openssl::md5(row.names(tax_NISr))))  # Check MD5 hashes
row.names(tax_NISr) <- row.names(tax_NISr)

NISr <- as(otu_table(phyloseq_NIS_18S), "matrix")
NISrdf <- as.data.frame(NISr)

# Extract the sample data from 18S phyloseq object
sample_data_opua18S <- sample_data(Opua18S.ps)

# Create final phyloseq object for marine NIS
phyloseq_NIS_18S <- phyloseq(otu_table(NISrdf, taxa_are_rows = FALSE),
                              sample_data_opua18S, 
                              tax_table(tax_NISr))

# Save the final filtered marine NIS phyloseq object
saveRDS(phyloseq_NIS_18S, file.path(output_rds_directory, "Opua18S_NIS.rds")) 

# Decide not to use rarefy_even_depth going forward as sample depth looks consistent and reached biodiversity.
# Instead, will use percentage or relative abundance.
# Processing samples: relative abundance
# Remove reads less than or equal to 1 (i.e., singletons)
Opua18S_NIS.ps <- prune_samples(sample_sums(phyloseq_NIS_18S) > 0, phyloseq_NIS_18S)
# Calculate relative abundance
Opua18S_NIS.ps <- transform_sample_counts(Opua18S_NIS.ps, function(x) x / sum(x))

# For visualization, look at Species level
# Subset taxa by removing NAs in the 'Species' column and aggregate by Species for visualization of marine NIS
Opua18S_NIS.ps <- tax_glom(Opua18S_NIS.ps, taxrank = "Species") %>%
  subset_taxa(!is.na(Species))
