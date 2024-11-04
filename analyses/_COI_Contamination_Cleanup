### Note: The same process applies for COI and 18S datasets.


# Load the installed packages
cat("Loading libraries")
suppressPackageStartupMessages({
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(gdata)
})

# Display package versions
cat("Package Versions:\n")
cat("ggplot2:", packageVersion("ggplot2"), "\n")
cat("phyloseq:", packageVersion("phyloseq"), "\n")

# Set the working directory to your R scripts directory
script_directory <- "C:/Users/miche/OneDrive - Cawthron/Google_Drive/PhD Project/Chapter 3_dispersal of eDNA_eRNA/Data Analyze/Metabarcoding_data_2/R_script"
setwd(script_directory)

# Define the directory paths for input and output RDS files
input_rds_directory <- "C:/Users/miche/OneDrive - Cawthron/Google_Drive/PhD Project/Chapter 3_dispersal of eDNA_eRNA/Data Analyze/Metabarcoding_data_2/input_ps"
output_rds_directory <- "C:/Users/miche/OneDrive - Cawthron/Google_Drive/PhD Project/Chapter 3_dispersal of eDNA_eRNA/Data Analyze/Metabarcoding_data_2/input_ps"

# Define the directory path for output figures
fig_output_directory <- "C:/Users/miche/OneDrive - Cawthron/Google_Drive/PhD Project/Chapter 3_dispersal of eDNA_eRNA/Data Analyze/Metabarcoding_data_2/Figures"

# Load Data Set
# Load the phyloseq object and metadata
ps_run <- readRDS(file.path(input_rds_directory, "18S_ps_run.rds"))
map18S1 <- read.csv(file.path(input_rds_directory, "18Smetadata.csv"), header = TRUE, row.names = 1)

# Ensure all metadata columns are factors
map18S1[] <- lapply(map18S1, function(x) if(is.character(x)) as.factor(x) else x)

# Create a phyloseq object with merged taxa, OTU table, and sample data
Opua18S.ps <- phyloseq(otu_table(ps_run, taxa_are_rows = FALSE),
                        sample_data(map18S1),
                        tax_table(ps_run))

# Save the phyloseq object
saveRDS(Opua18S.ps, file.path(output_rds_directory, "Opua18S.ps.rds"))

# Addressing contamination in your processing
# Subset controls and filter taxa not found in the controls
Controls <- subset_samples(Opua18S.ps, Type %in% c("Extraction blank", "Field Blank", "NK", "Water"))
Controls <- phyloseq::filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
sample_sums(Controls)

# Plot taxonomic makeup of the controls
plot_bar(Controls, "Type", fill = "phylum") + geom_bar(stat = "identity")

# Subtraction of ASVs
# Subset the extraction controls
Extraction_neg <- subset_samples(Opua18S.ps, Type %in% c("Extraction blank", "Field Blank"))

# Find max value for each ASVs
Extraction_neg_max <- apply(as.data.frame(as.matrix(t(otu_table(Extraction_neg)))), 1, max)
Extraction_neg_sums <- colSums(otu_table(Extraction_neg))

# Export the ASV table from phyloseq into a data.frame
Extractiondf <- as.data.frame(as(otu_table(Opua18S.ps), "matrix"))

# Subtract values stored in the Extraction_neg_sums
Extractiondf[] <- sweep(Extractiondf, 2, Extraction_neg_sums)
Extractiondf <- replace(Extractiondf, Extractiondf < 0, 0)

# Create a new phyloseq object
Opua18S.subtractextract.ps <- phyloseq(otu_table(Extractiondf, taxa_are_rows = FALSE),
                                         sample_data(map18S1),
                                         tax_table(ps_run))

# Subset PCR & sequencing controls
PCR_neg <- subset_samples(Opua18S.subtractextract.ps, Type == "NK")
PCR_neg_max <- apply(as.data.frame(as.matrix(t(otu_table(PCR_neg)))), 1, max)
PCR_neg_sums <- colSums(otu_table(PCR_neg))

# Subtract PCR negatives
PCRdf <- as.data.frame(as(otu_table(Opua18S.subtractextract.ps), "matrix"))
PCRdf[] <- sweep(PCRdf, 2, PCR_neg_sums)
PCRdf <- replace(PCRdf, PCRdf < 0, 0)

# Create a new phyloseq object for PCR
Opua18S.subtractextractnk.ps <- phyloseq(otu_table(PCRdf, taxa_are_rows = FALSE),
                                           sample_data(map18S1),
                                           tax_table(ps_run))

# Subset sequencing controls
Sequencing_neg <- subset_samples(Opua18S.subtractextractnk.ps, Type == "Water")
Sequencing_neg_max <- apply(as.data.frame(as.matrix(t(otu_table(Sequencing_neg)))), 1, max)
Sequencing_neg_sums <- colSums(otu_table(Sequencing_neg))

# Subtract sequencing negatives
Sequencingdf <- as.data.frame(as(otu_table(Opua18S.subtractextractnk.ps), "matrix"))
Sequencingdf[] <- sweep(Sequencingdf, 2, Sequencing_neg_sums)
Sequencingdf <- replace(Sequencingdf, Sequencingdf < 0, 0)

# Create a new phyloseq object for sequencing
Opua18S.subtract.ps <- phyloseq(otu_table(Sequencingdf, taxa_are_rows = FALSE),
                                  sample_data(map18S1),
                                  tax_table(ps_run))

# Calculate the percentage of reads retained after subtraction
org.ss <- as.data.frame(sample_sums(Opua18S.ps))
org.ss$Names <- rownames(org.ss)
new.ss <- as.data.frame(sample_sums(Opua18S.subtract.ps))
new.ss$Names <- rownames(new.ss)
control.rem.maxsub <- dplyr::left_join(org.ss, new.ss, by = "Names")
colnames(control.rem.maxsub) <- c("Original", "Names", "New")
control.rem.maxsub.final <- control.rem.maxsub %>%
  mutate(perc = New / Original * 100)
print(control.rem.maxsub.final)

# Remove control samples and singletons from the dataset
Opua18S.subtract.ps <- subset_samples(Opua18S.subtract.ps, Type == "Sample")
Opua18S.subtract.ps <- prune_taxa(taxa_sums(Opua18S.subtract.ps) > 1, Opua18S.subtract.ps)

# Summarize and save the results
summary(Opua18S.subtract.ps)
saveRDS(Opua18S.subtract.ps, file.path(output_rds_directory, "Opua18S.subtract.ps.rds"))

# Note: The same process applies for COI and 18S datasets.
