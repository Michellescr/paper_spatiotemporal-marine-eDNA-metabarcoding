# Load the installed packages
cat("Loading libraries")
library(dada2); packageVersion("dada2")
library(data.table)
library(phyloseq)
library(openssl)
library(theseus)
library(tidyverse)
library(Biostrings)
library(ggpubr)

# Environments
cutadapt = "/home/john/miniconda3/bin/cutadapt" # path to your cutadapt program. 
cores=18
path_main <- getwd()
path_fastq <- paste0(path_main,"/fastq_files")
path_results <- paste0(path_main,"/DADA2_results_18S")

# Create results directory if it doesn't exist
if(!dir.exists(path_results)) dir.create(path_results)

# 18S primer sequence: Zhang GK, Chain FJJ, Abbott CL, Cristescu ME. Metabarcoding using multiplexed markers increases species detection in complex zooplankton communities. Evol Appl. 2018 Sep 15;11(10):1901-1914. doi: 10.1111/eva.12694.
forward_primer=c("AGGGCAAKYCTGGTGCCAGC")
reverse_primer=c("GRCGGTATCTRATCGYCTT") 

# Taxonomic database path
taxo_db = "/srv/users/Oli/silva_132.18s.99_rep_set.dada2.fa.gz" # path to your reference database
norm = T # Whether the taxonomic assingment needs normalisation or not
sqlFile = "/home/john/database_files/accessionTaxa.sql"
ranks = c("Superkingdom", "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus", "Species")
keepSAR = T

#Truncation length and overlap
tL = c(225,216) # Truncation length for the forward and reverse reads
overlap = "--overlap 18" # Minimum overlap for primer matching with cutadapt

## Step 1
#  Demultiplexing and primer removal
cat("\nProcessing fastq files...\n")
fas_Fs_raw <- sort(list.files(path_fastq, pattern="R1_001.fastq.gz", full.names = TRUE))
fas_Rs_raw <- sort(list.files(path_fastq, pattern="R2_001.fastq.gz", full.names = TRUE))

fas_Fs_raw[2]
fas_Rs_raw[2]

# This is our set of primers
FWD <- forward_primer
REV <- reverse_primer
FWD_RC <- dada2:::rc(FWD)
REV_RC <- dada2:::rc(REV)

# Create cut directory
path_cut <- file.path(path_results, "cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)
fas_Fs_cut <- file.path(path_cut, basename(fas_Fs_raw))
fas_Rs_cut <- file.path(path_cut, basename(fas_Rs_raw))

R1_flags <- paste(paste("-g", FWD, collapse = " "), paste("-a", REV_RC, collapse = " "))
R2_flags <- paste(paste("-G", REV, collapse = " "), paste("-A", FWD_RC, collapse = " "))

# Run cutadapt for each sample
for(i in seq_along(fas_Fs_raw)) {
  cat("Processing", "-----------", i, "/", length(fas_Fs_raw), "-----------\n")
  system2(cutadapt, args = c(R1_flags, R2_flags,
                             "--discard-untrimmed",
                             "--max-n 0",
                             overlap,
                             "-o", fas_Fs_cut[i], "-p", fas_Rs_cut[i],
                             fas_Fs_raw[i], fas_Rs_raw[i]))
}

#Quality Assessment
out_1 <- ShortRead::qa(fas_Fs_raw)[["readCounts"]][,"read", drop = FALSE]

head(out_1)

#### Inspect read quality profiles ####
pF <- plotQualityProfile(sample(fas_Fs_cut, replace = FALSE, size = ifelse(length(fas_Fs_cut) < 100, length(fas_Fs_cut), 100)),aggregate = TRUE) + ggplot2::labs(title = "Forward")
pR <- plotQualityProfile(sample(fas_Rs_cut, replace = FALSE, size = ifelse(length(fas_Rs_cut) < 100, length(fas_Rs_cut), 100)),aggregate = TRUE)+ ggplot2::labs(title = "Reverse")
test = ggarrange(pF,pR, nrow = 2)
ggsave(filename = file.path(path_results, "Read_quality_profile_aggregated.pdf"), plot = test, width = 6, height = 8)


## Step 2: Quality filtering
# Place filtered files in filtered/ subdirectory
cat("\nPerforming quality filtering\n")
path_process <- path_cut # Path to processeed sequence
fnFs <- sort(list.files(path_process, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path_process, pattern="R2_001.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S\\d+"), `[`, 1)
filtFs <- file.path(path_results, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_results, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Quality Filtering
out_2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(0,0),
                       truncLen=tL,maxN=0, maxEE=c(2,2), truncQ=2,
                       rm.phix=TRUE,compress=TRUE, multithread=cores) # On Windows set multithread=FALSE
head(out_2,20)

## Step 3: Learn error rates
cat("\nLearning error rates\n")
filtFs <- paste0(file.path(path_results, "filtered", list.files(path=paste0(path_results,"/filtered"), pattern="_F_filt.fastq.gz")))
filtRs <- paste0(file.path(path_results, "filtered", list.files(path=paste0(path_results,"/filtered"), pattern="_R_filt.fastq.gz")))
errF <- learnErrors(filtFs, multithread=cores)
errR <- learnErrors(filtRs, multithread=cores)

#plot error rates
perrF <- plotErrors(errF, nominalQ = TRUE) + ggplot2::labs(title = "Error Forward")
perrR <- plotErrors(errR, nominalQ = TRUE) + ggplot2::labs(title = "Error Reverse")
fig = ggarrange(perrF,perrR, nrow = 2)

#save error rates
ggsave(filename = paste0(path_results,"/Error_rates_learning.pdf"), fig ,width = 6, height = 8)

## Step 4: Dereplication
exists <- file.exists(filtFs)
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(filtFs), "_F_filt"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(filtFs), "_F_filt"), `[`, 1)

## Step 5: Sample inference
cat("\nDenoising the data\n")
dadaFs <- dada(derepFs, err=errF, multithread=cores)
dadaRs <- dada(derepRs, err=errR, multithread=cores)

## Step 6: Merging forward and reverse reads
cat("\nMerging reads\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 0, verbose=TRUE)

## Step 7: Construct feature table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## Step 8: Remove chimeras
cat("\nRemoving chimeras\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=cores, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, file = paste0(path_results,"/seqtab.nochim.rds"))

# Save the non-chimeric sequence table
seqtab = readRDS(paste0(path_results,"/seqtab.nochim.rds"))

# Inspect distribution of sequence lengths
df = as.data.frame(table(nchar(getSequences(seqtab))))
df$Var1 = as.numeric(as.character(df$Var1))

# Basic density plot for sequence lengths
p <- ggplot(df, aes(x=Var1)) +
 geom_density()+
  #geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1)+
  ggthemes::theme_clean() +
  ylab("Density") +
  #geom_text(aes(x=mean(length), label=paste0("Mean\n",mean(length)), y=0.25)) +
  scale_x_continuous(name ="Sequence length (bp)",breaks = c(200,250,300,350,374,400,450), labels = c("200","250","300","350","374","400","450"))
cat("\nMean and median read length\n")
mean(df$Var1)
cat("\n")
median(df$Var1)
ggsave(filename =  paste0(path_results,"/Seq_lenght_density.pdf"), plot = p,device = "pdf")

## Step 9: Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- as.data.frame(cbind(out_2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)))

# Combine with initial counts
track = gdata::cbindX(out_1,track)
colnames(track) <- c("inputF", "demultiplexed", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
Sample_names <- as.data.frame(gsub("_S\\d+_L001.*","", rownames(track)))
head(track)
track = gdata::cbindX(Sample_names,track)
colnames(track)[1] = "Sample.name"

# Write track data to CSV
write.table(as.data.frame(track), paste0(path_results,"/Read_counts.csv"), sep=",", quote = F, row.names = F)

## Step 10: Assinging taxonomy
cat("\nAssingning taxonomy with NB classifier\n")
seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
# This uses a Naives-Bayes trained classifier
taxa <- assignTaxonomy(seqtab.nochim, taxo_db, multithread=cores)
taxa <- taxa[,!is.na(colnames(taxa))]
rownames(taxa) <-openssl::md5(rownames(taxa))
saveRDS(taxa, file = "taxaNB.rds")
taxa = taxa

## Step 11: Saving sequences
cat("\nSaving sequences\n")
uniquesToFasta(getUniques(seqtab.nochim), fout= paste0(path_results,"/uniqueSeqs.fasta"),ids=as.character(as.list(md5(names(getUniques(seqtab.nochim))))))
dna <- readDNAStringSet(paste0(path_results,"/uniqueSeqs.fasta")) # Create a DNAStringSet from the ASVs

## Step 12: Hands off to Phyloseq
cat("\nCreating a phyloseq object\n")
colnames(seqtab.nochim) <-md5(colnames(seqtab.nochim))
ps_run <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE),
  tax_table(taxa),
  refseq(dna))
  
  ps_run
  
# Save the phyloseq object
saveRDS(ps_run, file = paste0(path_results,"/ps_run.rds"))

# Prepare ASVs table for output  
otab = pstoveg_otu(ps_run) %>% t() %>% as.data.frame()
otab = cbind(ASVs = rownames(otab), otab)
write.table(otab, paste0(path_results,"/dada_table.txt"),quote=FALSE,sep="\t", row.names = FALSE)

# Save taxonomy table
taxa = cbind(ASVs = rownames(taxa), taxa)
write.table(taxa, paste0(path_results,"/tax_table.txt"),quote=FALSE,sep="\t", row.names = F)

## Taxo normalization
if(norm==T){
  ps_run = taxo_normalisation(obj = ps_run, sqlFile = sqlFile, ranks = ranks, keepSAR = keepSAR)
}

ps_run

saveRDS(ps_run, file = paste0(path_results,"/ps_run_norm.rds"))

# Save normalized ASVs table
otab = pstoveg_otu(ps_run) %>% t() %>% as.data.frame()
otab = cbind(ASVs = rownames(otab), otab)
write.table(otab, paste0(path_results,"/dada_table_norm.txt"),quote=FALSE,sep="\t", row.names = FALSE)

# Save normalized taxonomy table
taxa_df <- as.data.frame(tax_table(ps_run))
taxa_matrix <- as.matrix(taxa_df)
write.table(taxa_matrix, paste0(path_results,"/tax_table_norm.txt"),quote=FALSE,sep="\t", row.names = F)




