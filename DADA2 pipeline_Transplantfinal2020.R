#Dada2 pipeline for analysis of samples from Singapore P. acuta Reciprocol Transplantation Experiment

#Create your own 'Samle.csv' metadata files for the phyloseq object

##########################################################################
#Below is the code for running Dada2, convert to phyloseq, and decontamination with sample blanks


library(dada2); packageVersion("dada2")
path <- 'C:/Users/SET_LOCAL_PATH'
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
## It is common in Illumina sequencing that the reverse reads are worse than the forward reads

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## Based on the Quality Profile graphs, adjust the truncLen
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=10, truncLen=c(200,170),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

##LEARN Error rates for forward and reverse sequences
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#Check that error rate goes down as quality score increases
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## Dereplication combines all identical sequencing reads into into “unique sequences” 
#with a corresponding “abundance” equal to the number of reads with that unique sequence

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Apply the sample inference algorithm to the derep data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# then inspect it if you want to see
dadaFs[[1]]

## Merge forward and reverse sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Contruct sequence table
## This is where you decide the cutoff length of your sequnces.
# Choose based on what it should be and the length of the majority of the sequences

seqtab <- makeSequenceTable(mergers)
## The sequences being tabled vary in length.
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Here you see the distribution of the different sequence lengths
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(271,273)]


seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## TO confirm that all the steps worked,and you didn't lose too much in any one step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, file = "track_read_counts_at_each_step.csv", row.names = TRUE)

#########################################################
#Remore contaminants with Decontam

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam)

# import sample data ####
sample = read.csv("sample.csv")
row.names(sample) <- sample$Sample
row.names(seqtab.nochim)


# Find controlsamples (extraction negatives) ####
sample$controls <- sample$Type == "Blank"

# find contaminants
contams = isContaminant(seqtab.nochim, neg = sample$controls, normalize = TRUE)
table(contams$contaminant)
write.csv(contams, file = "likely_contaminants.csv", row.names = TRUE)


# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[sample$controls == FALSE,]
sample = sample[sample$controls == FALSE,]


#Assign taxonomy with IDTAXA algorithm (classification performance reported to be better than the long-time standard set by the naive Bayesian classifier)
library(DECIPHER); packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/tax/IDTaxa/silva_rn99_v138_train_set.fa.gz") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

taxa <- taxid

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


#To save a copy of the data 
write.csv(seqtab.nochim, "seqtab.nochim.csv")
write.csv(taxa, "taxa.csv")


#Move to phyloseq to remove untargetted sequences (mito, chloroplasts, unassigned sequences)
#Call the files up with the sample metadata csv file that I created

library(phyloseq)

## This loads/organizes/prepares the files that will be put in the phyloseq object
LarvaeOTU_table <- read.table('seqtab.nochim.csv', header = T, row.names = 1, check.names = F, sep = ",")
Larvaetaxa <- as.matrix(read.table('taxa.csv', header = T, row.names = 1, check.names = F, sep = ","))
Larvesample <- read.table('Sample.csv', header = T, row.names = 1, check.names = F, sep = ",")

## Make each phyloseq appropriate object, match col/row names
LarvaeOTU_table_ps <- otu_table(LarvaeOTU_table, taxa_are_rows = F)
Larvesample_ps <- sample_data(Larvesample)
Larvaetaxa_ps <- tax_table(Larvaetaxa)

#Can view the different row and col headings to make sure everything is matching
colnames(LarvaeOTU_table_ps)
rownames(Larvesample_ps)
rownames(LarvaeOTU_table_ps)
rownames(Larvaetaxa_ps)

## Merge all of the three prepared files into a single phyloseq object
Transplant_ps <- phyloseq(LarvaeOTU_table_ps, Larvaetaxa_ps, Larvesample_ps)

#Clean up with the ASV names and replaces the long sequence reads from the rows/cols 
dna <- Biostrings::DNAStringSet(taxa_names(Transplant_ps))
names(dna) <- taxa_names(Transplant_ps)
Transplant_ps <- merge_phyloseq(Transplant_ps, dna)
taxa_names(Transplant_ps) <- paste0("ASV", seq(ntaxa(Transplant_ps)))

##Check your result
Transplant_ps


table(tax_table(Transplant_ps)[, "phylum"], exclude = NULL)
table(tax_table(Transplant_ps)[, "family"], exclude = NULL)

# subset ps to remove mitochondrial, chloroplast, and unassigned sequences####
ps <- subset_taxa(Transplant_ps, !is.na(domain) & !domain %in% c("", "NA"))

# Remove "mitochondria" taxa
filterFamily = c("Mitochondria")
ps1 = subset_taxa(ps, !family %in% filterFamily)

# Remove "Chloroplast" taxa
filterOrder = c("Chloroplast")
ps2 = subset_taxa(ps1, !order %in% filterOrder)


#To save a copy of the data files

ps2_asvtable <- otu_table(ps2)
write.csv(ps2_asvtable, file='ps2_asvtable_finalnocontam.csv')
ps2_asvtax <- tax_table(ps2)
write.csv(ps2_asvtax, file='ps2_asvtax_finalnocontam.csv')

#############################################################################
#Run the DMM models (used to test for temporal bias in sample processing) 

library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)

#This loads/organizes/prepares the files that will be put in the phyloseq object
TransplantOTU_table <- read.table('ps2_asvtable_finalnocontam.csv', header = T, row.names = 1, check.names = F, sep = ",")
Transplanttaxa <- as.matrix(read.table('ps2_asvtax_finalnocontam.csv', header = T, row.names = 1, check.names = F, sep = ","))
Transplantsample <- read.table('sample.csv', header = T, row.names = 1, check.names = F, sep = ",")

#Make each one the phyloseq appropriate object, match col/row names
TransplantOTU_table_ps <- otu_table(TransplantOTU_table, taxa_are_rows = F)
Transplantsample_ps <- sample_data(Transplantsample)
Transplanttaxa_ps <- tax_table(Transplanttaxa)

#Merge all of the three prepared files into a single phyloseq object
Transplant_ps_final <- phyloseq(TransplantOTU_table_ps, Transplanttaxa_ps, Transplantsample_ps)

#Subset to specific data to analyze
Time00subset_ps_final = subset_samples(Transplant_ps_final, Timepoint=="0_day")
Time00subset_ps_final
Timesubset_coral_ps_final = subset_samples(Time00subset_ps_final, Type=="Coral")
Timesubset_coral_ps_final

pseq <- Timesubset_coral_ps_final

pseq.comp <- microbiome::transform(pseq, "compositional")
taxa <- core_members(pseq.comp, detection = 0.1/100, prevalence = 50/100)
pseq <- prune_taxa(taxa, pseq)


# Pick the OTU count matrix
# and convert it into samples x taxa format
dat <- abundances(pseq)
count <- as.matrix(t(dat))
fit <- mclapply(1:6, dmn, count = count, verbose=TRUE)
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
best <- fit[[which.min(lplc)]]
ass <- apply(mixture(best), 1, which.max)
ass
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}


###########################################
#Start analysis with final dataset
#Create rarefaction curves to determine where to set rarefy minimum

sample_sums(Transplant_ps_final)
rarecurve((otu_table(Transplant_ps_final)), step=50, cex=0.5)

Transplant_ps_final_rare <- rarefy_even_depth(Transplant_ps_final, sample.size = 8416,
                 replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

sample_sums(Transplant_ps_final_rare)

#Save the final rarefied data as csv files
Transplant_asv_final_rare <- otu_table(Transplant_ps_final_rare)
write.csv(Transplant_asv_final_rare, file='Transplant_asv_final_rare.csv')

Transplant_asvtax_final_rare <- tax_table(Transplant_ps_final_rare)
write.csv(Transplant_asvtax_final_rare, file='Transplant_asvtax_final_rare.csv')


##########################################
#Calculate and plot alpha diversity metrics of rarefied data

Transplant_alpha <- estimate_richness(Transplant_ps_final_rare, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
Transplant_alpha


write.csv(Transplant_alpha, file='Transplant_alphadiversity.csv')

alpha <- read.csv("Transplant_alphadiversity.csv", header = T)
str(alpha)


results<-aov(alpha$Richness ~ alpha$Transplantation)
results<-aov(alpha$Chao1 ~ alpha$Transplantation)
results<-aov(alpha$Shannon ~ alpha$Transplantation)
results<-aov(alpha$InvSimpson ~ alpha$Transplantation)

summary(results)
TukeyHSD(results)


library(ggplot2)

Richness <- ggplot(alpha, aes(x=Transplantation, y=Richness)) + 
  geom_boxplot()

Chao1 <- ggplot(alpha, aes(x=Transplantation, y=Chao1)) + 
  geom_boxplot()

Shannon <- ggplot(alpha, aes(x=Transplantation, y=Shannon)) + 
  geom_boxplot()

InvSimpson <- ggplot(alpha, aes(x=Transplantation, y=InvSimpson)) + 
  geom_boxplot()

library(cowplot)
cowplot::plot_grid(Richness, Chao1, Shannon, InvSimpson, labels = c("A", "B", "C","D"))


##########################################
#Calculate relative abundance by phylum

psrelabund <- tax_glom(Transplant_ps_final_rare, "phylum")
psrelabund0 <- transform_sample_counts(psrelabund, function(x) x / sum(x))
psrelabund1 <- merge_samples(psrelabund0, "Group")
psrelabund2 <- transform_sample_counts(psrelabund1, function(x) x / sum(x))
plot_bar(psrelabund2, fill="phylum")

summarize_taxa(psrelabund2, "phylum")

psrelabund2_list <- otu_table(psrelabund2)
write.csv(psrelabund2_list, file='RelativeAbundance.csv')
psrelabund2_asvtax_final_rare <- tax_table(psrelabund2)
write.csv(psrelabund2_asvtax_final_rare, file='RelativeAbundance_asvtax.csv')

#Make Relative abundance plot
library(ggplot2)
library(reshape2)

RA <- read.csv("RelativeAbundance.csv", header = T)
RA

RA$Sample <- factor(RA$Sample, levels = c("Resident Raffles", "Resident Kusu", "Transplant KR", "Transplant RK", "Seawater"))

df_long <- melt(RA, id.vars = "Sample", variable.name = "Phyla")

ggplot(df_long, aes(x = Sample, y = value, fill = Phyla)) + 
  geom_bar(stat = "identity", width=0.6)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())

