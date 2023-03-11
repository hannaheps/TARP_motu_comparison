### Pre-Processing for Epstein et al. Guano Pilot Experiment, Tetiaroa French Polynesia ###

#set working directory to procedure folder

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/core_analysis/island_survey_june2021_current/procedure/")

library(qiime2R)
library(phyloseq)
library(decontam)
library(ggplot2)

##Pre-processing for phyloseq
#To start with provided RDS file, go to line 81

#Upload all files from qiime2 into phyloseq
physeq <- qza_to_phyloseq("../output/bioinformatics/table-trim-filtered.qza", 
                          "../output/bioinformatics/rooted-tree-trim.qza", 
                          "../output/bioinformatics/taxonomy-trim.qza", 
                          "../metadata/tet_transect_2021_metadata.txt")


#Check the taxonomic classification is correct
rank_names(physeq)
tax_table(physeq) #10285 ASVs by 7 taxonomic ranks

#Remove eukaryotes and Unassigned reads
#And chloroplasts - maybe check bioinformatic pipeline for this!
physeq <- subset_taxa(physeq, Kingdom != "d_0__Eukaryota")
physeq <- subset_taxa(physeq, Kingdom != "Eukaryota")
physeq <- subset_taxa(physeq, Kingdom != "Unassigned")
physeq <- subset_taxa(physeq, Order != "Chloroplast")

#Remove singletons before further filtering to minimize the manual work on removing mitochondria
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq) #9290 ASVs

#Check your summary stats:
print(microbiome::summarize_phyloseq(physeq)) 
#total reads: 3,215,030
#avg reads/sample: 22,482.7
#median reads: 18281
#Min # of reads: 69
#Max # of reads: 145,577
#Sparsity: 0.99

#Check for mitochondria & blast seqs to make sure you are not removing bacterial taxa
mito <- subset_taxa(physeq, Family == "Mitochondria")
mito.taxa <- as.data.frame(tax_table(mito))
View(mito.taxa)
#All rickettsiales

#Need sample data for following steps to remove contaminants
sample.data <- as(sample_data(physeq), "data.frame")

#Remove contaminants
#Inspect library size
sample.data$LibrarySize <- sample_sums(physeq)
sample.data <- sample.data[order(sample.data$LibrarySize),]
sample.data$Index <- seq(nrow(sample.data))
#visualize
ggplot(data = sample.data, aes(x=Index, y=LibrarySize, color = sample.type)) +
  geom_point()

#Next check for contaminants using prevalence and threshold 0.5 (more conservative)
sample_data(physeq)$is.neg <- sample_data(physeq)$sample.type == "blank"
contamdf.prev <- isContaminant(physeq, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant)) #only six contaminants

contamdf.prev.contamsonly <- contamdf.prev %>% filter(contaminant == TRUE)
View(contamdf.prev.contamsonly)

#Identify the contaminants for reporting
bad.taxa <- rownames(contamdf.prev.contamsonly)
bad.taxa <- c("3488aa431e41813ca40c335c5ff6ff36","23b654cbbea7606890cffdf88b65101f","65cabee7f6aacfe085f9e0d30ebc32ca",
              "694df3c7f8b6b66c922ed51a965d75d0", "3ebe761bfb1238c87195d431f41bf976", "ff9d93d7b7e46787568f2d241caeaf3b")

all.taxa <- taxa_names(physeq)
contaminant.taxa <- all.taxa[(all.taxa %in% bad.taxa)]
physeq.contams.only <- prune_taxa(contaminant.taxa, physeq)
physeq.contams.only.df <- as.data.frame(tax_table(physeq.contams.only))
print(physeq.contams.only.df)
#3488aa431e41813ca40c335c5ff6ff36	d__Bacteria	Verrucomicrobiota	Chlamydiae	Chlamydiales	Simkaniaceae	uncultured	jellyfish_metagenome
#23b654cbbea7606890cffdf88b65101f	d__Bacteria	Proteobacteria	Alphaproteobacteria	SAR11_clade	Clade_I Clade_Ia
#65cabee7f6aacfe085f9e0d30ebc32ca	d__Bacteria	Proteobacteria	Gammaproteobacteria	Oceanospirillales	Endozoicomonadaceae	Endozoicomonas
#694df3c7f8b6b66c922ed51a965d75d0	d__Bacteria	Proteobacteria	Gammaproteobacteria	Oceanospirillales	Endozoicomonadaceae	Endozoicomonas
#3ebe761bfb1238c87195d431f41bf976	d__Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales	Moraxellaceae	Acinetobacter
#ff9d93d7b7e46787568f2d241caeaf3b	d__Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales	Pseudomonadaceae	Pseudomonas

#Remove contaminants
physeq.noncont <- prune_taxa(!contamdf.prev$contaminant, physeq)

#Check final numbers
print(microbiome::summarize_phyloseq(physeq.noncont)) 
#Total reads: 996,529
#Average reads/sample: 6,968.7
#Min # of reads: 6
#Max # of reads: 57,302
tax_table(physeq.noncont)  #Total ASVs:9,284 
sample_data(physeq.noncont) #Total samples:143 samples
#Count Sequences/sample
totalreads <- sample_sums(physeq.noncont) #lowest workable read count is 1,263 - lost five samples to failed sequences: TC88 TC109 TC145 TC169 and TC196

#Remove blanks & mock because we already dealt with potential contaminants
physeq.noncont <- subset_samples(physeq.noncont, is.neg != "TRUE")

#Prune singletons (these are reads that are only found once) after all filtering steps 
physeq.prune <- prune_taxa(taxa_sums(physeq.noncont) > 1, physeq.noncont)
print(microbiome::summarize_phyloseq(physeq.prune))
#Total reads: 990,348
#Avg reads/sample: 7,176.43

#Save the pruned file as an RDS so that all subsequent R scripts can call this pre-processed data
saveRDS(physeq.prune, "../output/downstream_analyses/island-survey-2021-pruned.RDS")
