### Pre-Processing for Tetiaroa Atoll Restoration Project: 16S rRNA gene sequencing data November 2021
#Contributing to the manuscript on coral reef conservation across scales 

#This script includes both the water and coral microbial data

#set working directory to procedure folder

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/TARP_motu_comparison/microbiome_analyses/downstream_analyses/pre-processing/procedure/")

library(qiime2R)
library(phyloseq)
library(decontam)
library(ggplot2)

##Pre-processing for phyloseq
##We start with the coral data

#Upload all files from qiime2 into phyloseq
physeq.coral <- qza_to_phyloseq("../../../bioinformatics/coral_nov21/output/table-trim-filtered.qza",
                          "../../../bioinformatics/coral_nov21/output/rooted-tree-trim.qza", 
                          "../../../bioinformatics/coral_nov21/output/taxonomy-trim.qza", 
                          "../../../../metadata/coral_nov21_metadata.txt")

physeq.water <- qza_to_phyloseq("../../../bioinformatics/water_nov21/output/table-trim-filtered.qza",
                                "../../../bioinformatics/water_nov21/output/rooted-tree-trim.qza", 
                                "../../../bioinformatics/water_nov21/output/taxonomy-trim.qza", 
                                "../../../../metadata/water_nov21_metadata.txt")

#Check the taxonomic classification is correct
rank_names(physeq.coral)
rank_names(physeq.water)
tax_table(physeq.coral) #6752 ASVs by 7 taxonomic ranks
tax_table(physeq.water) #26893 ASVs by 7 taxonomic ranks

#Remove eukaryotes and Unassigned reads
#And chloroplasts - maybe check bioinformatic pipeline for this!
physeq.coral <- subset_taxa(physeq.coral, Kingdom != "d_0__Eukaryota")
physeq.coral <- subset_taxa(physeq.coral, Kingdom != "Eukaryota")
physeq.coral <- subset_taxa(physeq.coral, Kingdom != "Unassigned")
physeq.coral <- subset_taxa(physeq.coral, Order != "Chloroplast")

physeq.water <- subset_taxa(physeq.water, Kingdom != "d_0__Eukaryota")
physeq.water <- subset_taxa(physeq.water, Kingdom != "Eukaryota")
physeq.water <- subset_taxa(physeq.water, Kingdom != "Unassigned")
physeq.water <- subset_taxa(physeq.water, Order != "Chloroplast")

#Remove singletons before further filtering to minimize the manual work on removing mitochondria
physeq.coral <- prune_taxa(taxa_sums(physeq.coral) > 1, physeq.coral) 
physeq.water <- prune_taxa(taxa_sums(physeq.water) > 1, physeq.water)
#Check your summary stats:
print(microbiome::summarize_phyloseq(physeq.coral)) 
#total reads: 1,986,697
#avg reads/sample: 26,489.293
#median reads: 20,474
#Min # of reads: 69
#Max # of reads: 145,694
#Sparsity: 0.98
print(microbiome::summarize_phyloseq(physeq.water))
#total reads: 3,261,758
#avg reads/sample: 41,817.41
#median reads: 43,767
#Min # of reads: 1
#Max # of reads: 78,496
#Sparsity: 0.97

#Check for mitochondria & blast seqs to make sure you are not removing bacterial taxa
##Family == Mitochondria are only in the order Rickettsiales, so we check on family & genus levels
#And remove accordingly
mito <- subset_taxa(physeq.coral, Order == "Rickettsiales")
mito.taxa <- as.data.frame(tax_table(mito))
View(mito.taxa)
#need to remove species: Nitzschia_sp., Berkeleya_fennica, Ceramium_japonicum, NA
physeq.coral <- subset_taxa(physeq.coral, Genus != "Mitochondria")

#For water
mito <- subset_taxa(physeq.water, Order == "Rickettsiales")
mito.taxa <- as.data.frame(tax_table(mito))
View(mito.taxa)

physeq.water <- subset_taxa(physeq.water, Genus != "Mitochondria")

#Need sample data for following steps to remove contaminants 
####First for coral####
sample.data <- as(sample_data(physeq.coral), "data.frame")

#Remove contaminants
#Inspect library size
sample.data$LibrarySize <- sample_sums(physeq.coral)
sample.data <- sample.data[order(sample.data$LibrarySize),]
sample.data$Index <- seq(nrow(sample.data))
#visualize
ggplot(data = sample.data, aes(x=Index, y=LibrarySize, color = sample.type)) +
  geom_point()

#Next check for contaminants using prevalence and threshold 0.5 (more conservative)
sample_data(physeq.coral)$is.neg <- sample_data(physeq.coral)$sample.type == "blank"
contamdf.prev <- isContaminant(physeq.coral, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant)) #only six contaminants

contamdf.prev.contamsonly <- contamdf.prev %>% filter(contaminant == TRUE)
View(contamdf.prev.contamsonly)

#Identify the contaminants for reporting
bad.taxa <- rownames(contamdf.prev.contamsonly)
bad.taxa <- c("3488aa431e41813ca40c335c5ff6ff36","694df3c7f8b6b66c922ed51a965d75d0")

all.taxa <- taxa_names(physeq.coral)
contaminant.taxa <- all.taxa[(all.taxa %in% bad.taxa)]
physeq.contams.only <- prune_taxa(contaminant.taxa, physeq.coral)
physeq.contams.only.df <- as.data.frame(tax_table(physeq.contams.only))
print(physeq.contams.only.df)
#3488aa431e41813ca40c335c5ff6ff36	d__Bacteria	Verrucomicrobiota	Chlamydiae	Chlamydiales	Simkaniaceae	uncultured	jellyfish_metagenome
#694df3c7f8b6b66c922ed51a965d75d0	d__Bacteria	Proteobacteria	Gammaproteobacteria	Oceanospirillales	Endozoicomonadaceae	Endozoicomonas

#Remove contaminants
physeq.noncont.coral <- prune_taxa(!contamdf.prev$contaminant, physeq.coral)

#Check final numbers
print(microbiome::summarize_phyloseq(physeq.noncont.coral)) 
#Total reads: 412,443
#Average reads/sample: 5,499.24
#Min # of reads: 0
#Max # of reads: 27,915
tax_table(physeq.noncont.coral)  #Total taxa: 3617
sample_data(physeq.noncont.coral) #Total samples:143 samples
#Count Sequences/sample
print(totalreads <- sort(sample_sums(physeq.noncont.coral))) #lowest workable read count is 1,453 - lost three samples to failed sequencing: TC169, TC196, TC214

#Remove blanks & mock because we already dealt with potential contaminants
physeq.noncont.coral <- subset_samples(physeq.noncont.coral, is.neg != "TRUE")

#Prune singletons (these are reads that are only found once) after all filtering steps 
physeq.prune.coral <- prune_taxa(taxa_sums(physeq.noncont.coral) > 1, physeq.noncont.coral)
print(microbiome::summarize_phyloseq(physeq.prune.coral))
#Total reads: 412,424
#Avg reads/sample: 5,891.8

#Save the pruned file as an RDS so that all subsequent R scripts can call this pre-processed data
saveRDS(physeq.prune.coral, "../output/physeq-coral-nov21.RDS")



####Second for water####
sample.data <- as(sample_data(physeq.water), "data.frame")
head(sample.data)
#Remove contaminants
#Inspect library size
sample.data$LibrarySize <- sample_sums(physeq.water)
sample.data <- sample.data[order(sample.data$LibrarySize),]
sample.data$Index <- seq(nrow(sample.data))
#visualize
ggplot(data = sample.data, aes(x=Index, y=LibrarySize, color = sample_type)) +
  geom_point()

#Next check for contaminants using prevalence and threshold 0.5 (more conservative)
sample_data(physeq.water)$is.neg <- sample_data(physeq.water)$sample_type == "blank"
contamdf.prev <- isContaminant(physeq.water, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant)) #only six contaminants

contamdf.prev.contamsonly <- contamdf.prev %>% filter(contaminant == TRUE)
View(contamdf.prev.contamsonly)

#Identify the contaminants for reporting
bad.taxa <- rownames(contamdf.prev.contamsonly)
bad.taxa <- c("12fb659ec96872c9175ccfa0e8f43b0c","e57c5df6a9b3b982472e7754ed31f313", "ff9d93d7b7e46787568f2d241caeaf3b", "65d43491988bfe557da4d86a5ba25dae")

all.taxa <- taxa_names(physeq.water)
contaminant.taxa <- all.taxa[(all.taxa %in% bad.taxa)]
physeq.contams.only <- prune_taxa(contaminant.taxa, physeq.water)
physeq.contams.only.df <- as.data.frame(tax_table(physeq.contams.only))
print(physeq.contams.only.df)
#12fb659ec96872c9175ccfa0e8f43b0c d__Bacteria Bdellovibrionota     Bdellovibrionia Bacteriovoracales  Bacteriovoraceae  Halobacteriovorax
#e57c5df6a9b3b982472e7754ed31f313 d__Bacteria   Proteobacteria Gammaproteobacteria   Burkholderiales  Oxalobacteraceae  Massilia
#ff9d93d7b7e46787568f2d241caeaf3b d__Bacteria   Proteobacteria Gammaproteobacteria   Pseudomonadales  Pseudomonadaceae  Pseudomonas
#65d43491988bfe557da4d86a5ba25dae d__Bacteria       Firmicutes             Bacilli  Staphylococcales  Staphylococcaceae Staphylococcus


#Remove contaminants
physeq.noncont.water <- prune_taxa(!contamdf.prev$contaminant, physeq.water)

#Check final numbers
print(microbiome::summarize_phyloseq(physeq.noncont.water)) 
#Total reads: 2,980,787
#Average reads/sample: 38,215.22
#Min # of reads: 1
#Max # of reads: 73,154
tax_table(physeq.noncont.water)  #19,798 taxa by 7 taxonomic ranks
sample_data(physeq.noncont.water)
#Count Sequences/sample
print(totalreads <- sort(sample_sums(physeq.noncont.water))) #lowest workable read count is 17,959 - lost one sample as a result of same sequence #s as blanks: TW205 & TW210

#Remove blanks & mock because we already dealt with potential contaminants
physeq.noncont.water <- subset_samples(physeq.noncont.water, is.neg != "TRUE")

#Prune singletons (these are reads that are only found once) after all filtering steps 
physeq.prune.water <- prune_taxa(taxa_sums(physeq.noncont.water) > 1, physeq.noncont.water)
print(microbiome::summarize_phyloseq(physeq.prune.water))
#Total reads: 2,961,014
#Avg reads/sample: 40,561

#Save the pruned file as an RDS so that all subsequent R scripts can call this pre-processed data
saveRDS(physeq.prune.water, "../output/physeq-water-nov21.RDS")

