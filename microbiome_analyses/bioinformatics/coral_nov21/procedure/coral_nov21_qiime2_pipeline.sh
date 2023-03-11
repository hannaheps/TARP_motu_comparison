##Note: Rawdata are not kept on github due to file size limitations
##These will be available publicly via the SRA on publication

#Find & change to your working directory

#Need this other filepath for the first qiime2 code so that we can access raw data and save files too large for github elsewhere
cd ~/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/core_analysis/island_survey_june2021_current


conda activate qiime2-2022.2

#Import into qiime2

##Need to create a manifest file to direct qiime to each data file - view this in "input" directory
 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../../TARP_motu_comparison/microbiome_analyses/bioinformatics/coral_nov21/input/coral_nov21_manifest.txt \
  --output-path output/bioinformatics/coral-nov21-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
  
qiime demux summarize \
  --i-data output/bioinformatics/coral-nov21-paired-end-demux.qza  \
  --o-visualization ../../TARP_motu_comparison/microbiome_analyses/bioinformatics/coral_nov21/output/coral-nov21-paired-end-demux.qzv
  
#qiime tools view paired-end-demux.qzv

#Having an issue with differences in length in the reads (probably bc there are host sequences in here too, but it's a big pain!)
#First trial cutadapt to remove the exact sequences of the primers and go from there
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences output/bioinformatics/coral-nov21-paired-end-demux.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --o-trimmed-sequences output/bioinformatics/coral-nov21-paired-end-demux-trim.qza

qiime demux summarize \
  --i-data output/bioinformatics/coral-nov21-paired-end-demux-trim.qza  \
  --o-visualization ../../TARP_motu_comparison/microbiome_analyses/bioinformatics/coral_nov21/output/coral-nov21-paired-end-demux-trim.qzv

##Here we can now change to the current directory in which these files are located for the rest of the pipeline
cd ~/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/TARP_motu_comparison/microbiome_analyses/bioinformatics/coral_nov21/procedure


##Denoising using dada2
#The trim length removes primers 
#The truncation length is based on the visualization of the quality scores
#by basepair location. This can be viewed using paired-end-demux.qza in view.qiime2.org
#This step can take a long time/a lot of computing power depending on sample number & read depth
  #--p-trim-left-f 19 \
  #--p-trim-left-r 20 \

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ~/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/core_analysis/island_survey_june2021_current/output/bioinformatics/coral-nov21-paired-end-demux-trim.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --o-table ../output/table-trim.qza \
  --o-representative-sequences ../output/rep-seqs-trim.qza \
  --o-denoising-stats ../output/denoising-stats-trim.qza
  
#Visualize the DADA2 output: 
#In Qiime view you can save each of these as readable versions e.g., tsv or fasta
qiime metadata tabulate \
  --m-input-file ../output/denoising-stats-trim.qza \
  --o-visualization ../output/denoising-stats-trim.qzv
  
qiime feature-table summarize \
  --i-table ../output/table-trim.qza  \
  --o-visualization ../output/table-trim.qzv  \
  --m-sample-metadata-file ../../../../metadata/coral_nov21_metadata.txt

qiime feature-table tabulate-seqs \
  --i-data ../output/rep-seqs-trim.qza \
  --o-visualization ../output/rep-seqs-trim.qzv

#Taxonomic assignment
##SILVA database too large to keep on GitHub. Please download latest version & use following code
#with respective fasta and taxonomy files. 

#Here I make a silva classifier using SILVA v 138 (Quast et al.)

##The v 138.1 silva is apparently complicated to get into Qiime format.
#You can download qiime2 formatted SILVA 138.1 database here: https://docs.qiime2.org/2022.2/data-resources/
#These documents are already in .qza format, so we start by making a classifier and then training on our own data

#These files are large, and are not kept on github repository. Please adjust code for your specific file structure
#qiime feature-classifier fit-classifier-naive-bayes \
#  --i-reference-reads ../../../../../../../../../Biocomputing/databases/SILVA_138/silva-138-99-seqs-515-806.qza \
#  --i-reference-taxonomy ../../../../../../../../../Biocomputing/databases/SILVA_138/silva-138-99-tax.qza \
#  --o-classifier ../../../../../../../../../Biocomputing/databases/SILVA_138/classifier-silva-138.qza

#Train the classifier on your own data
qiime feature-classifier classify-sklearn \
  --i-classifier ../../../../../../../../../Biocomputing/databases/SILVA_138/classifier-silva-138.qza \
  --i-reads ../output/rep-seqs-trim.qza \
  --o-classification ../output/taxonomy-trim.qza
  
#Trial the following code to see if the taxonomy.qza file bug is fixed in QIIME2-2022.2
#It was!

#This code spits out a table of features, or ASVs, present in sample set
#You can view by dragging the .qzv file into view.qiime2.org
#Or use the "qiime tools view" command
qiime metadata tabulate \
  --m-input-file ../output/taxonomy-trim.qza \
  --o-visualization ../output/taxonomy-trim.qzv

### Filtering your table: removing mitochondria and chloroplast reads, PCR/sequencing errors, etc###

##Filter table
##I don't filter out mitochondria here bc it is a sticky category - SILVA database likes
#to classify some bacteria as Family:Mitochondria (e.g., Rickettsiales) so I like to leave it in &
#remove it manually using BLASTn on all mitochondria reads
#But I do remove eukaryotes and chloroplasts! 

qiime taxa filter-table \
  --i-table ../output/table-trim.qza \
  --i-taxonomy ../output/taxonomy-trim.qza \
  --p-exclude Eukaryota \
  --o-filtered-table ../output/table-trim-filtered.qza
   
qiime feature-table summarize \
  --i-table ../output/table-trim-filtered.qza \
  --o-visualization ../output/table-trim-filtered.qzv
  
##Create a phylogenetic tree (With large datasets this takes ALOT of RAM - give it 50G)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ../output/rep-seqs-trim.qza \
  --o-alignment ../output/aligned-rep-seqs-trim.qza \
  --o-masked-alignment ../output/masked-aligned-rep-seqs-trim.qza \
  --o-tree ../output/unrooted-tree-trim.qza \
  --o-rooted-tree ../output/rooted-tree-trim.qza

##Check for alpha rarefaction
##NOTE:make sure your metadata file does NOT have a variable column named "depth"
##Otherwise it will return an error: "Plugin error from diversity: cannot insert depth, already exists"
##as the code will add a depth variable for rarefaction
##Thus, I changed the metadata column to "water_depth" so that it does not interfere

qiime diversity alpha-rarefaction \
  --i-table ../output/table-trim-filtered.qza \
  --i-phylogeny ../output/rooted-tree-trim.qza \
  --p-max-depth 10000 \
  --m-metadata-file ../../../../metadata/coral_nov21_metadata.txt \
  --o-visualization ../output/alpha-rarefaction-trim.qzv

qiime taxa barplot \
  --i-table ../output/table-trim-filtered.qza \
  --i-taxonomy ../output/taxonomy-trim.qza \
  --m-metadata-file ../../../../metadata/coral_nov21_metadata.txt \
  --o-visualization ../output/barplot-trim.qzv

##For downstream analyses you need the following four files:
#1. table-filtered-trim.qza
#2. rooted-tree-trim.qza
#3. taxonomy-trim.qza
#4. coral_nov21_transect_metadata.txt


