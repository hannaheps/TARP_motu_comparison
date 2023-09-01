###Using the breeding and adult biomass of birds, can we look at how this relates to algal diversity?


seabirds <- read.csv("../output/n15_seabirds_combined_no_iti.csv", strip.white = T, header = T)
microbes <- read.csv("../output/microbes_summary_bysite.csv", strip.white = T, header = T)
algae <- read.csv("../../algae_diversity_surveys/output/algae_div_summary_bysite.csv", strip.white = T, header = T)

merge.data <- merge(microbes, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
merge.data <- merge(merge.data, algae, by = "site.name", all = TRUE, no.dups = TRUE)
colnames(merge.data) 

merge.data.comp <- merge.data[, -c(6:12,14,20:26,35,36, 38:44,46:50)]
colnames(merge.data.comp)
#Reorder names 
col_order <- c("site.name", "breeding_biomass_kgha_side", "adult_biomass_kgha_side", "N.15_at_10m", "N.15_at_20m", "N.15_at_30m",
               "N.15_at_40m", "mean.algae.richness", "mean.algae.evenness", "mean.algae.shannon", "mean.algae.turbinaria",
               "mean.algae.lobophora", "mean.algae.halimeda", "mean.algae.chl.fastigita", "mean.algae.caul.serrulata",
               "coral.mean.richness", "coral.mean.evenness", "coral.mean.shannon", "coral.mean.faithPD", "coral.mean.betadisp",
               "coral.mean.endozoicomonas", "water.mean.richness", "water.mean.evenness", "water.mean.shannon", "water.mean.faithPD",
               "water.mean.betadisp", "water.mean.litoricola", "water.mean.Synechococcus", "water.mean.Prochlorococcus")  
merge.data.comp <- merge.data.comp[, col_order]
colnames(merge.data.comp)


##Create a matrix
data.matrix <- as.data.frame(merge.data.comp[,2:29]) #Use only the numerical values



#Run a correlation test using the library corrplot
library(corrplot)

cor.mtest(data.matrix)
correlation.matrix <- cor(data.matrix, use = "pairwise.complete.obs")
write.csv(correlation.matrix, "../output/seabirds_algae_microbes_corrmatrix_no_iti.csv")

pdf(file = "../output/seabirds_microbies_algae_no_iti.pdf")

corrplot(cor(data.matrix, use = "pairwise.complete.obs"), type = "upper",
         addCoef.col = NULL, addCoefasPercent = FALSE, tl.col = "black", tl.cex = 0.5, title = "Seabirds Algae Microbes")

dev.off()

