##Microbiome data vs. seabirds and N15
library(plyr)
library(tidyverse)

#Bring in Microbes and seabird data

#seabirds <- read.csv("../output/n15_seabirds_combined_with_iti.csv", strip.white = T, header = T)
seabirds <- read.csv("../output/n15_seabirds_combined_no_iti.csv", strip.white = T, header = T)
microbes <- read.csv("../../microbiome_analyses/downstream_analyses/integration/output/nov2021_microbiome_metrics.csv")

#combine microbes into site averages

microbes <- microbes %>% mutate(site.name = recode(site.name, "A1" = "Aie_Protected", "A2" = 'Aie_Exposed',
                                                   "Re1" = "Reiono_Protected", "Re2" = "Reiono_Exposed", 
                                                   "Rm1" = "Rimatuu_Protected", "Rm2" = "Rimatuu_Exposed"))


microbes.sum <- ddply(microbes, c("sample.type", "site.name"), summarise,
                                   mean.richness = mean(Observed),
                      mean.shannon = mean(Shannon),
                      mean.faithPD = mean(FaithPD),
                      mean.evenness = mean(Evenness),
                      mean.Rhodospirillales_AEGEAN_169 = mean(RelAbund_Rhodospirillales_AEGEAN_169_marine_group),
                      mean.SAR11.CladeIa = mean(RelAbund_SAR11_CladeIa),
                      mean.SAR11.CladeIb = mean(RelAbund_SAR11_CladeIb),
                      mean.SAR11.CladeII = mean(RelAbund_SAR11_CladeII),
                      mean.SAR116 = mean(RelAbund_SAR116_Clade),
                      mean.SAR86 = mean(RelAbund_SAR86_Clade),
                      mean.Cellvibrionales_Aestuariicella = mean(RelAbund_Cellvibrionales_Aestuariicella),
                      mean.Ca.Ameobophilus = mean(RelAbund_Cytophagales_Candidatus_Amoebophilus),
                      mean.Glaciecola = mean(RelAbund_Alteromonadales_Glaciecola),
                      mean.Meptuniibacter = mean(RelAbund_Alteromonadales_Glaciecola),
                      mean.Porticoccus = mean(RelAbund_Cellvibrionales_Porticoccus),
                      mean.Woesearchaeales = mean(RelAbund_Woesearchaeales_SCGC_AAA286E23),
                      mean.Unc.Alteromonadaceae = mean(RelAbund_Uncultured_Alteromonadaceaeae),
                      mean.Endozoicomonas = mean(RelAbund_Endozoicomonas),
                      mean.Litoricola = mean(RelAbund_Litoricola),
                      mean.Synechococcus = mean(RelAbund_Synechococcus),
                      mean.Prochlorococcus = mean(RelAbund_Prochlorococcus),
                      mean.betadisp = mean(beta_dispersion_motu_islandside)
                                   
)

microbes.sum.coral <- microbes.sum %>% filter(sample.type == "coral")
microbes.sum.coral <- microbes.sum.coral %>% rename(c("coral.mean.richness" = "mean.richness", "coral.mean.shannon" = "mean.shannon",
                                                      "coral.mean.faithPD" = "mean.faithPD", "coral.mean.evenness" = "mean.evenness",
                                                      "coral.mean.aestuariicella" = "mean.Cellvibrionales_Aestuariicella", "coral.mean.Ca.amoebophilus" = "mean.Ca.Ameobophilus",
                                                      "coral.mean.glaciecola" = "mean.Glaciecola", "coral.mean.neptuniibacter" = "mean.Meptuniibacter", 
                                                      "coral.mean.porticoccus" = "mean.Porticoccus", "coral.mean.woesearchaeales" = "mean.Woesearchaeales", 
                                                      "coral.mean.alteromonadaceae" = "mean.Unc.Alteromonadaceae",
                                                      "coral.mean.endozoicomonas" = "mean.Endozoicomonas", "coral.mean.litoricola" = "mean.Litoricola",
                                                      "coral.mean.betadisp" = "mean.betadisp"))
microbes.sum.coral <- microbes.sum.coral[, -c(1,7:12, 22, 23)]

microbes.sum.water <- microbes.sum %>% filter(sample.type == "water")
microbes.sum.water <- microbes.sum.water %>% rename(c("water.mean.richness" = "mean.richness", "water.mean.shannon" = "mean.shannon",
                                                      "water.mean.faithPD" = "mean.faithPD", "water.mean.evenness" = "mean.evenness",
                                                      "water.mean.rhodospirillales_aegean169" = "mean.Rhodospirillales_AEGEAN_169", 
                                                      "water.mean.SAR11.cladeIa" = "mean.SAR11.CladeIa",
                                                      "water.mean.SAR11.CladeIb" = "mean.SAR11.CladeIb",
                                                      "water.mean.SAR11.CladeII" = "mean.SAR11.CladeII", "water.mean.SAR116" = "mean.SAR116", "water.mean.SAR86" = "mean.SAR86",
                                                      "water.mean.alteromonadaceae" = "mean.Unc.Alteromonadaceae",
                                                      "water.mean.litoricola" = "mean.Litoricola",
                                                      "water.mean.Synechococcus" = "mean.Synechococcus", "water.mean.Prochlorococcus" = "mean.Prochlorococcus",
                                                      "water.mean.betadisp" = "mean.betadisp"))
microbes.sum.water <- microbes.sum.water[, -c(1,13:18,20)]

microbes.all <- merge(microbes.sum.coral, microbes.sum.water, by = "site.name", all = TRUE, no.dups = TRUE)
View(microbes.all)
write.csv(microbes.all, "../output/microbes_summary_bysite.csv", row.names = FALSE)

##Combine microbes with seabirds from the seabirds_vs_n15 script *** This is a change in the code - double check & come back

microbes.seabirds <- merge(microbes.all, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
colnames(microbes.seabirds) 

##Create a matrix
data.matrix <- as.data.frame(microbes.seabirds[,2:50]) #Use only the numerical values


#Run a correlation test using the library corrplot
library(corrplot)

cor.mtest(data.matrix)
correlation.matrix <- cor(data.matrix, use = "pairwise.complete.obs")
write.csv(correlation.matrix, "../output/seabirds_mmicrobes_corrmatrix_with_iti.csv") #change to "_no_iti" if removing iti

pdf(file = "../output/seabirds_v_microbes_with_iti.pdf")

corrplot(cor(data.matrix, use = "pairwise.complete.obs"), type = "upper",
         addCoef.col = NULL, addCoefasPercent = FALSE, tl.col = "black", tl.cex = 0.5, title = "seabirds vs. Microbes")

dev.off()




