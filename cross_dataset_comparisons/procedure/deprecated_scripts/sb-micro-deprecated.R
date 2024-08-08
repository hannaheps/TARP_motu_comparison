##Deprecated microbe vs. seabird code
#Saving JUST IN CASE


##Observed Richness of coral mirobiomes
mod.glmm1 <- glmer(Observed~seabird_level*Distance_to_shore + Exposure + (1|site.name), family = "poisson", data = coral.microbes.seabirds)
summary(mod.glmm1)
anova(mod.glmm1)
Anova(mod.glmm1)
car::Anova(mod.glmm1)
#Analysis of Variance Table
#                                npar Sum Sq Mean Sq F value
#seabird_level                      2  20.90  10.451 10.4514
#Distance_to_shore                  3 103.38  34.460 34.4603
#Exposure                           1   6.36   6.363  6.3635
#seabird_level:Distance_to_shore    6 544.10  90.683 90.6829


#Analysis of Deviance Table (Type II Wald chisquare tests)

#Response: Observed
#                                   Chisq Df Pr(>Chisq)    
#seabird_level                    21.6564  2  1.983e-05 ***
#Distance_to_shore               103.0863  3  < 2.2e-16 ***
#Exposure                          6.5432  1    0.01053 *  
#seabird_level:Distance_to_shore 544.4811  6  < 2.2e-16 ***

#Fixed effects:
#                                      Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                            4.13043    0.27251  15.157  < 2e-16 ***
#seabird_levelmid                       0.24318    0.33346   0.729  0.46583    
#seabird_levelhigh                     -0.14842    0.36376  -0.408  0.68327    
#Distance_to_shore20                   -0.40370    0.10162  -3.973 7.11e-05 ***
#Distance_to_shore30                   -0.16461    0.09202  -1.789  0.07364 .  
#Distance_to_shore40                   -0.69699    0.11931  -5.842 5.17e-09 ***
#ExposureProtected                     -0.97991    0.38308  -2.558  0.01053 *  
#seabird_levelmid:Distance_to_shore20   0.37837    0.11925   3.173  0.00151 ** 
#seabird_levelhigh:Distance_to_shore20  2.33251    0.17799  13.104  < 2e-16 ***
#seabird_levelmid:Distance_to_shore30   0.56877    0.10992   5.175 2.28e-07 ***
#seabird_levelhigh:Distance_to_shore30  0.55734    0.18482   3.016  0.00256 ** 
#seabird_levelmid:Distance_to_shore40   0.82496    0.13360   6.175 6.62e-10 ***
#seabird_levelhigh:Distance_to_shore40  1.81011    0.19259   9.399  < 2e-16 ***

mod.glmm.r <- glmer(Observed~Distance_to_shore + Exposure + (1|site.name), family = "poisson", data = coral.microbes.seabirds)
anova(mod.glmm1, mod.glmm.r)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#mod.glmm.r    6 3841.3 3854.4 -1914.6   3829.3                         
#mod.glmm1    14 3234.3 3265.0 -1603.2   3206.3 622.97  8  < 2.2e-16 ***

mod.glmm2 <- glmer(Observed~algae.N15*Distance_to_shore + Exposure  + (1|site.name), data=coral.microbes.seabirds, family="poisson")

#Fixed effects:
#                               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    3.746877   0.336037  11.150  < 2e-16 ***
#algae.N15                      0.067630   0.028947   2.336   0.0195 *  
#Distance_to_shore20           -3.116892   0.277343 -11.238  < 2e-16 ***
#Distance_to_shore30            0.218754   0.308577   0.709   0.4784    
#Distance_to_shore40            1.807384   0.273079   6.619 3.63e-11 ***
#ExposureProtected             -0.714900   0.368326  -1.941   0.0523 .  
#algae.N15:Distance_to_shore20  0.515105   0.038044  13.540  < 2e-16 ***
#algae.N15:Distance_to_shore30 -0.001467   0.050233  -0.029   0.9767    
#algae.N15:Distance_to_shore40 -0.312308   0.046022  -6.786 1.15e-11 ***

summary(mod.glmm2)
Anova(mod.glmm2)
#Analysis of Deviance Table (Type II tests)#

#Response: Observed
#                              Chisq Df Pr(>Chisq)    
#algae.N15                    89.8691  1    < 2e-16 ***
#Distance_to_shore           172.7879  3    < 2e-16 ***
#Exposure                      3.7673  1    0.05227 .  
#algae.N15:Distance_to_shore 308.0568  3    < 2e-16 ***


#Okay these are significant????!!!

emmeans(mod.glmm1, list(pairwise ~ seabird_level*Distance_to_shore), adjust = "fdr")


emmeans(mod.glmm2, list(pairwise~ algae.N15*Distance_to_shore), adjust = "fdr")



mod.glmm1 %>% 
  emmeans(~ seabird_level,
          type = "response") 

emmeans(mod.glmm1, list(pairwise ~ seabird_level*Distance_to_shore), type = "response", adjust = "fdr")


plot_data_richness_sb<-emmip(mod.glmm1, ~ seabird_level*Distance_to_shore,
                             type = "response", CIs = TRUE, plotit=FALSE)

plot_data_richness_sb

response_plot_richness<-
  ggplot(data = plot_data_richness_sb, aes(x = Distance_to_shore, y = yvar, colour = seabird_level, fill = seabird_level)) +#, group = seabird_level))+
  #ggplot(data = plot_data_richness_sb, aes(x = seabird_level, y = yvar, colour = seabird_level, fill = seabird_level))+
  geom_pointrange(aes(ymin = UCL, ymax = LCL), alpha = .3, linewidth = 2, size = 0.2, position = position_dodge(.1)) +
  geom_line(alpha = .5 ,  position = position_dodge(.1)) +
  scale_fill_manual(values = c( "low" = "#CD1913", "mid" = "#F2BB05", "high" ="#2F9D3E")) + 
  scale_colour_manual(values = c( "low" = "#CD1913", "mid" = "#F2BB05", "high" = "#2F9D3E")) +
  ylab("Coral Microbiome Observed ASV Richness") +
  xlab("Distance From Shore (m)")+
  labs(color='Seabird Level', fill = 'Seabird Level') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        legend.position = c(.85,.85))

response_plot_richness
ggsave(response_plot_richness, file = "../output/seabird-microbes/coral_micro_richness_model.jpg", width = 7, height = 5)



##N15
mod.glmm.n15 <-glmer(Observed~algae.N15*Distance_to_shore + Exposure + (1|site.name), family = "poisson", data = coral.microbes.seabirds)
summary(mod.glmm.n15)
#Fixed effects:
#                               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    3.746877   0.336037  11.150  < 2e-16 ***
#algae.N15                      0.067630   0.028947   2.336   0.0195 *  
#Distance_to_shore20           -3.116892   0.277343 -11.238  < 2e-16 ***
#Distance_to_shore30            0.218754   0.308577   0.709   0.4784    
#Distance_to_shore40            1.807384   0.273079   6.619 3.63e-11 ***
#ExposureProtected             -0.714900   0.368326  -1.941   0.0523 .  
#algae.N15:Distance_to_shore20  0.515105   0.038044  13.540  < 2e-16 ***
#algae.N15:Distance_to_shore30 -0.001467   0.050233  -0.029   0.9767    
#algae.N15:Distance_to_shore40 -0.312308   0.046022  -6.786 1.15e-11 ***


anova(mod.glmm.n15)
#                            npar  Sum Sq Mean Sq  F value
#algae.N15                      1  93.016  93.016  93.0164
#Distance_to_shore              3 171.078  57.026  57.0261
#Exposure                       1   1.012   1.012   1.0123
#algae.N15:Distance_to_shore    3 309.200 103.067 103.0668

effect_plot(mod.glmm.n15, pred = algae.N15, interval = TRUE, plot.points = TRUE, data = coral.microbes.seabirds)

gauss.pd.lm<-lmerTest::lmer(Evenness~ seabird.level*Distance_to_shore + Exposure + (1|site.name), data=coral.microbes.seabirds)
summary(gauss.pd.lm)
Anova(gauss.pd.lm)

anova(gauss.pd.lm)

library(aod)
?wald.test()
wald.test(b=fixef(mod.glmm),
          Sigma=vcov(mod.glmm), Terms = 2:7 )
intervals(mod.glmm)
library(gmodels)
ci(mod.glmm)
#No significant parameters. 

#What about a quasi-R2?
totalss <- var(resid(mod.glmm,type='pearson')+predict(mod.glmm, type='link'))
1-var(residuals(mod.glmm, type='pearson'))/(totalss)
#0.2502419

#type = response emmeans  for beta as long as family is specified in model 

mod.glmm2 <- glmmPQL(Observed~algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family="poisson")
summary(mod.glmm2)
Anova(mod.glmm2)















