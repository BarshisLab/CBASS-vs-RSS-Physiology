#Statistical analyses for univariate responses
library(lmerTest)
library(emmeans)
library(sjPlot)

#Load file
Eilat_2019<-read.csv("Eilat_Full_responses_data_clean.csv")
View(Eilat_2019)
Eilat_2019$Geno<- as.factor(Eilat_2019$Geno)
Eilat_2019$Tank<- as.factor(Eilat_2019$Tank)
Eilat_2019$Temp<- as.factor(Eilat_2019$Temp)

str(Eilat_2019)

#Subset data
Control_check1<-subset(Eilat_2019, Temp=="Field" | Temp=="Start" | Temp=="27")
Control_check<-subset(Control_check1,Timepoint=="Hold")

Eilat_no_field<-subset(Eilat_2019, Temp=="27" | Temp=="29.5" | Temp=="32" | Temp=="34.5")

Hold_all<-subset(Eilat_no_field, Timepoint=="Hold")
Recovery_all<-subset(Eilat_no_field, Timepoint=="Recovery")

#### GLMMs for Protein, Chl a, and Sym counts by timepoint (T1 and T2)

#Field/T0 stats - comparing 27 to field and start of experiment values
Protein_control<-lmer(Protein_mg_cm2 ~ Experiment*Temp + (1|Geno),data=Control_check)
step(Protein_control,reduce.random=FALSE)
anova(Protein_control)

Chl_control<-lmer(Chla_cm2 ~ Experiment*Temp + (1|Geno),data=Control_check)
step(Chl_control)
anova(Chl_control)

Sym_control<-lmer(Sym_density ~ Experiment*Temp + (1|Geno),data=Control_check)
step(Sym_control)
anova(Sym_control)

#T1 - Hold statistical models

#Protein
Protein_hold<-lmer(Protein_mg_cm2 ~ Temp*Experiment + (1|Tank) + (1|Geno),data=Hold_all) #Full model
sjPlot::plot_model(Protein_hold, type="diag") # Model diagnostics - checking model assumptions
step(Protein_hold) #Model selection

Protein_hold_final<-lmer(Protein_mg_cm2 ~ Temp + (1|Tank) + (1|Geno),data=Hold_all) #Final model based on step output
summary(Protein_hold_final)
anova(Protein_hold_final)
residuals(Protein_hold_final)

print(emmeans(Protein_hold_final, list(pairwise ~ Temp)), adjust = c("tukey")) #Post hoc test depending on significant main effects

#Chl a
Chla_hold<-lmer(Chla_cm2 ~ Experiment*Temp + (1|Tank) + (1|Geno), data=Hold_all)
sjPlot::plot_model(Chla_hold, type="diag")
step(Chla_hold)

Chla_hold_final<-lmer(Chla_cm2 ~ Temp + (1|Tank) + (1|Geno),data=Hold_all)
anova(Chla_hold_final)

print(emmeans(Chla_hold_final, list(pairwise ~ Temp)), adjust = c("tukey"))
print(emmeans(Chla_hold, list(pairwise ~ Experiment | Temp)), adjust = c("tukey"))

#Symbiont density
Sym_log<-log(Hold_all$Sym_density)
Sym_hold<-lmer(Sym_log ~ Experiment*Temp + (1|Tank) + (1|Geno), data=Hold_all)
sjPlot::plot_model(Sym_hold, type="diag")
step(Sym_hold)

Sym_hold_final<-lmer(Sym_log ~ Experiment*Temp + (1|Geno),data=Hold_all)
anova(Sym_hold_final)

print(emmeans(Sym_hold_final, list(pairwise ~ Temp|Experiment)), adjust = c("tukey"))
print(emmeans(Sym_hold_final, list(pairwise ~ Experiment|Temp)), adjust = c("tukey"))

############################################################################################################
############################################################################################################
############################################################################################################

#T2 - Recovery time point stats

#Protein
Protein_recovery<-lmer(Protein_mg_cm2 ~ Experiment*Temp + (1|Tank) + (1|Geno), data=Recovery_all)
sjPlot::plot_model(Protein_recovery, type="diag")
step(Protein_recovery)

Protein_recovery_final<-lmer(Protein_mg_cm2 ~ Temp + (1 | Geno),data=Recovery_all)
anova(Protein_recovery_final)

print(emmeans(Protein_recovery_final, list(pairwise ~ Temp)), adjust = c("tukey"))

#Chl a
Chla_recovery<-lmer(Chla_cm2 ~ Experiment*Temp + (1|Tank) + (1|Geno), data=Recovery_all)
sjPlot::plot_model(Chla_recovery, type="diag")
step(Chla_recovery)

Chla_recovery_final<-lmer(Chla_cm2 ~ Experiment*Temp + (1 | Geno),data=Recovery_all)
anova(Chla_recovery_final)

print(emmeans(Chla_recovery_final, list(pairwise ~ Temp|Experiment)), adjust = c("tukey"))
print(emmeans(Chla_recovery_final, list(pairwise ~ Experiment|Temp)), adjust = c("tukey"))

#Symbiont density
Sym_recovery<-lmer(Sym_density ~ Experiment*Temp + (1|Tank) + (1|Geno), data=Recovery_all)
sjPlot::plot_model(Sym_recovery, type="diag")
step(Sym_recovery)

Sym_recovery_final<-lmer(Sym_density ~ Experiment*Temp + (1 | Geno),data=Recovery_all)
anova(Sym_recovery_final)

print(emmeans(Sym_recovery_final, list(pairwise ~ Temp|Experiment)), adjust = c("tukey"))
print(emmeans(Sym_recovery_final, list(pairwise ~ Experiment|Temp)), adjust = c("tukey"))

#End of script