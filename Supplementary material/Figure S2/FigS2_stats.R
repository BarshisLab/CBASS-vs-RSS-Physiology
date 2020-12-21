#Statistical analyses for univariate responses
library(lmerTest)
library(emmeans)
library(sjPlot)

#Load file
Eilat_2019<-read.csv("Eilat_Full_responses_data_clean.csv")
View(Eilat_2019)
Eilat_2019$Geno<- as.factor(Eilat_2019$Geno)
Eilat_2019$Tank<- as.factor(Eilat_2019$Tank)
str(Eilat_2019)

#Subset data
Control_check1<-subset(Eilat_2019, Temp=="Field" | Temp=="Start" | Temp=="27")
Control_check<-subset(Control_check1,Timepoint=="Hold")

Eilat_no_field<-subset(Eilat_2019, Temp=="27" | Temp=="29.5" | Temp=="32" | Temp=="34.5")

Hold_all<-subset(Eilat_no_field, Timepoint=="Hold")
Recovery_all<-subset(Eilat_no_field, Timepoint=="Recovery")

#### GLMMs for Chl a per cell

#Field/T0 stats - comparing 27 to field and start of experiment values
Chla.cell_control<-lmer(Chla_cell ~ Temp*Experiment + (1|Tank) + (1|Geno),data=Control_check) #Full model
sjPlot::plot_model(Chla.cell_control, type="diag") # Model diagnostics - checking model assumptions
step(Chla.cell_control,reduce.random=FALSE)
anova(Chla.cell_control)

#T1 - Hold stats
#Protein
Chla.cell_hold<-lmer(Chla_cell ~ Temp*Experiment + (1|Tank) + (1|Geno),data=Hold_all) #Full model
sjPlot::plot_model(Chla.cell_hold, type="diag") # Model diagnostics - checking model assumptions
step(Chla.cell_hold,reduce.random=FALSE) #Model selection

Chla.cell_hold_final<-lmer(Chla_cell ~ Experiment + (1|Geno),data=Hold_all) #Final model based on step output
anova(Chla.cell_hold_final)

############################################################################################################
############################################################################################################
############################################################################################################

#T2 - Recovery time point stats

#Chl a per cell
Chla.cell_recovery<-lmer(Chla_cell ~ Experiment*Temp + (1|Tank) + (1|Geno), data=Recovery_all)
sjPlot::plot_model(Chla.cell_recovery, type="diag")
step(Chla.cell_recovery,reduce.random=FALSE)

Chla.cell_recovery_final<-lmer(Chla_cell ~ Temp + (1 | Geno),data=Recovery_all)
anova(Chla.cell_recovery_final)

print(emmeans(Chla.cell_recovery_final, list(pairwise ~ Temp)), adjust = c("tukey"))
