#Statistical analyses for univariate responses
library(lmerTest)
library(emmeans)
library(sjPlot)
library(car)

#Load file
Eilat_2019<-read.csv("Eilat_Full_responses_data_clean.csv")
View(Eilat_2019)
Eilat_2019$Geno<- as.factor(Eilat_2019$Geno)
Eilat_2019$Tank<- as.factor(Eilat_2019$Tank)
str(Eilat_2019)

#Subset data
Control_check<-subset(Eilat_2019, Temp=="Field" | Temp=="Start" | Temp=="27")
Control_check<-subset(Control_check,Timepoint=="Hold")
Control_check<-subset(Control_check,Experiment=="CBASS")

Eilat_no_field<-subset(Eilat_2019, Temp=="27" | Temp=="29.5" | Temp=="32" | Temp=="34.5")

Hold_all<-subset(Eilat_no_field, Timepoint=="Hold")
Recovery_all<-subset(Eilat_no_field, Timepoint=="Recovery")

#### GLMMs for Respiration and Net photosynthesis by timepoint (T1 and T2)

#Field/T0 stats - comparing 27 to field and start (T0) of experiment values, only for CBASS as no field or T0 values for RSS
#Genotype also not included as only 3 genotypes are present across treatments
Resp_control<-aov(Resp ~ Temp,data=Control_check)
summary(Resp_control)
print(emmeans(Resp_control, list(pairwise ~ Temp)), adjust = c("tukey"))

NetPS_control<-aov(NetPS ~ Temp,data=Control_check)
summary(NetPS_control)
print(emmeans(NetPS_control, list(pairwise ~ Temp)), adjust = c("tukey"))

#T1 - Hold statistical models
#Tank removed due to a lack of replication for respiration and photosynthesis across tanks
#Respiration
Resp_hold<-lmer(Resp ~ Temp*Experiment +(1|Geno),data=Hold_all) #Full model
sjPlot::plot_model(Resp_hold, type="diag") # Model diagnostics - checking model assumptions
step(Resp_hold,reduce.random=FALSE) #Model selection

Resp_hold_final<-lmer(Resp ~ Temp*Experiment +(1|Geno),data=Hold_all)  #Final model based on step output
anova(Resp_hold_final)

print(emmeans(Resp_hold_final, list(pairwise ~ Temp|Experiment)), adjust = c("tukey")) #Post hoc test depending on significant main effects
print(emmeans(Resp_hold_final, list(pairwise ~ Experiment|Temp)), adjust = c("tukey")) 

#Net photosynthesis
NetPS_hold<-lmer(NetPS ~ Experiment*Temp + (1|Geno), data=Hold_all)
sjPlot::plot_model(NetPS_hold, type="diag")
step(NetPS_hold, reduce.random=FALSE)

NetPS_hold_final<-lmer(NetPS ~ Experiment*Temp + (1|Geno),data=Hold_all)
anova(NetPS_hold_final)

print(emmeans(NetPS_hold_final, list(pairwise ~ Experiment|Temp)), adjust = c("tukey")) 
print(emmeans(NetPS_hold_final, list(pairwise ~ Temp|Experiment)), adjust = c("tukey"))

############################################################################################################
############################################################################################################
############################################################################################################

#T2 - Recovery time point stats

#Resp
Resp_recovery<-lmer(Resp ~ Experiment*Temp + (1|Geno), data=Recovery_all)
sjPlot::plot_model(Resp_recovery, type="diag")
step(Resp_recovery, reduce.random=FALSE)

Resp_recovery_final<-lmer(Resp ~ Experiment+Temp + (1 | Geno),data=Recovery_all)
anova(Resp_recovery_final)

print(emmeans(Resp_recovery_final, list(pairwise ~ Temp)), adjust = c("tukey"))

#Net photosynthesis
NetPS_recovery<-lmer(NetPS ~ Experiment*Temp + (1|Geno), data=Recovery_all)
sjPlot::plot_model(NetPS_recovery, type="diag")
step(NetPS_recovery,reduce.random=FALSE)

NetPS_recovery_final<-lmer(NetPS ~ Temp + (1|Geno),data=Recovery_all)
anova(NetPS_recovery_final)

print(emmeans(NetPS_recovery_final, list(pairwise ~ Temp)), adjust = c("tukey"))
