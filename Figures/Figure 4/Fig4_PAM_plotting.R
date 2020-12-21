##Statistical analyses for PAM at T1
library(lmerTest)
library(emmeans)
library(sjPlot)
library(drc)
library(ggplot2)
library(Rmisc)

#need to set wd to source file location

#Load file
Eilat_2019<-read.csv("Eilat_Full_responses_data_clean.csv")
View(Eilat_2019)
Eilat_2019$Geno<- as.factor(Eilat_2019$Geno)
str(Eilat_2019)

#Subset data
Control_check1<-subset(Eilat_2019, Temp=="Field" | Temp=="Start" | Temp=="27")
Control_check<-subset(Control_check1,Timepoint=="Hold")

Eilat_no_field<-subset(Eilat_2019, Temp=="27" | Temp=="29.5" | Temp=="32" | Temp=="34.5")
Eilat_no_field$Tank<- as.factor(Eilat_no_field$Tank)

Hold_all<-subset(Eilat_no_field, Timepoint=="Hold")

#Stats model
PAM_hold<-lmer(PAM ~ Experiment*Temp + (1|Tank) + (1|Geno), data=Hold_all)
sjPlot::plot_model(PAM_hold, type="diag")
step(PAM_hold, reduce.random=FALSE)

PAM_hold_final<-lmer(PAM ~ Temp + (1|Tank) + (1|Geno),data=Hold_all)
anova(PAM_hold_final)

print(emmeans(PAM_hold_final, list(pairwise ~ Temp)), adjust = c("tukey"))

############################################################################################################
############################################################################################################
############################################################################################################

###### Compute and plot T1 Fv/Fm ED50 for acute (CBASS) AND chronic (RSS) experiments

### Temp as numeric for drc fitting
Hold_all$Temp<-as.character(Hold_all$Temp)
Hold_all$Temp<-as.numeric(Hold_all$Temp)

str(Hold_all)

#### Compare curves and ED50s ####

Eilat_DRC <- drm(PAM ~ Temp, data = Hold_all, curveid = Experiment, fct = LL.3(names = c('hill', 'max', 'ed50')))
mselect(Eilat_DRC, list(LL.2(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.4(), LL2.5(), AR.2(), AR.3(), EXD.2(), EXD.3()), icfct = AIC)
modelFit(Eilat_DRC)
summary(Eilat_DRC)
plot(Eilat_DRC)

#extract ED50
Eilat_RSS_coeff<-Eilat_DRC$coefficients[5]
Eilat_CBASS_coeff<-Eilat_DRC$coefficients[6]

#### Run individually for plotting ####

#### CBASS ####

#Run model
Eilat_CBASS <- drm(PAM ~ Temp, data = Hold_all[Hold_all$Experiment=="CBASS",], fct = LL.3())
summary(Eilat_CBASS)
plot(Eilat_CBASS)

#### RSS ####

#Run model
Eilat_RSS <- drm(PAM ~ Temp, data = Hold_all[Hold_all$Experiment=="RSS",],fct = LL.3())
summary(Eilat_RSS)
plot(Eilat_RSS)

############################################################################################
#### Combine ED50 data plus predict curves from models for plotting ####
###########################################################################################

Eilat_coeffs<-data.frame(Eilat_CBASS_coeff, Eilat_RSS_coeff)

Eilat_CBASS_preddata = data.frame(temp = seq(27,39, length.out = 100))
Eilat_CBASS_pred = as.data.frame(predict(Eilat_CBASS, newdata = Eilat_CBASS_preddata, interval = 'confidence'))
Eilat_CBASS_preddata = data.frame(Eilat_CBASS_preddata, fvfm = Eilat_CBASS_pred$Prediction, Lower = Eilat_CBASS_pred$Lower, Upper = Eilat_CBASS_pred$Upper)

Eilat_RSS_preddata = data.frame(temp = seq(27,39, length.out = 100))
Eilat_RSS_pred = as.data.frame(predict(Eilat_RSS, newdata = Eilat_RSS_preddata, interval = 'confidence'))
Eilat_RSS_preddata = data.frame(Eilat_RSS_preddata, fvfm = Eilat_RSS_pred$Prediction, Lower = Eilat_RSS_pred$Lower, Upper = Eilat_RSS_pred$Upper)

#### PLOT  ####
levels(Hold_all$Experiment)

Eilat_DRC_plot<- ggplot() +
  geom_jitter(data = Hold_all, aes(x = Temp, y = PAM, color = Experiment), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(26,40), breaks=c(26,28,30,32,34,36,38)) +
  scale_y_continuous(limits=c(-0.1, 0.7), breaks=c(0, 0.2, 0.4, 0.6)) +
  
  geom_line(data = Eilat_CBASS_preddata, aes(x = temp, y = fvfm), color = '#66ccfe', show.legend = FALSE) +
  geom_ribbon(data = Eilat_CBASS_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#66ccfe', linetype=2, alpha = 0.2) +
  geom_vline(data = Eilat_coeffs, aes(xintercept = Eilat_CBASS_coeff), color = '#66ccfe', show.legend = FALSE) +
  geom_text(data = Eilat_coeffs, aes(label = Eilat_CBASS_coeff), x = 30, y = 0.2, show.legend = FALSE, color = '#66ccfe') +
 
  geom_line(data = Eilat_RSS_preddata, aes(x = temp, y = fvfm), color = '#e92000', show.legend = FALSE) +
  geom_ribbon(data = Eilat_RSS_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#e92000', linetype=2, alpha = 0.2) +
  geom_vline(data = Eilat_coeffs, aes(xintercept = Eilat_RSS_coeff), color = '#e92000', show.legend = FALSE) +
  geom_text(data = Eilat_coeffs, aes(label = Eilat_RSS_coeff), x = 30, y = 0.1, show.legend = FALSE, color = '#e92000') +

  scale_color_manual(values=c('#66ccfe','#e92000')) +
  ylab("Fv/Fm") +
  xlab("Temperature (°C)") +
  theme_bw()

Eilat_DRC_plot

ggsave(Eilat_DRC_plot, height = 6 , width = 10, filename = "Eilat_CBASS_vs_RSS_DRCs.pdf", useDingbats=FALSE)

