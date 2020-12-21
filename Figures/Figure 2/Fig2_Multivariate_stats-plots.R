## ###################################### ##
## ###################################### ##
#### Multivariate stats and plotting ####
## ##################################### ##
## ###################################### ##

#### Load and view file
Eilat_2019<- read.csv("Eilat_Full_responses_data_clean.csv")
View(Eilat_2019)
str(Eilat_2019)

#Required packages
library(vegan)
library(dplyr)

#Define data
MVdata<-transform(Eilat_2019, 
                  Experiment=as.factor(Experiment), 
                  Timepoint=as.factor(Timepoint),
                  Temp=as.factor(Temp),Tank=as.factor(Tank))
str(MVdata)
names(MVdata)

#Select columns of interest + Remove any rows that have NAs to be able to compute complete distance matrix
MVdata_T1<-select(MVdata,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS,FvFm,Experiment,Geno,Timepoint,Temp,Tank)
MVdata_T1_complete <- na.omit(MVdata_T1)

MVdata_T2<-select(MVdata,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS,Experiment,Geno,Timepoint,Temp,Tank)
MVdata_T2_complete <- na.omit(MVdata_T2)

#Standardize the data
Standard_T1_complete<-MVdata_T1_complete %>% mutate_if(is.numeric, scale)
Standard_T2_complete<-MVdata_T2_complete %>% mutate_if(is.numeric, scale)

#Final step to separate T1 and T2 data for analysis and plotting
T1_complete<-subset(Standard_T1_complete, Timepoint=="Hold")
T2_complete<-subset(Standard_T2_complete, Timepoint=="Recovery")

# Creating separate response (com) and predictor (meta) variable files and compute distance matrices for each
# T1 full
T1_com<-select(Standard_T1_complete,FvFm,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS)
T1_meta<-select(Standard_T1_complete,Experiment,Temp,Tank)
str(T1_com)
str(T1_meta)

T1_dist<- vegdist(T1_com, method="euclidean",na.rm = TRUE)

# PERMANOVA T1
Perm_T1<-adonis(T1_dist~Experiment*Temp, data=T1_meta)
Perm_T1

# T2 full
T2_com<-select(T2_complete,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS)
T2_meta<-select(T2_complete,Experiment,Temp,Tank)
str(T2_com)
str(T2_meta)

T2_dist<- vegdist(T2_com, method="euclidean",na.rm = TRUE)

# PERMANOVA T2
Perm_T2<-adonis(T2_dist~Experiment*Temp, data=T2_meta)
Perm_T2

#Post hoc function - CBASS and RSS need to be separated within time points to compute pairwise comparisons 
#                    (hence the subsetting options above)

##start copy here for function pairwise.adonis() 
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='holm')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 
## end copy here

pairwise.adonis(T1_com,factors=T1_meta$Temp,sim.function='vegdist',sim.method='euclidean',p.adjust.m='bonferroni')

####Extra subsetting to run PERMANOVA, and thus be able to run full pairwise comparisons

#T1 field control
Field_control<-subset(Standard_T1_complete, Experiment =="CBASS" & Temp=="27" | Temp=="Start" | Temp=="Field")

Field_com<-select(Field_control,FvFm,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS)
Field_meta<-select(Field_control,Temp)

Field_dist<- vegdist(Field_com, method="euclidean",na.rm = TRUE)
Field_Perm<-adonis(Field_dist~Temp, data=Field_meta)
Field_Perm

Field_post<-pairwise.adonis(Field_com,factors=Field_meta$Temp,sim.function='vegdist',sim.method='euclidean',p.adjust.m='bonferroni')
Field_post

# T1 - CBASS
T1_CBASS<-subset(T1_complete, Experiment=="CBASS" & Temp=="27" | Temp=="29.5" | Temp=="32" | Temp=="34.5") #without field and T0 (they're compared to 27°C treatment above)

T1_CBASS_com<-select(T1_CBASS,FvFm,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS)
T1_CBASS_meta<-select(T1_CBASS,Temp)

T1_CBASS_dist<- vegdist(T1_CBASS_com, method="euclidean",na.rm = TRUE)
T1_CBASS_Perm<-adonis(T1_CBASS_dist~Temp, data=T1_CBASS_meta)
T1_CBASS_Perm

T1_CBASS_post<-pairwise.adonis(T1_CBASS_com,factors=T1_CBASS_meta$Temp,sim.function='vegdist',sim.method='euclidean',p.adjust.m='bonferroni')
T1_CBASS_post

# T1 - RSS
T1_RSS<-subset(T1_complete, Experiment=="RSS")

T1_RSS_com<-select(T1_RSS,FvFm,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS)
T1_RSS_meta<-select(T1_RSS,Temp)

T1_RSS_dist<- vegdist(T1_RSS_com, method="euclidean",na.rm = TRUE)
T1_RSS_Perm<-adonis(T1_RSS_dist~Temp, data=T1_RSS_meta)
T1_RSS_Perm

T1_RSS_post<-pairwise.adonis(T1_RSS_com,factors=T1_RSS_meta$Temp,sim.function='vegdist',sim.method='euclidean',p.adjust.m='bonferroni')
T1_RSS_post

# T2 - CBASS
T2_CBASS<-subset(T2_complete, Experiment=="CBASS")

T2_CBASS_com<-select(T2_CBASS,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS)
T2_CBASS_meta<-select(T2_CBASS,Temp)

T2_CBASS_dist<- vegdist(T2_CBASS_com, method="euclidean",na.rm = TRUE)
T2_CBASS_Perm<-adonis(T2_CBASS_dist~Temp, data=T2_CBASS_meta)
T2_CBASS_Perm

T2_CBASS_post<-pairwise.adonis(T2_CBASS_com,factors=T2_CBASS_meta$Temp,sim.function='vegdist',sim.method='euclidean',p.adjust.m='bonferroni')
T2_CBASS_post

# T2 - RSS
T2_RSS<-subset(T2_complete, Experiment=="RSS")

T2_RSS_com<-select(T2_RSS,Sym_density,Protein_mg_cm2,Chla_cm2,Chla_cell,Resp,NetPS)
T2_RSS_meta<-select(T2_RSS,Temp)

T2_RSS_dist<- vegdist(T2_RSS_com, method="euclidean",na.rm = TRUE)
T2_RSS_Perm<-adonis(T2_RSS_dist~Temp, data=T2_RSS_meta)
T2_RSS_Perm

T2_RSS_post<-pairwise.adonis(T2_RSS_com,factors=T2_RSS_meta$Temp,sim.function='vegdist',sim.method='euclidean',p.adjust.m='bonferroni')
T2_RSS_post

##### PREPARING SCORES FOR PLOTTING

#Extract NMDS scores and separate data scores by time point for plotting
Sol1 <- metaMDS(T1_dist)
Sol1
stressplot(Sol1)

t1.data.scores <- as.data.frame(scores(Sol1))
t1.data.scores <- cbind(t1.data.scores,T1_meta)

###
Sol2 <- metaMDS(T2_dist)
Sol2
stressplot(Sol2)

t2.data.scores <- as.data.frame(scores(Sol2))
t2.data.scores <- cbind(t2.data.scores,T2_meta)

#create t1 treatment column to plot together
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="CBASS") & (t1.data.scores$Temp=="Field")]<-"Field-CBASS"
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="CBASS") & (t1.data.scores$Temp=="Start")]<-"T0-CBASS"
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="CBASS") & (t1.data.scores$Temp=="27")]<-"27-CBASS"
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="CBASS") & (t1.data.scores$Temp=="29.5")]<-"29.5-CBASS"
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="CBASS") & (t1.data.scores$Temp=="32")]<-"32-CBASS"
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="CBASS") & (t1.data.scores$Temp=="34.5")]<-"34.5-CBASS"

t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="RSS") & (t1.data.scores$Temp=="27")]<-"27-RSS"
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="RSS") & (t1.data.scores$Temp=="29.5")]<-"29.5-RSS"
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="RSS") & (t1.data.scores$Temp=="32")]<-"32-RSS"
t1.data.scores$treatment[Name=(t1.data.scores$Experiment=="RSS") & (t1.data.scores$Temp=="34.5")]<-"34.5-RSS"

#create t2 treatment column to plot together
t2.data.scores$treatment[Name=(t2.data.scores$Experiment=="CBASS") & (t2.data.scores$Temp=="27")]<-"27-CBASS"
t2.data.scores$treatment[Name=(t2.data.scores$Experiment=="CBASS") & (t2.data.scores$Temp=="29.5")]<-"29.5-CBASS"
t2.data.scores$treatment[Name=(t2.data.scores$Experiment=="CBASS") & (t2.data.scores$Temp=="32")]<-"32-CBASS"
t2.data.scores$treatment[Name=(t2.data.scores$Experiment=="CBASS") & (t2.data.scores$Temp=="34.5")]<-"34.5-CBASS"

t2.data.scores$treatment[Name=(t2.data.scores$Experiment=="RSS") & (t2.data.scores$Temp=="27")]<-"27-RSS"
t2.data.scores$treatment[Name=(t2.data.scores$Experiment=="RSS") & (t2.data.scores$Temp=="29.5")]<-"29.5-RSS"
t2.data.scores$treatment[Name=(t2.data.scores$Experiment=="RSS") & (t2.data.scores$Temp=="32")]<-"32-RSS"
t2.data.scores$treatment[Name=(t2.data.scores$Experiment=="RSS") & (t2.data.scores$Temp=="34.5")]<-"34.5-RSS"

#### PLOTTING
library(ggplot2)

#colour assignments for T1
cols <- c("27-RSS" = "navyblue", "29.5-RSS" = "yellow1", "32-RSS" = "darkorange2", "34.5-RSS" = "red3",
          "Field-CBASS" = "black", "T0-CBASS" = "royalblue1","27-CBASS" = "navyblue", "29.5-CBASS" = "yellow1", "32-CBASS" = "darkorange2", "34.5-CBASS" = "red3")

#T1 ggplot - added colours to scale fill manual for field and T0 + values to scale shape manual
ggplot(t1.data.scores, aes(x = NMDS1, y = NMDS2, col=treatment, shape=treatment)) + scale_colour_manual(values=cols) +
  stat_ellipse(geom = "polygon", aes(fill=treatment), alpha = 0.50, size= 0.5, level = 0.66) + 
  scale_fill_manual(values=c("transparent","navyblue","transparent","yellow1","transparent","darkorange2","transparent","red3","transparent", "transparent")) +  # ellipses encompassing 2/3 of replicates
  geom_point(size = 2,stroke = 0.8, aes(shape=treatment)) +  scale_shape_manual(values=c(21,16, 21, 16, 21, 16, 21, 16,21,21)) +                                     # sample scores
  coord_fixed() +                                              # same axis scaling
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  geom_vline(xintercept=0, linetype="dashed") + geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  scale_y_continuous(limits=c(-4.5,4.5), breaks = c(-4,-2,0,2,4)) + scale_x_continuous(limits=c(-4.6,4.5), breaks = c(-4,-2,0,2,4))

#### T2 plot
#colour assignments for T2
cols <- c("27-RSS" = "navyblue", "29.5-RSS" = "yellow1", "32-RSS" = "darkorange2", "34.5-RSS" = "red3",
          "27-CBASS" = "navyblue", "29.5-CBASS" = "yellow1", "32-CBASS" = "darkorange2", "34.5-CBASS" = "red3")

#T2 ggplot (fewer shapes/colours without field and T0)
ggplot(t2.data.scores, aes(x = NMDS1, y = NMDS2, col=treatment, shape=treatment)) + scale_colour_manual(values=cols) +
  stat_ellipse(geom = "polygon", aes(fill=treatment), alpha = 0.50, size= 0.5, level = 0.66) + 
  scale_fill_manual(values=c("transparent","navyblue","transparent","yellow1","transparent","darkorange2","transparent","red3")) +
  geom_point(size = 2, aes(shape=treatment)) +  scale_shape_manual(values=c(21,16, 21, 16, 21, 16, 21, 16)) +  # sample scores
  coord_fixed() + # same axis scaling
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  geom_vline(xintercept=0, linetype="dashed") + geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  scale_y_continuous(limits=c(-4.5,5), breaks = c(-4,-2,0,2,4)) + scale_x_continuous(limits=c(-4.6,4.5), breaks = c(-4,-2,0,2,4))

#Plotting vectors for T1
ord_IR<-metaMDS(T1_dist,distance = "euclidean",
                k=2,trymax=1000,autotransform=TRUE,expand=FALSE, plot=FALSE)
ord_IR
fit<-envfit(ord_IR,T1_com, perm=999)
fit
scores(fit, "vectors")
plot(ord_IR,col="white",axes=FALSE, ylim=c(-4.5,4.5), xlim=c(-4.6,4.5))
plot(fit, p.max = 0.95, col = "black") + abline(h = 0, v = 0, col = "black", lty = 2)
axis(side=1, at=seq(-4, 2, by=2))
axis(side=2, at=seq(-2, 4, by=2))
box()

#Plotting vectors for T2
ord_IR2<-metaMDS(T2_dist,distance = "euclidean",
                k=2,trymax=1000,autotransform=TRUE,expand=FALSE, plot=FALSE)
ord_IR2
fit2<-envfit(ord_IR2,T2_com, perm=999)
fit2
scores(fit2, "vectors")
plot(ord_IR2,col="white",axes=FALSE, ylim=c(-4.5,4.5), xlim=c(-4.6,4.5))
plot(fit2, p.max = 0.95, col = "black") + abline(h = 0, v = 0, col = "black", lty = 2)
axis(side=1, at=seq(-4, 2, by=2))
axis(side=2, at=seq(-2, 4, by=2))
box()

### Plots for each time point (along with vectors) were then combined and cleaned up using Illustrator.
