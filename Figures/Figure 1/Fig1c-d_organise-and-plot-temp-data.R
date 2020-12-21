############################################# 
#### Sorting and plotting RSS temp. data #### 
#############################################

#bind and merge data
RSS1<-read.csv("RSS1.csv")
RSS1$Date<-as.POSIXct(RSS1$Date, format="%Y-%m-%d %H:%M:%S")
RSS2<-read.csv("RSS2.csv")
RSS2$Date<-as.POSIXct(RSS2$Date, format="%Y-%m-%d %H:%M:%S")
RSS_main<-rbind(RSS1, RSS2)

library(openair)
RSS_main$Date<-as.POSIXct(RSS_main$Date)
RSS_main$Tank<-as.factor(RSS_main$Tank)
names(RSS_main)[names(RSS_main) == "Date"] <- "date"
levels(RSS_main$Tank)[levels(RSS_main$Tank)=="AQUA_125_TEMP_PV"] <- "34.5A"
levels(RSS_main$Tank)[levels(RSS_main$Tank)=="AQUA_126_TEMP_PV"] <- "34.5B"
levels(RSS_main$Tank)[levels(RSS_main$Tank)=="AQUA_127_TEMP_PV"] <- "32A"
levels(RSS_main$Tank)[levels(RSS_main$Tank)=="AQUA_128_TEMP_PV"] <- "32B"

levels(RSS_main$Tank)[levels(RSS_main$Tank)=="AQUA_225_TEMPERATURE_2"] <- "29.5A"
levels(RSS_main$Tank)[levels(RSS_main$Tank)=="AQUA_226_TEMPERATURE_2"] <- "29.5B"
levels(RSS_main$Tank)[levels(RSS_main$Tank)=="AQUA_227_TEMPERATURE_2"] <- "27A"
levels(RSS_main$Tank)[levels(RSS_main$Tank)=="AQUA_228_TEMPERATURE_2"] <- "27B"

#Convert to hourly average temperatures for plotting
RSS_av<-timeAverage(RSS_main, avg.time = "hour", type = "Tank")

#select dates of the experiment within the dataframe and remove NAs
RSS_cut<-selectByDate(RSS_av, start = "2019-01-20", end = "2019-02-03")
RSS_cut <- na.omit(RSS_cut)

#Calculate means during hold - separate calculation for 34.5°C tanks due to shorter hold period
RSS_hold<-selectByDate(RSS_main, start = "2019-01-24", end = "2019-01-29")

library(plyr)
RSS_summary_data <- ddply(RSS_hold, c("Tank"), summarise, N= length(Temp), mean = mean(Temp),
                      sd   = sd(Temp),
                      se   = sd / sqrt(N))

RSS_hold_34<-selectByDate(RSS_main, start = "2019-01-24", end = "2019-01-25")
RSS_34_summary_data <- ddply(RSS_hold_34, c("Tank"), summarise, N= length(Temp), mean = mean(Temp),
                          sd   = sd(Temp),
                          se   = sd / sqrt(N))

#plotting Fig. 1d
library(ggplot2)
cols <- c("27A" = "navyblue", "29.5A" = "yellow2", "32A" = "darkorange2","34.5A" = "red3",
          "27B" = "navyblue", "29.5B" = "yellow2", "32B" = "darkorange2", "34.5B" = "red3")
ggplot(RSS_cut, aes(x=date, y=Temp, col=Tank)) + scale_colour_manual(values=cols) +
  geom_line(aes(linetype=Tank)) +
  scale_linetype_manual(values=c("solid", "dashed","solid", "dashed","solid", "dashed","solid","dashed")) +
  scale_y_continuous(limits=c(25.5,36), breaks = c(27,29,31,33,35)) + 
  scale_x_datetime(name = "Date", date_labels = "%b %d", date_breaks = "2 day") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.3,"cm")) +
  geom_hline(yintercept=c(27,29.5,32,34.5), linetype="dashed", color = c("navyblue","yellow2","darkorange2","red3"), size=0.5) +
  geom_vline(xintercept=c(as.numeric(RSS_cut$date[81])), linetype="longdash", color = "black", size=0.5) +
  geom_vline(xintercept=c(as.numeric(RSS_cut$date[1120])), linetype="longdash", color = "black", size=0.5) +
  geom_vline(xintercept=c(as.numeric(RSS_cut$date[1788])), linetype="longdash", color = "black", size=0.5) +
  geom_vline(xintercept=c(as.numeric(RSS_cut$date[2268])), linetype="longdash", color = "black", size=0.5)

######################################
#### CBASS plotting in the same style
######################################

#merge and clean temp files
boing1<-read.csv("CBASS1.txt")
boing1$time<-strptime(paste(boing1$Date,paste(boing1$Th,boing1$Tm,boing1$Ts,sep=":"),sep=" "),format="%Y_%B_%d %H:%M:%S")
boing1<-boing1[order(boing1$time),]
boing1<-subset(boing1, select=c(24,8,12,16,20))

boing2<-read.csv("CBASS2.txt")
boing2$time<-strptime(paste(boing2$Date,paste(boing2$Th,boing2$Tm,boing2$Ts,sep=":"),sep=" "),format="%Y_%B_%d %H:%M:%S")
boing2<-boing2[order(boing2$time),]
boing2<-subset(boing2, select=c(24,8,12,16,20))

library(plyr)
boing1<-rename(boing1, c("T1inT"="27A", "T2inT"="29.5A", "T3inT"="32A", "T4inT"="34.5A"))
boing2<-rename(boing2, c("T1inT"="27B", "T2inT"="29.5B", "T3inT"="32B", "T4inT"="34.5B"))

library(reshape2)
boing1$time<-as.POSIXct(boing1$time)
boing1.1<-melt(boing1, id.vars=c("time"))

boing2$time<-as.POSIXct(boing2$time)
boing2.1<-melt(boing2, id.vars=c("time"))

spoing<-merge(boing1, boing2, all=T)
CBASS_full <- melt(spoing, id="time")

names(CBASS_full)[names(CBASS_full) == "time"] <- "date"
names(CBASS_full)[names(CBASS_full) == "variable"] <- "Tank"
names(CBASS_full)[names(CBASS_full) == "value"] <- "Temp"

library(openair)
CBASS_full$date<-as.POSIXct(CBASS_full$date)

#Select for experimental period within larger file, average data to 15-min intervals for plotting
CBASS_full <- na.omit(CBASS_full)
CBASS_av<-timeAverage(CBASS_full, avg.time = "15 min", type = "Tank")

#calculate mean temp during the hold
CBASS_hold<-selectByDate(CBASS_full, start = "2019-01-27", end = "2019-01-27",hour = 16:18)
CBASS_summary_data <- ddply(CBASS_hold, c("Tank"), summarise, N= length(Temp), mean = mean(Temp),
                          sd   = sd(Temp),
                          se   = sd / sqrt(N))
#plotting Fig. 1c
library(ggplot2)
cols <- c("27A" = "navyblue", "27B" = "navyblue", "29.5A" = "yellow2", "29.5B" = "yellow2", 
          "32A" = "darkorange2", "32B" = "darkorange2", "34.5A" = "red3", "34.5B" = "red3")
ggplot(CBASS_av, aes(x=date, y=Temp, col=Tank)) + scale_colour_manual(values=cols) +
  geom_line(aes(linetype=Tank)) +
  scale_linetype_manual(values=c("solid","solid","solid","solid","dashed","dashed","dashed","dashed")) +
  scale_y_continuous(limits=c(25.5,36), breaks = c(27,29,31,33,35)) + 
  scale_x_datetime(name = "Date", date_labels = "%b %d %H-h", date_breaks = "1 hour") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.3,"cm")) +
  geom_hline(yintercept=c(27,29.5,32,34.5), linetype="dashed", color = c("navyblue","yellow2","darkorange2","red3"), size=0.5) +
  geom_vline(xintercept=c(as.numeric(CBASS_av$date[5])), linetype="longdash", color = "black", size=0.5) +
  geom_vline(xintercept=c(as.numeric(CBASS_av$date[29])), linetype="longdash", color = "black", size=0.5) +
  geom_vline(xintercept=c(as.numeric(CBASS_av$date[81])), linetype="longdash", color = "black", size=0.5)

### Panels c and d were then cleaned up and combined with panels a and b using Illustrator.

