#Load packages
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

library(nlme)
library(dplyr)
library(vioplot)
library(tidyr)
library(ggplot2)
library(RInSp)
library(GLDEX)
library(truncnorm)
library(gtools)
library(RColorBrewer)

#Clear envt and import data
rm(list=ls())
dat<-read.csv(file="Herrmann_et_al_MAIN.csv",header=TRUE)

#Clean and prep data
colnames(dat)
str(dat)
colnames(dat)[colnames(dat)=="X"] <- "Observation"
#Add column that combines site and stage, generating 6 unique site/stage combos
dat$SiteSession <- paste(dat$Site,dat$Session)
#add column to hold all perch diameters with ground observations included
dat$DiamTEST=dat$Diam 
dat$DiamTEST=replace(dat$DiamTEST,is.na(dat$Diam),max(dat$Diam, na.rm=T)) #replace all NAs
#(which were ground observations, with the maximum diameter recorded - in this case 24 cm)


######################################################
#################  POPULATION-LEVEL  #################
######################################################

###PERCH HEIGHT

#Visualize

#Make a violin plot of perch heights
x1=dat$Height[dat$Site=="SL13C" & dat$Session=="first"]
x2=dat$Height[dat$Site=="SL1B" & dat$Session=="first"]
x3=dat$Height[dat$Site=="SL13C" & dat$Session=="second"]
x4=dat$Height[dat$Site=="SL1B" & dat$Session=="second"]
x5=dat$Height[dat$Site=="SL13C" & dat$Session=="third"]
x6=dat$Height[dat$Site=="SL1B" & dat$Session=="third"]

catego1=c("Pre-removal","Post-removal", "Delayed post-removal")
catego2=c("n=120","122","166","162","187","72")

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,570),
     axes=FALSE,ann=FALSE)
vioplot(x1, x2, x3, x4, x5, x6,
        col=rep(c("forestgreen","saddlebrown"),3),
        rectCol = rep(c("black","forestgreen"),3),
        add=TRUE
)
axis(1, at=c(1.5,3.5,5.5),labels=catego1,cex.axis=1.4, tick=FALSE)
axis(2, cex.axis=1)
title(ylab="Perch height (cm)", xlab="",cex.lab=1.4)
text(x=1:6,y=565, labels=catego2, cex=1)


#Population-level model for perch height mean

modelPHmean <- lme(Height~ SiteSession, data = dat, random = ~1|ID, na.action = "na.omit")
summary(modelPHmean)


#Population-level model for perch height variance

datVAR<-data.frame(cbind(dat$SiteSession,dat$Height))
colnames(datVAR)<-c("SiteSession","Height")
datVAR$Height<-as.numeric(datVAR$Height)
str(datVAR)

#calculate perch height residuals manually
datVAR <- datVAR %>% group_by(SiteSession) %>% 
  mutate(Height.avg = mean(Height))
datVAR$Height.avg.res<-abs(datVAR$Height-datVAR$Height.avg)

# Then we run an ANOVA, and post-hoc:
levene.aov.PH<-aov(Height.avg.res~SiteSession,data=datVAR)
summary(levene.aov.PH)
TukeyHSD(levene.aov.PH)


###PERCH DIAMETER

#Visualize

#Make a violin plot of perch diameters
x1=dat$Diam[dat$Site=="SL13C" & dat$Session=="first"]
x2=dat$Diam[dat$Site=="SL1B" & dat$Session=="first"]
x3=dat$Diam[dat$Site=="SL13C" & dat$Session=="second"]
x4=dat$Diam[dat$Site=="SL1B" & dat$Session=="second"]
x5=dat$Diam[dat$Site=="SL13C" & dat$Session=="third"]
x6=dat$Diam[dat$Site=="SL1B" & dat$Session=="third"]

catego1=c("Pre-removal","Post-removal", "Delayed post-removal")
catego2=c("on ground = 4","5","3","7","6","0")

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(0,27),
     axes=FALSE,ann=FALSE)
vioplot(x1, x2, x3, x4, x5, x6,
        col=rep(c("forestgreen","saddlebrown"),3),
        rectCol = rep(c("black","forestgreen"),3),
        add=TRUE
)
axis(1, at=c(1.5,3.5,5.5),labels=catego1,cex.axis=1.4, tick=FALSE)
axis(2, cex.axis=1)
title(ylab="Perch diameter (cm)", xlab="",cex.lab=1.4)
text(x=1:6,y=26, labels=catego2, cex=1)


#Population-level model for perch diameter mean

modePDmean1 <- lme(DiamTEST~ SiteSession, data = dat, random = ~1|ID, na.action = "na.omit")
summary(modePDmean1)


#Population-level model for perch diameter variance

datVAR<-data.frame(cbind(dat$SiteSession,dat$DiamTEST))
colnames(datVAR)<-c("SiteSession","Diam")
datVAR$Diam<-as.numeric(datVAR$Diam)
str(datVAR)

#calculate perch height residuals manually
datVAR <- datVAR %>% group_by(SiteSession) %>% 
  mutate(Diam.avg = mean(Diam))
datVAR$Diam.avg.res<-abs(datVAR$Diam-datVAR$Diam.avg)

# Then we run an ANOVA, and post-hoc:
levene.aov.PD<-aov(Diam.avg.res~SiteSession,data=datVAR)
summary(levene.aov.PD)
TukeyHSD(levene.aov.PD)



######################################################
#################  INDIVIDUAL-LEVEL  #################
######################################################

#####################################################################
#STEP ONE IS TO GENERATE A NEW DATA FRAME WITH ONE ROW PER INDIVIDUAL THAT DOCUMENTS
#OBSERVATION COUNT, ph MEAN, pd MEAN AND LATERAL MOVEMENT FOR ALL SAMPLING SESSIONS


#Generate data frame indivs1 which summarizes individual level habitat use across sessions
IDs=as.character(unique(dat$ID))
Site=rep(NA,length(IDs))
Sex=rep(NA,length(IDs))

PHmeansPRE=rep(NA,length(IDs))
PDmeansPRE=rep(NA,length(IDs))
sizePRE=rep(NA,length(IDs))

PHmeansPOST=rep(NA,length(IDs))
PDmeansPOST=rep(NA,length(IDs))
sizePOST=rep(NA,length(IDs))

PHmeansDPOST=rep(NA,length(IDs))
PDmeansDPOST=rep(NA,length(IDs))
sizeDPOST=rep(NA,length(IDs))


for (i in 1:length(IDs)){
  heightsPRE=dat$Height[dat$ID==IDs[i] & dat$Session=="first"]
  heightsPOST=dat$Height[dat$ID==IDs[i] & dat$Session=="second"]
  heightsDPOST=dat$Height[dat$ID==IDs[i] & dat$Session=="third"]
  PHmeansPRE[i]=mean(heightsPRE)
  PHmeansPOST[i]=mean(heightsPOST)
  PHmeansDPOST[i]=mean(heightsDPOST)
  
  diamsPRE=dat$DiamTEST[dat$ID==IDs[i] & dat$Session=="first"]
  diamsPOST=dat$DiamTEST[dat$ID==IDs[i] & dat$Session=="second"]
  diamsDPOST=dat$DiamTEST[dat$ID==IDs[i] & dat$Session=="third"]
  PDmeansPRE[i]=mean(diamsPRE)
  PDmeansPOST[i]=mean(diamsPOST)
  PDmeansDPOST[i]=mean(diamsDPOST)
  
  
  sizePRE[i]=length(heightsPRE)
  sizePOST[i]=length(heightsPOST)
  sizeDPOST[i]=length(heightsDPOST)
  Site[i]=names(which.max(table(dat$Site[dat$ID==IDs[i]])))
  Sex[i]=names(which.max(table(dat$Sex[dat$ID==IDs[i]])))
  
}


indivs1=as.data.frame(cbind(IDs,Site,Sex, PHmeansPRE, PHmeansPOST, PHmeansDPOST,
                            PDmeansPRE, PDmeansPOST, PDmeansDPOST,
                            sizePRE,sizePOST,sizeDPOST))

indivs1$PHmeansPRE=as.numeric(PHmeansPRE)
indivs1$PHmeansPOST=as.numeric(PHmeansPOST)
indivs1$PHmeansDPOST=as.numeric(PHmeansDPOST)
indivs1$PDmeansPRE=as.numeric(PDmeansPRE)
indivs1$PDmeansPOST=as.numeric(PDmeansPOST)
indivs1$PDmeansDPOST=as.numeric(PDmeansDPOST)
indivs1$sizePRE=as.numeric(sizePRE)
indivs1$sizePOST=as.numeric(sizePOST)
indivs1$sizeDPOST=as.numeric(sizeDPOST)

#so far so good. now to add a measure of lateral movement -
#minimum number of grid squares traversed (calculated manually)

#add lateral movement metrics to indivs1
lateral<-read.csv(file="Herrmann_et_al_LateralMovement.csv",header=TRUE)
indivs1<-left_join(indivs1, lateral,by="IDs")

#these extra rows are needed to look at differences across pre and post sessions
indivs1$Height_Mean_Shift=indivs1$PHmeansPOST-indivs1$PHmeansPRE
indivs1$Diameter_Mean_Shift=indivs1$PDmeansPOST-indivs1$PDmeansPRE
indivs1$LateralMovement_Shift=indivs1$latPOST-indivs1$latPRE
indivs1$SummerObservations=indivs1$sizePRE+indivs1$sizePOST


str(indivs1)

############################
#####  VISUALIZE!!!!  !#####
############################

###
#ALL INDIVIDUALS ACROSS ALL SESSIONS
###

#AvgHeight

indivsHEIGHT=indivs1
indivsHEIGHT=gather(indivsHEIGHT,'PHmeansPRE','PHmeansPOST','PHmeansDPOST',key = "Session1",value="MeanHeight")
indivsHEIGHT=gather(indivsHEIGHT,'sizePRE','sizePOST','sizeDPOST',key = "Session2",value="Observations")


A=indivsHEIGHT$Session1=="PHmeansPRE" & indivsHEIGHT$Session2=="sizePRE"
B=indivsHEIGHT$Session1=="PHmeansPOST" & indivsHEIGHT$Session2=="sizePOST"
C=indivsHEIGHT$Session1=="PHmeansDPOST" & indivsHEIGHT$Session2=="sizeDPOST"
indivsHEIGHT=subset(indivsHEIGHT, A|B|C)
indivsHEIGHT$Session1<-factor(indivsHEIGHT$Session1,ordered=TRUE,
                         levels = as.factor(c("PHmeansPRE","PHmeansPOST","PHmeansDPOST")))
levels(indivsHEIGHT$Session1) <- c('Pre-removal', 'Post-removal', ' Delayed post-removal')


dodge <- position_dodge(width = .2)
indivPLOTheight <- ggplot(indivsHEIGHT,aes(colour=Site,group=IDs,alpha=Site,shape=Sex))+
  geom_point(aes(x=Session1,y=MeanHeight,size=Observations),position = dodge)+
  geom_line(aes(x=Session1,y=MeanHeight),position = dodge)

indivPLOTheight + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                values=c("forestgreen", "saddlebrown")) +
  scale_alpha_manual(breaks = c("SL13C", "SL1B"),
                     values=c(0.3,0.3),guide=F)+
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=14),
        legend.position = "none")+
  theme(axis.title.y=element_text(margin=margin(r=10)))+
  coord_cartesian(ylim = c(0, 400))


#AvgDiam

indivsDIAM=indivs1
indivsDIAM=gather(indivsDIAM,'PDmeansPRE','PDmeansPOST','PDmeansDPOST',key = "Session1",value="MeanDiam")
indivsDIAM=gather(indivsDIAM,'sizePRE','sizePOST','sizeDPOST',key = "Session2",value="Observations")


A=indivsDIAM$Session1=="PDmeansPRE" & indivsDIAM$Session2=="sizePRE"
B=indivsDIAM$Session1=="PDmeansPOST" & indivsDIAM$Session2=="sizePOST"
C=indivsDIAM$Session1=="PDmeansDPOST" & indivsDIAM$Session2=="sizeDPOST"
indivsDIAM=subset(indivsDIAM, A|B|C)
indivsDIAM$Session1<-factor(indivsDIAM$Session1,ordered=TRUE,
                         levels = as.factor(c("PDmeansPRE","PDmeansPOST","PDmeansDPOST")))
levels(indivsDIAM$Session1) <- c('Pre-removal', 'Post-removal', ' Delayed post-removal')


dodge <- position_dodge(width = .25)
indivPLOTdiam <- ggplot(indivsDIAM,aes(colour=Site,group=IDs,alpha=Site,shape=Sex))+
  geom_point(aes(x=Session1,y=MeanDiam,size=Observations),position = dodge)+
  geom_line(aes(x=Session1,y=MeanDiam),position = dodge)


indivPLOTdiam + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                values=c("forestgreen", "saddlebrown")) +
  scale_alpha_manual(breaks = c("SL13C", "SL1B"),
                     values=c(0.3,0.3),guide=F)+
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=14),
        legend.position = "none")+
  theme(axis.title.y=element_text(margin=margin(r=10)))+
  coord_cartesian(ylim = c(0, 12))


#Lateral movement

indivsLAT=indivs1
indivsLAT=gather(indivsLAT,'latPRE','latPOST','latDPOST',key = "Session1",value="LateralMovement")
indivsLAT=gather(indivsLAT,'sizePRE','sizePOST','sizeDPOST',key = "Session2",value="Observations")


A=indivsLAT$Session1=="latPRE" & indivsLAT$Session2=="sizePRE"
B=indivsLAT$Session1=="latPOST" & indivsLAT$Session2=="sizePOST"
C=indivsLAT$Session1=="latDPOST" & indivsLAT$Session2=="sizeDPOST"
indivsLAT=subset(indivsLAT, A|B|C)
indivsLAT$Session1<-factor(indivsLAT$Session1,ordered=TRUE,
                         levels = as.factor(c("latPRE","latPOST","latDPOST")))
levels(indivsLAT$Session1) <- c('Pre-removal', 'Post-removal', ' Delayed post-removal')


dodge <- position_dodge(width = .2)
indivPLOTlat <- ggplot(indivsLAT,aes(colour=Site,group=IDs,alpha=Site,shape=Sex))+
  geom_point(aes(x=Session1,y=LateralMovement,size=Observations),position = dodge)+
  geom_line(aes(x=Session1,y=LateralMovement),position = dodge)


indivPLOTlat + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                values=c("forestgreen", "saddlebrown")) +
  scale_alpha_manual(breaks = c("SL13C", "SL1B"),
                     values=c(0.3,0.3),guide=F)+
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=14),
        legend.position = "none")+
  theme(axis.title.y=element_text(margin=margin(r=10)))+
  coord_cartesian(ylim = c(0, 16))


###
###INDIVIDUAL SHIFTS FROM JULY-AUGUST
###

#HEIGHT

indivs2=subset(indivs1,indivs1$sizePRE>=1 & indivs1$sizePOST>=1) #min of 8 observations total
indivs2$IDs<-factor(indivs2$IDs,ordered=TRUE, levels = indivs2$IDs[order(indivs2$Height_Mean_Shift)])

indivPLOT_heightshifts <- ggplot(indivs2,aes(shape=Sex,colour=Site))+
  geom_point(aes(x=IDs,y=Height_Mean_Shift,size=SummerObservations))+
  geom_hline(aes(yintercept=0),linetype = 2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
indivPLOT_heightshifts + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                values=c("forestgreen", "saddlebrown")) +
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(axis.title.y=element_text(margin=margin(r=10)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#PD

indivPLOT_diamshifts <- ggplot(indivs2,aes(shape=Sex,colour=Site))+
  geom_point(aes(x=IDs,y=Diameter_Mean_Shift,size=SummerObservations))+
  geom_hline(aes(yintercept=0),linetype = 2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
indivPLOT_diamshifts + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                 values=c("forestgreen", "saddlebrown")) +
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(axis.title.y=element_text(margin=margin(r=10)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Lateral movement

indivPLOT_latshifts <- ggplot(indivs2,aes(shape=Sex,colour=Site))+
  geom_point(aes(x=IDs,y=LateralMovement_Shift,size=SummerObservations))+
  geom_hline(aes(yintercept=0),linetype = 2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
indivPLOT_latshifts + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                 values=c("forestgreen", "saddlebrown")) +
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))+
  theme(axis.title.y=element_text(margin=margin(r=10)))



###
###OLD VS NEW INDIVIDUALS DURING DELAYED POST-REMOVAL
###

indivs3=subset(indivs1,!is.na(indivs1$PHmeansDPOST) & indivs1$sizeDPOST>=1)

for(i in 1:length(indivs3$IDs)){
  if(is.na(indivs3$PHmeansPRE[i])){indivs3$Age[i]="New Lizards"}
  else{indivs3$Age[i]="Old Lizards"}
}
indivs3$Age<-factor(indivs3$Age,ordered=TRUE,
                    levels = as.factor(c("Old Lizards","New Lizards")))

colnames(indivs3)[which(colnames(indivs3)=="PHmeansDPOST")]="MeanHeight"
colnames(indivs3)[which(colnames(indivs3)=="PDmeansDPOST")]="MeanDiam"
colnames(indivs3)[which(colnames(indivs3)=="PHvarsDPOST")]="HeightSD"
colnames(indivs3)[which(colnames(indivs3)=="PDvarsDPOST")]="DiamSD"
colnames(indivs3)[which(colnames(indivs3)=="latDPOST")]="LateralMovement"
colnames(indivs3)[which(colnames(indivs3)=="sizeDPOST")]="Observations"

#Avg PH

dodge <- position_dodge(width = .3)

indivPLOT_OldVsNewHeight <- ggplot(indivs3,aes(x=Age,y=MeanHeight,
                                 colour=Site, shape=Sex, size=Observations))+
  geom_point(alpha=0.3,position = dodge)


indivPLOT_OldVsNewHeight + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                values=c("forestgreen", "saddlebrown")) +
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=14),
        legend.position = "none")+
  theme(axis.title.y=element_text(margin=margin(r=10)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  coord_cartesian(ylim = c(0, 360)) 



#Avg PD

indivPLOT_OldVsNewDiam <- ggplot(indivs3,aes(x=Age,y=MeanDiam, shape=Sex,
                                  colour=Site, size=Observations))+
  geom_point(alpha=0.3,position = dodge)


indivPLOT_OldVsNewDiam + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                 values=c("forestgreen", "saddlebrown")) +
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=14),
        legend.position = "none")+
  theme(axis.title.y=element_text(margin=margin(r=10)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  coord_cartesian(ylim = c(0, 8)) 


#Lateral movement

indivPLOT_OldVsNewLat <- ggplot(indivs3,aes(x=Age,y=LateralMovement,shape=Sex,
                                  colour=Site, size=Observations))+
  geom_point(alpha=0.3,position = dodge)


indivPLOT_OldVsNewLat + scale_color_manual(breaks = c("SL13C", "SL1B"),
                                 values=c("forestgreen", "saddlebrown")) +
  scale_shape_manual(breaks = c("M", "F"),
                     values=c(16,17)) +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=14),
        legend.position = "none")+
  theme(axis.title.y=element_text(margin=margin(r=10)))+
  coord_cartesian(ylim = c(0, 15)) 


###
#RANK ORDER TESTS FOR Avg PH, Avg, PD, and Lat
###

indivs13C=subset(indivs1,indivs1$Site=="SL13C")
indivs1B=subset(indivs1,indivs1$Site=="SL1B")

a=indivs13C$PHmeansPRE[!is.na(indivs13C$PHmeansPRE)]
b=indivs1B$PHmeansPRE[!is.na(indivs1B$PHmeansPRE)]
c=indivs13C$PHmeansPOST[!is.na(indivs13C$PHmeansPOST)]
d=indivs1B$PHmeansPOST[!is.na(indivs1B$PHmeansPOST)]
e=indivs13C$PHmeansDPOST[!is.na(indivs13C$PHmeansDPOST)]
f=indivs1B$PHmeansDPOST[!is.na(indivs1B$PHmeansDPOST)]

wilcox.test(a,b)
wilcox.test(c,d)
wilcox.test(e,f)


a=indivs13C$PDmeansPRE[!is.na(indivs13C$PDmeansPRE)]
b=indivs1B$PDmeansPRE[!is.na(indivs1B$PDmeansPRE)]
c=indivs13C$PDmeansPOST[!is.na(indivs13C$PDmeansPOST)]
d=indivs1B$PDmeansPOST[!is.na(indivs1B$PDmeansPOST)]
e=indivs13C$PDmeansDPOST[!is.na(indivs13C$PDmeansDPOST)]
f=indivs1B$PDmeansDPOST[!is.na(indivs1B$PDmeansDPOST)]

wilcox.test(a,b)
wilcox.test(c,d)
wilcox.test(e,f)


a=indivs13C$latPRE[!is.na(indivs13C$latPRE)]
b=indivs1B$latPRE[!is.na(indivs1B$latPRE)]
c=indivs13C$latPOST[!is.na(indivs13C$latPOST)]
d=indivs1B$latPOST[!is.na(indivs1B$latPOST)]
e=indivs13C$latDPOST[!is.na(indivs13C$latDPOST)]
f=indivs1B$latDPOST[!is.na(indivs1B$latDPOST)]

wilcox.test(a,b)
wilcox.test(c,d)
wilcox.test(e,f)



###
#Individual Perch Height Specialization
###


#First reformat data - need to adjust because zeroes don't represent missing data
#subset data into sites and sessions
july13C=subset(dat,dat$SiteSession=="SL13C first")
aug13C=subset(dat,dat$SiteSession=="SL13C second")
mar13C=subset(dat,dat$SiteSession=="SL13C third")
july1B=subset(dat,dat$SiteSession=="SL1B first")
aug1B=subset(dat,dat$SiteSession=="SL1B second")
mar1B=subset(dat,dat$SiteSession=="SL1B third")


#here's a function that will cbind data frames of uneven length...
#and automatically pad the shorter ones with NAs

cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix(NA, mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}

#and here's a loop to reformat the data so it will play nicely with RInSp

###JULY13C (Pre-removal on the one-species island) ###

datPREP=july13C
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE HEIGHTS
  storeheights=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storeheights=c(storeheights, datPREP$Height[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storeheights=data.frame(storeheights)
  datNEW=cbindPad(datNEW,storeheights)
}
datNEW=data.frame(t(datNEW))
datNEW=datNEW+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
#...from missing values
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once


RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
july13C_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")


###AUGUST13C (Post-removal on the one-species island)###

datPREP=aug13C
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE HEIGHTS
  storeheights=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storeheights=c(storeheights, datPREP$Height[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storeheights=data.frame(storeheights)
  datNEW=cbindPad(datNEW,storeheights)
}
datNEW=data.frame(t(datNEW))
datNEW=datNEW+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
#...from missing values
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once


RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
aug13C_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")


###MARCH13C (Delayed Post-removal on the one-species island)###

datPREP=mar13C
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE HEIGHTS
  storeheights=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storeheights=c(storeheights, datPREP$Height[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storeheights=data.frame(storeheights)
  datNEW=cbindPad(datNEW,storeheights)
}
datNEW=data.frame(t(datNEW))
datNEW=datNEW+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
#...from missing values
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once

RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
mar13C_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")


###JULY1B (Pre-removal on the two-species island)###

datPREP=july1B
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE HEIGHTS
  storeheights=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storeheights=c(storeheights, datPREP$Height[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storeheights=data.frame(storeheights)
  datNEW=cbindPad(datNEW,storeheights)
}
datNEW=data.frame(t(datNEW))
datNEW=datNEW+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
#...from missing values
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once

RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
july1B_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")


###AUGUST1B (Post-removal on the two-species island)###

datPREP=aug1B
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE HEIGHTS
  storeheights=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storeheights=c(storeheights, datPREP$Height[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storeheights=data.frame(storeheights)
  datNEW=cbindPad(datNEW,storeheights)
}
datNEW=data.frame(t(datNEW))
datNEW=datNEW+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
#...from missing values
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once

RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
aug1B_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")


###MARCH1B (Delayed post-removal on the two-species island)###

datPREP=mar1B
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE HEIGHTS
  storeheights=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storeheights=c(storeheights, datPREP$Height[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
  storeheights=data.frame(storeheights)
  datNEW=cbindPad(datNEW,storeheights)
}
datNEW=data.frame(t(datNEW))
datNEW=datNEW+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
#...from missing values
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once

RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
mar1B_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")


#SO ARE THESE DIFFERENCES IN INDIVIDUAL SPECIALIZATION MEANINGFUL???...
#..or can they be explained my sampling error?

#Generate a null distribution of observed WIC/TNW differences across sampling sessions,
#based on our own sample sizes and with different known true WIC/TNW differences,
#to determine how reliable our empirically observed metrics are.

str(indivs1) #this dataframe was generated earlier

#Choose populations/sessions to compare using empirical sample sizes
popn1<-subset(indivs1,indivs1$Site=="SL13C")$sizePRE
popn2<-subset(indivs1,indivs1$Site=="SL13C")$sizePOST
popn3<-subset(indivs1,indivs1$Site=="SL13C")$sizeDPOST
popn4<-subset(indivs1,indivs1$Site=="SL1B")$sizePRE
popn5<-subset(indivs1,indivs1$Site=="SL1B")$sizePOST
popn6<-subset(indivs1,indivs1$Site=="SL1B")$sizeDPOST


#drop zeroes
popn1<-as.vector(fun.zero.omit(popn1))
popn2<-as.vector(fun.zero.omit(popn2))
popn3<-as.vector(fun.zero.omit(popn3))
popn4<-as.vector(fun.zero.omit(popn4))
popn5<-as.vector(fun.zero.omit(popn5))
popn6<-as.vector(fun.zero.omit(popn6))


#mu_s is the vector of empirical perch heights across sampling sessions
mu_s <- c(180,170,149,169,151,279)
#sig_BIC is the vector of SDs determining the perch height divergence BETWEEN individuals...
#within a sampling session
sig_BIC <- c(70,70,70,70,70,70)
#and assign a within-individual SD to each sampling session
sig_WIC <- c(100,100,100,100,100,100)


#START HERE WITH LOOPING EVERYTHING TOGETHER TO RUN FULL SET OF SIMULATIONS

N=1000 #set number of simulations
storeSIMS<-matrix(nrow = N,ncol=6)

for (z in 1:N){
  
  #for each individual, draw it's "baseline" perch height appropriate sampling session mean and sd
  popn1baseline <- rep(0,length(popn1))
  popn2baseline <- rep(0,length(popn2))
  popn3baseline <- rep(0,length(popn3))
  popn4baseline <- rep(0,length(popn4))
  popn5baseline <- rep(0,length(popn5))
  popn6baseline <- rep(0,length(popn6))
  
  L <-0 # lower bound on perch height response variable
  for (j in 1:length(popn1)){
    popn1baseline[j] <- rtruncnorm(1,a=L,b=Inf,mu_s[1], sig_BIC[1])
  }
  
  for (k in 1:length(popn2)){
    popn2baseline[k] <- rtruncnorm(1,a=L,b=Inf,mu_s[2], sig_BIC[2])
  }
  
  for (j in 1:length(popn3)){
    popn3baseline[j] <- rtruncnorm(1,a=L,b=Inf,mu_s[3], sig_BIC[3])
  }
  
  for (k in 1:length(popn4)){
    popn4baseline[k] <- rtruncnorm(1,a=L,b=Inf,mu_s[4], sig_BIC[4])
  }
  
  for (j in 1:length(popn5)){
    popn5baseline[j] <- rtruncnorm(1,a=L,b=Inf,mu_s[5], sig_BIC[5])
  }
  
  for (k in 1:length(popn6)){
    popn6baseline[k] <- rtruncnorm(1,a=L,b=Inf,mu_s[6], sig_BIC[6])
  }
  
  #Now simulate the fake data for each population and store
  samplings1<-data.frame()
  for(i in 1:length(popn1)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn1[i],a=L,b=Inf,popn1baseline[i], sig_WIC[1])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings1=cbindPad(samplings1,storeheights)
  }
  
  samplings1=data.frame(t(samplings1))
  samplings1=samplings1+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
  #...from missing values
  samplings1[is.na(samplings1)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings1 <- samplings1[samplings1$X2!=0,] #remove individuals who were observed only once
  
  
  samplings2<-data.frame()
  for(i in 1:length(popn2)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn2[i],a=L,b=Inf,popn2baseline[i], sig_WIC[2])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings2=cbindPad(samplings2,storeheights)
  }
  
  samplings2=data.frame(t(samplings2))
  samplings2=samplings2+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
  #...from missing values
  samplings2[is.na(samplings2)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings2 <- samplings2[samplings2$X2!=0,] #remove individuals who were observed only once
  
  
  samplings3<-data.frame()
  for(i in 1:length(popn3)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn3[i],a=L,b=Inf,popn3baseline[i], sig_WIC[3])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings3=cbindPad(samplings3,storeheights)
  }
  
  samplings3=data.frame(t(samplings3))
  samplings3=samplings3+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
  #...from missing values
  samplings3[is.na(samplings3)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings3 <- samplings3[samplings3$X2!=0,] #remove individuals who were observed only once
  
  
  samplings4<-data.frame()
  for(i in 1:length(popn4)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn4[i],a=L,b=Inf,popn4baseline[i], sig_WIC[4])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings4=cbindPad(samplings4,storeheights)
  }
  
  samplings4=data.frame(t(samplings4))
  samplings4=samplings4+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
  #...from missing values
  samplings4[is.na(samplings4)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings4 <- samplings4[samplings4$X2!=0,] #remove individuals who were observed only once
  
  
  samplings5<-data.frame()
  for(i in 1:length(popn5)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn5[i],a=L,b=Inf,popn5baseline[i], sig_WIC[5])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings5=cbindPad(samplings5,storeheights)
  }
  
  samplings5=data.frame(t(samplings5))
  samplings5=samplings5+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
  #...from missing values
  samplings5[is.na(samplings5)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings5 <- samplings5[samplings5$X2!=0,] #remove individuals who were observed only once
  
  
  samplings6<-data.frame()
  for(i in 1:length(popn6)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn6[i],a=L,b=Inf,popn6baseline[i], sig_WIC[6])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings6=cbindPad(samplings6,storeheights)
  }
  
  samplings6=data.frame(t(samplings6))
  samplings6=samplings6+1 #as previously stated, this is to distinguish PH=0 (which are converted to 1s)...
  #...from missing values
  samplings6[is.na(samplings6)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings6 <- samplings6[samplings6$X2!=0,] #remove individuals who were observed only once
  
  
  
  #Calculate and store WIC/TNW for each population
  RIS1= import.RInSp(samplings1, col.header=TRUE, data.type="double")
  RIS2= import.RInSp(samplings2, col.header=TRUE, data.type="double")
  RIS3= import.RInSp(samplings3, col.header=TRUE, data.type="double")
  RIS4= import.RInSp(samplings4, col.header=TRUE, data.type="double")
  RIS5= import.RInSp(samplings5, col.header=TRUE, data.type="double")
  RIS6= import.RInSp(samplings6, col.header=TRUE, data.type="double")
  
  RIS1calc=WTcMC(RIS1, replicates = 999, weight = "N_items")
  RIS2calc=WTcMC(RIS2, replicates = 999, weight = "N_items")
  RIS3calc=WTcMC(RIS3, replicates = 999, weight = "N_items")
  RIS4calc=WTcMC(RIS4, replicates = 999, weight = "N_items")
  RIS5calc=WTcMC(RIS5, replicates = 999, weight = "N_items")
  RIS6calc=WTcMC(RIS6, replicates = 999, weight = "N_items")
  
  storeSIMS[z,1]=as.numeric(RIS1calc$WonT)
  storeSIMS[z,2]=as.numeric(RIS2calc$WonT)
  storeSIMS[z,3]=as.numeric(RIS3calc$WonT)
  storeSIMS[z,4]=as.numeric(RIS4calc$WonT)
  storeSIMS[z,5]=as.numeric(RIS5calc$WonT)
  storeSIMS[z,6]=as.numeric(RIS6calc$WonT)
  
}

storeSIMS=as.data.frame(storeSIMS)
colnames(storeSIMS)=c("July13C","Aug13C","Mar13C","July1B","Aug1B","Mar1B")
#write.csv(storeSIMS,"PH_WIC_sims.csv")

#So now the data frame storeSIMS contains WIC/TNW for N simulated samples from EACH sampling period...
#...and we can compare differences to our empirical WIC/TNW values to evaluate whether our empirical...
#...values are more different than we would expect by chance (if the underlying niche parameters were...
#...in reality identical, as they were here in these simulations [i.e. all sampling sessions had...
#..the same true BIC and WIC])

hist(storeSIMS$July13C)
hist(storeSIMS$Aug13C)
hist(storeSIMS$Mar13C)
hist(storeSIMS$July1B)
hist(storeSIMS$Aug1B)
hist(storeSIMS$Mar1B)

#Here I'm checking the largest empirical differences in WIC/TNW between island-stage combos...
#First with a histrogram and then, if it looks close, actually counting to see if the empirical...
#...difference does indeed fall in the outer 5% of simulated differences

hist(storeSIMS$Mar13C-storeSIMS$July13C, breaks=50) #nope
hist(storeSIMS$Aug1B-storeSIMS$July1B, breaks=50) #oooo, maybe...
sum((storeSIMS$Aug1B-storeSIMS$July1B)>=.3) #yup!!
hist(storeSIMS$Aug1B-storeSIMS$July13C, breaks=50) #nope
hist(storeSIMS$Mar13C-storeSIMS$July1B, breaks=50) #maayyyyybe....
sum((storeSIMS$Mar13C-storeSIMS$July1B)>=.18) #CLOSE, but not quite


###
#Individual Perch Diameter Specialization
###

#here's a loop to reformat the data so it will play nicely with RInSp

###JULY13C (Pre-removal on the one-species island)###

datPREP=july13C
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE DIAMS
  storediams=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storediams=c(storediams, datPREP$DiamTEST[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storediams=data.frame(storediams)
  datNEW=cbindPad(datNEW,storediams)
}
datNEW=data.frame(t(datNEW))
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once


RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
july13C_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")

###AUGUST13C (Post-removal on the one-species island) ###

datPREP=aug13C
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE DIAMS
  storediams=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storediams=c(storediams, datPREP$DiamTEST[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storediams=data.frame(storediams)
  datNEW=cbindPad(datNEW,storediams)
}
datNEW=data.frame(t(datNEW))
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once


RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
aug13C_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")

###MARCH13C (Delayed post-removal on the one-species island)###

datPREP=mar13C
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE DIAMS
  storediams=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storediams=c(storediams, datPREP$DiamTEST[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storediams=data.frame(storediams)
  datNEW=cbindPad(datNEW,storediams)
}
datNEW=data.frame(t(datNEW))
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once

RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
mar13C_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")


###JULY1B (Pre-removal on the two-species island)###

datPREP=july1B
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE DIAMS
  storediams=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storediams=c(storediams, datPREP$DiamTEST[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storediams=data.frame(storediams)
  datNEW=cbindPad(datNEW,storediams)
}
datNEW=data.frame(t(datNEW))
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once

RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
july1B_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")



###AUGUST1B (Post-removal on the two-species island)###

datPREP=aug1B
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE DIAMS
  storediams=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storediams=c(storediams, datPREP$DiamTEST[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CBINDPAD
  storediams=data.frame(storediams)
  datNEW=cbindPad(datNEW,storediams)
}
datNEW=data.frame(t(datNEW))
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once

RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
aug1B_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")


###MARCH1B (Delayed post-removal on the two-species island)###

datPREP=mar1B
individuals=unique(datPREP$ID)
datNEW=data.frame() #This data frame will collect all the reformatted perch heights

for(i in 1:length(individuals)){
  ##HERE SET UP A VECTOR TO STORE DIAMS
  storediams=as.numeric()
  for(j in 1:nrow(datPREP)){
    if(individuals[i]==datPREP$ID[j]){
      storediams=c(storediams, datPREP$DiamTEST[j])
    }else{}
  }
  ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
  storediams=data.frame(storediams)
  datNEW=cbindPad(datNEW,storediams)
}
datNEW=data.frame(t(datNEW))
datNEW[is.na(datNEW)] <- 0 #rewrite all NAs (missing values) as 0s
datNEW <- datNEW[datNEW$X2!=0,] #remove individuals who were observed only once

RIS= import.RInSp(datNEW, col.header=TRUE, data.type="double")
#Calculating WIC/TNW
mar1B_RIS_MC_weighted = WTcMC(RIS, replicates = 999, weight = "N_items")



#SO ARE THESE DIFFERENCES IN INDIVIDUAL SPECIALIZATION MEANINGFUL???

#Generate a null distribution of observed WIC/TNW differences across sampling sessions,
#based on our own sample sizes and with different known true WIC/TNW differences,
#to determine how reliable our empirically observed metrics are.

str(indivs1) #this dataframe ws generated earlier based on raw data

#Choose populations/sessions to compare using empirical sample sizes
popn1<-subset(indivs1,indivs1$Site=="SL13C")$sizePRE
popn2<-subset(indivs1,indivs1$Site=="SL13C")$sizePOST
popn3<-subset(indivs1,indivs1$Site=="SL13C")$sizeDPOST
popn4<-subset(indivs1,indivs1$Site=="SL1B")$sizePRE
popn5<-subset(indivs1,indivs1$Site=="SL1B")$sizePOST
popn6<-subset(indivs1,indivs1$Site=="SL1B")$sizeDPOST


#drop zeroes
popn1<-as.vector(fun.zero.omit(popn1))
popn2<-as.vector(fun.zero.omit(popn2))
popn3<-as.vector(fun.zero.omit(popn3))
popn4<-as.vector(fun.zero.omit(popn4))
popn5<-as.vector(fun.zero.omit(popn5))
popn6<-as.vector(fun.zero.omit(popn6))


#mu_s is the vector of empirical mean perch diameters across sampling sessions
mu_s <- c(4.5,3.6,5,5.7,6.1,2.6)
#sig_BIC is the vector of SDs determining the perch diam divergence BETWEEN individuals...
#within a sampling session
sig_BIC <- c(2.5,2.5,2.5,2.5,2.5,2.5)
#and assign a within-individual SD to each sampling session
sig_WIC <- c(6.7,6.7,6.7,6.7,6.7,6.7)


#START HERE WITH LOOPING EVERYTHING TOGETHER TO RUN FULL SET OF SIMULATIONS

N=1000 #set number of simulations
storeSIMS<-matrix(nrow = N,ncol=6)

for (z in 1:N){
  
  #for each individual, draw it's "baseline" perch height appropriate sampling session mean and sd
  popn1baseline <- rep(0,length(popn1))
  popn2baseline <- rep(0,length(popn2))
  popn3baseline <- rep(0,length(popn3))
  popn4baseline <- rep(0,length(popn4))
  popn5baseline <- rep(0,length(popn5))
  popn6baseline <- rep(0,length(popn6))
  
  L <-0 # lower bound on perch height response variable
  for (j in 1:length(popn1)){
    popn1baseline[j] <- rtruncnorm(1,a=L,b=Inf,mu_s[1], sig_BIC[1])
  }
  
  for (k in 1:length(popn2)){
    popn2baseline[k] <- rtruncnorm(1,a=L,b=Inf,mu_s[2], sig_BIC[2])
  }
  
  for (j in 1:length(popn3)){
    popn3baseline[j] <- rtruncnorm(1,a=L,b=Inf,mu_s[3], sig_BIC[3])
  }
  
  for (k in 1:length(popn4)){
    popn4baseline[k] <- rtruncnorm(1,a=L,b=Inf,mu_s[4], sig_BIC[4])
  }
  
  for (j in 1:length(popn5)){
    popn5baseline[j] <- rtruncnorm(1,a=L,b=Inf,mu_s[5], sig_BIC[5])
  }
  
  for (k in 1:length(popn6)){
    popn6baseline[k] <- rtruncnorm(1,a=L,b=Inf,mu_s[6], sig_BIC[6])
  }
  
  #Now simulate the fake data for each population and store
  samplings1<-data.frame()
  for(i in 1:length(popn1)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn1[i],a=L,b=Inf,popn1baseline[i], sig_WIC[1])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings1=cbindPad(samplings1,storeheights)
  }
  
  samplings1=data.frame(t(samplings1))
  samplings1[is.na(samplings1)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings1 <- samplings1[samplings1$X2!=0,] #remove individuals who were observed only once
  
  
  samplings2<-data.frame()
  for(i in 1:length(popn2)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn2[i],a=L,b=Inf,popn2baseline[i], sig_WIC[2])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings2=cbindPad(samplings2,storeheights)
  }
  
  samplings2=data.frame(t(samplings2))
  samplings2[is.na(samplings2)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings2 <- samplings2[samplings2$X2!=0,] #remove individuals who were observed only once
  
  
  samplings3<-data.frame()
  for(i in 1:length(popn3)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn3[i],a=L,b=Inf,popn3baseline[i], sig_WIC[3])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings3=cbindPad(samplings3,storeheights)
  }
  
  samplings3=data.frame(t(samplings3))
  samplings3[is.na(samplings3)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings3 <- samplings3[samplings3$X2!=0,] #remove individuals who were observed only once
  
  
  samplings4<-data.frame()
  for(i in 1:length(popn4)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn4[i],a=L,b=Inf,popn4baseline[i], sig_WIC[4])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings4=cbindPad(samplings4,storeheights)
  }
  
  samplings4=data.frame(t(samplings4))
  samplings4[is.na(samplings4)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings4 <- samplings4[samplings4$X2!=0,] #remove individuals who were observed only once
  
  
  samplings5<-data.frame()
  for(i in 1:length(popn5)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn5[i],a=L,b=Inf,popn5baseline[i], sig_WIC[5])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings5=cbindPad(samplings5,storeheights)
  }
  
  samplings5=data.frame(t(samplings5))
  samplings5[is.na(samplings5)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings5 <- samplings5[samplings5$X2!=0,] #remove individuals who were observed only once
  
  
  samplings6<-data.frame()
  for(i in 1:length(popn6)){
    ##HERE SET UP A VECTOR TO STORE HEIGHTS
    storeheights=as.numeric()
    storeheights=rtruncnorm(popn6[i],a=L,b=Inf,popn6baseline[i], sig_WIC[6])
    ###HERE ADD THE NEW ROW TO ALL OLD ROWS USIND CPINDPAD
    storeheights=data.frame(storeheights)
    samplings6=cbindPad(samplings6,storeheights)
  }
  
  samplings6=data.frame(t(samplings6))
  samplings6[is.na(samplings6)] <- 0 #rewrite all NAs (missing values) as 0s
  samplings6 <- samplings6[samplings6$X2!=0,] #remove individuals who were observed only once
  
  
  
  #Calculate and store WIC/TNW for each population
  RIS1= import.RInSp(samplings1, col.header=TRUE, data.type="double")
  RIS2= import.RInSp(samplings2, col.header=TRUE, data.type="double")
  RIS3= import.RInSp(samplings3, col.header=TRUE, data.type="double")
  RIS4= import.RInSp(samplings4, col.header=TRUE, data.type="double")
  RIS5= import.RInSp(samplings5, col.header=TRUE, data.type="double")
  RIS6= import.RInSp(samplings6, col.header=TRUE, data.type="double")
  
  RIS1calc=WTcMC(RIS1, replicates = 999, weight = "N_items")
  RIS2calc=WTcMC(RIS2, replicates = 999, weight = "N_items")
  RIS3calc=WTcMC(RIS3, replicates = 999, weight = "N_items")
  RIS4calc=WTcMC(RIS4, replicates = 999, weight = "N_items")
  RIS5calc=WTcMC(RIS5, replicates = 999, weight = "N_items")
  RIS6calc=WTcMC(RIS6, replicates = 999, weight = "N_items")
  
  storeSIMS[z,1]=as.numeric(RIS1calc$WonT)
  storeSIMS[z,2]=as.numeric(RIS2calc$WonT)
  storeSIMS[z,3]=as.numeric(RIS3calc$WonT)
  storeSIMS[z,4]=as.numeric(RIS4calc$WonT)
  storeSIMS[z,5]=as.numeric(RIS5calc$WonT)
  storeSIMS[z,6]=as.numeric(RIS6calc$WonT)
  
}

storeSIMS=as.data.frame(storeSIMS)
colnames(storeSIMS)=c("July13C","Aug13C","Mar13C","July1B","Aug1B","Mar1B")
write.csv(storeSIMS,"PD_WIC_sims.csv")

#So now the data frame storeSIMS contains WIC/TNW for N simulated samples from EACH sampling period...
#...and we can compare differences to our empirical WIC/TNW values to evaluate whether our empirical...
#...values are more different than we would expect by chance (if the underlying niche parameters were...
#...in reality identical, as they were here in these simulations [i.e. all sampling sessions had...
#..the same true BIC and WIC])

hist(storeSIMS$July13C)
hist(storeSIMS$Aug13C)
hist(storeSIMS$Mar13C)
hist(storeSIMS$July1B)
hist(storeSIMS$Aug1B)
hist(storeSIMS$Mar1B)

#Here I'm checking the largest empirical differences in WIC/TNW between island-stage combos...
#First with a histrogram and then, if it looks close, actually counting to see if the empirical...
#...difference does indeed fall in the outer 5% of simulated differences

hist(storeSIMS$Mar13C-storeSIMS$July13C, breaks=50) #nope
hist(storeSIMS$Mar1B-storeSIMS$Aug1B, breaks=50) #maybe???...
sum((storeSIMS$Mar1B-storeSIMS$Aug1B)<=-.14) #nope
hist(storeSIMS$Aug13C-storeSIMS$Mar1B, breaks=50) #maybe???...
sum((storeSIMS$Aug13C-storeSIMS$Mar1B)>=.14) #close, but not quite


####################################################################
#################  SUPPLEMENTAL FIGS AND ANALYSES  #################
####################################################################

###
#Are any individual-level niche measurements correlated with...
#the number of observations per sampling stage?
###

plot(PHmeansPRE~sizePRE, data=indivs1)
summary(lm(PHmeansPRE~sizePRE, data=indivs1))
plot(PHmeansPOST~sizePOST, data=indivs1)
summary(lm(PHmeansPOST~sizePOST, data=indivs1))
plot(PHmeansDPOST~sizeDPOST, data=indivs1)
summary(lm(PHmeansDPOST~sizeDPOST, data=indivs1))

plot(PDmeansPRE~sizePRE, data=indivs1)
summary(lm(PDmeansPRE~sizePRE, data=indivs1))
plot(PDmeansPOST~sizePOST, data=indivs1)
summary(lm(PDmeansPOST~sizePOST, data=indivs1))
plot(PDmeansDPOST~sizeDPOST, data=indivs1)
summary(lm(PDmeansDPOST~sizeDPOST, data=indivs1))


plot(latPRE~sizePRE, data=indivs1)
summary(lm(latPRE~sizePRE, data=indivs1))
plot(latPOST~sizePOST, data=indivs1)
summary(lm(latPOST~sizePOST, data=indivs1))
plot(latDPOST~sizeDPOST, data=indivs1)
summary(lm(latDPOST~sizeDPOST, data=indivs1))



###
#MAKING MAPS WITH OBSERVATION COUNTS
###

#here are the Cartesian coordinates
maps<-read.csv(file="Herrmann_et_al_Maps.csv",header=TRUE)
colnames(maps)
colnames(maps)[1] <- "Site"
map13C=subset(maps,maps$Site=="SL13C")
map1B=subset(maps,maps$Site=="SL1B")

#Now I will assign these to each row of the main data frame in a single line...
#...because dplyr is DOPE
dat<-left_join(dat,maps,by=c("Site","Grid"))

#Now subset and plot
july13C=subset(dat,dat$Site=="SL13C" & dat$Session=="first")
aug13C=subset(dat,dat$Site=="SL13C" & dat$Session=="second")
mar13C=subset(dat,dat$Site=="SL13C" & dat$Session=="third")
july1B=subset(dat,dat$Site=="SL1B" & dat$Session=="first")
aug1B=subset(dat,dat$Site=="SL1B" & dat$Session=="second")
mar1B=subset(dat,dat$Site=="SL1B" & dat$Session=="third")

#add one fake observation for each grid square (needed for plotting)
july13CtoPLOT<-smartbind(july13C,map13C)
aug13CtoPLOT<-smartbind(aug13C,map13C)
mar13CtoPLOT<-smartbind(mar13C,map13C)
july1BtoPLOT<-smartbind(july1B,map1B)
aug1BtoPLOT<-smartbind(aug1B,map1B)
mar1BtoPLOT<-smartbind(mar1B,map1B)

rf <- colorRampPalette(brewer.pal(9,'YlGn'))
r <- rf(32)

#July on 13C

xBREAK=seq(-1,19,2)
yBREAK=seq(-1,19,2)

july13Chisto <- ggplot(july13CtoPLOT, aes(X,Y)) + stat_bin2d(breaks = list(x = xBREAK, y = yBREAK)) + 
  scale_fill_gradientn(colours=r,limits=c(0,27)) + 
  scale_x_continuous(limits=c(-1,19)) +
  scale_y_continuous(limits=c(-1,19)) +
  ggtitle("One-species island: Pre-removal")

# Get data - this includes counts and x,y coordinates 
newdat <- ggplot_build(july13Chisto)$data[[1]]
newdat$count=newdat$count-1

# add in text labels
july13Chisto + geom_text(data=newdat, aes((xmin + xmax)/2, (ymin + ymax)/2, 
                                          label=count), col="black")


#Aug on 13C
aug13Chisto <- ggplot(aug13CtoPLOT, aes(X,Y)) + stat_bin2d(breaks = list(x = xBREAK, y = yBREAK)) + 
  scale_fill_gradientn(colours=r, limits=c(0,27)) + 
  scale_x_continuous(limits=c(-1,19)) +
  scale_y_continuous(limits=c(-1,19)) +
  ggtitle("One-species island: Post-removal")

# Get data - this includes counts and x,y coordinates 
newdat <- ggplot_build(aug13Chisto)$data[[1]]
newdat$count=newdat$count-1

# add in text labels
aug13Chisto + geom_text(data=newdat, aes((xmin + xmax)/2, (ymin + ymax)/2, 
                                         label=count), col="black")

#March on 13C
mar13Chisto <- ggplot(mar13CtoPLOT, aes(X,Y)) + stat_bin2d(breaks = list(x = xBREAK, y = yBREAK)) + 
  scale_fill_gradientn(colours=r,limits=c(0,27)) + 
  scale_x_continuous(limits=c(-1,19)) +
  scale_y_continuous(limits=c(-1,19)) +
  ggtitle("One-species island: Delayed post-removal")

# Get data - this includes counts and x,y coordinates 
newdat <- ggplot_build(mar13Chisto)$data[[1]]
newdat$count=newdat$count-1

# add in text labels
mar13Chisto + geom_text(data=newdat, aes((xmin + xmax)/2, (ymin + ymax)/2, 
                                         label=count), col="black")


#July on 1B

xBREAK=seq(-1,17,2)
yBREAK=seq(-1,21,2)

july1Bhisto <- ggplot(july1BtoPLOT, aes(X,Y)) + stat_bin2d(breaks = list(x = xBREAK, y = yBREAK)) + 
  scale_fill_gradientn(colours=r,limits=c(0,27)) + 
  scale_x_continuous(limits=c(-1,17)) +
  scale_y_continuous(limits=c(-1,21)) +
  ggtitle("Two-species island: Pre-removal")

# Get data - this includes counts and x,y coordinates 
newdat <- ggplot_build(july1Bhisto)$data[[1]]
newdat$count=newdat$count-1

# add in text labels
july1Bhisto + geom_text(data=newdat, aes((xmin + xmax)/2, (ymin + ymax)/2, 
                                         label=count), col="black")


#Aug on 1B
aug1Bhisto <- ggplot(aug1BtoPLOT, aes(X,Y)) + stat_bin2d(breaks = list(x = xBREAK, y = yBREAK)) + 
  scale_fill_gradientn(colours=r, limits=c(0,27)) + 
  scale_x_continuous(limits=c(-1,17)) +
  scale_y_continuous(limits=c(-1,21)) +
  ggtitle("Two-species island: Post-removal")

# Get data - this includes counts and x,y coordinates 
newdat <- ggplot_build(aug1Bhisto)$data[[1]]
newdat$count=newdat$count-1

# add in text labels
aug1Bhisto + geom_text(data=newdat, aes((xmin + xmax)/2, (ymin + ymax)/2, 
                                        label=count), col="black")

#March on 1B
mar1Bhisto <- ggplot(mar1BtoPLOT, aes(X,Y)) + stat_bin2d(breaks = list(x = xBREAK, y = yBREAK)) + 
  scale_fill_gradientn(colours=r,limits=c(0,27)) + 
  scale_x_continuous(limits=c(-1,17)) +
  scale_y_continuous(limits=c(-1,21)) +
  ggtitle("Two-species island: Delayed post-removal")

# Get data - this includes counts and x,y coordinates 
newdat <- ggplot_build(mar1Bhisto)$data[[1]]
newdat$count=newdat$count-1

# add in text labels
mar1Bhisto + geom_text(data=newdat, aes((xmin + xmax)/2, (ymin + ymax)/2, 
                                        label=count), col="black")



###
#COMPARING THERMAL ENVIRONMENTS
###


#Compare the thermal environments across sites and sessions by comparing temperature models...
#that were left out at the sites.

#these are all model temps (both sites) from session 1 ONLY during the lizard observation window...
#which during Session 1 was 7:30-18:00
envTEMPS1<-read.csv(file="iBUTTS1cut.csv",header=TRUE)
colnames(envTEMPS1)[1]="Date"
str(envTEMPS1)

#Now to separate these iButtons among the two sites

#here is the location info for each iButton
modLOCATIONS1<-read.csv(file="modLOCATIONS1.csv",header=TRUE)
colnames(modLOCATIONS1)[1]="Site"

envTEMPS1<-envTEMPS1[,4:59] #removes dates and times from data frame

#vector to store temps
July13C<-as.numeric()
July1B<-as.numeric()


#Loop to pull the appropriate model temps from envTEMPS into the site-specific vectors
for (i in 1:nrow(modLOCATIONS1)){
  if(modLOCATIONS1$Site[i]=="SL13C"){
    July13C<-c(July13C,envTEMPS1[,i])
  }else{
    July1B<-c(July1B,envTEMPS1[,i])
  }
}

July13C<-July13C[!is.na(July13C)]
July1B<-July1B[!is.na(July1B)]


#Do the same for session 2 - only difference here is that the lizard observation window...
#was 7:30-15:30


envTEMPS2<-read.csv(file="iBUTTS2cut.csv",header=TRUE)
colnames(envTEMPS2)[1]="Date"
str(envTEMPS2)

#Now to separate these iButtons among the two sites

#here is the location info for each iButton
modLOCATIONS2<-read.csv(file="modLOCATIONS2.csv",header=TRUE)
colnames(modLOCATIONS2)[1]="Site"

envTEMPS2<-envTEMPS2[,4:59] #removes dates and times from data frame

#vector to store temps
Aug13C<-as.numeric()
Aug1B<-as.numeric()


#Loop to pull the appropriate model temps from envTEMPS into the site-specific vectors
for (i in 1:nrow(modLOCATIONS2)){
  if(modLOCATIONS2$Site[i]=="SL13C"){
    Aug13C<-c(Aug13C,envTEMPS2[,i])
  }else{
    Aug1B<-c(Aug1B,envTEMPS2[,i])
  }
}

Aug13C<-Aug13C[!is.na(Aug13C)]
Aug1B<-Aug1B[!is.na(Aug1B)]



#Do the same for session 3 - observation window was 8:30-17:00


envTEMPS3<-read.csv(file="iBUTTS3cut.csv",header=TRUE)
colnames(envTEMPS3)[1]="Date"
str(envTEMPS3)

#Now to separate these iButtons among the two sites

#here is the location info for each iButton
modLOCATIONS3<-read.csv(file="modLOCATIONS3.csv",header=TRUE)
colnames(modLOCATIONS3)[1]="Site"

envTEMPS3<-envTEMPS3[,4:59] #removes dates and times from data frame

#vector to store temps
Mar13C<-as.numeric()
Mar1B<-as.numeric()


#Loop to pull the appropriate model temps from envTEMPS into the site-specific vectors
for (i in 1:nrow(modLOCATIONS3)){
  if(modLOCATIONS3$Site[i]=="SL13C"){
    Mar13C<-c(Mar13C,envTEMPS3[,i])
  }else{
    Mar1B<-c(Mar1B,envTEMPS3[,i])
  }
}

Mar13C<-Mar13C[!is.na(Mar13C)]
Mar1B<-Mar1B[!is.na(Mar1B)]

#VIOLIN PLOT

x1=July13C
x2=July1B
x3=Aug13C
x4=Aug1B
x5=Mar13C
x6=Mar1B

catego1=c("Pre-removal","Post-removal","Delayed post-removal")
catego2=c("n=25,191","n=26,124","n=16,302","n=17,787","n=18,949","n=19,691")

plot(0:1,0:1,type="n",xlim=c(0.5,6.5),ylim=c(10,60),
     axes=FALSE,ann=FALSE)
vioplot(x1, x2, x3, x4,x5,x6,
        col=rep(c("forestgreen","saddlebrown"),2),
        rectCol = rep(c("black","forestgreen"),2),
        add=TRUE
)
axis(1, at=c(1.5,3.5,5.5),labels=catego1,cex.axis=1.4, tick=FALSE)
axis(2, cex.axis=1)
title(ylab="Model Temp (celsius)", xlab="",cex.lab=1.4)
text(x=1:6,y=58, labels=catego2, cex=1)


#COMPARE MEANS and VARIANCES

mean(July13C) #32.5
mean(July1B) #32.9
mean(Aug13C) #32.0
mean(Aug1B) #32.5
mean(Mar13C) #20.7
mean(Mar1B) #20.9

sd(July13C) #4.22
sd(July1B) #4.22
sd(Aug13C) #4.36
sd(Aug1B) #4.58
sd(Mar13C) #3.7
sd(Mar1B) #3.6

