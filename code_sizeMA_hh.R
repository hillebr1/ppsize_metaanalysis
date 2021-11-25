## Functional consequences of cell size
## Helmut Hillebrand, 24.11.2021

# remove remnants of previous analyses
rm(list=ls())
graphics.off()

# the following libraries are needed
library(calibrate)
library(maps)
library(Hmisc)
library(psych)
library(ggplot2)
library(grid)
library(gridExtra)
library(psych)
library(reshape2)
require(plyr)
library(agricolae)
library(ggExtra)
library(metafor)
library(cowplot)
library(ggridges)
library(multcomp)
library(RColorBrewer)
library(weights)
library(readr)
library(dplyr)
library(broom)
library(ggrepel)
library(ggalluvial)
library(tidyverse)
library(egg)

#set your working directory for this analysis
setwd("~/R/ppsize_metaanalysis")


## Section 1: alluvial plot of studies
studies <- read.csv("studies.csv", sep=";", dec=".")

# Select only studies that met criteria
studies<-studies[studies$Included==1,]
summary(studies)


#obtain number of studies per study system, type pf study, level of organisasation and aim
studies<-ddply(studies, .(System, Type,Level, Driver.Response, Aim), summarise, 
                  Freq= length(Year))

#create alluvial plot
alluvial<-ggplot(studies,
       aes(y = Freq, axis1 = System, axis2 = Level,axis3 = Type,   axis4 = Driver.Response,axis5 = Aim)) +
  geom_alluvium(aes(fill = as.factor(Aim))) +
  geom_stratum(fill = "grey", color = "black") +
  geom_label_repel(stat = "stratum", aes(label = ..stratum..)) +
  theme_minimal()+ylab("Number of studies")+
  theme(legend.position = "none")+
  theme(axis.title.y=element_text(size=16, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black"))+
  theme(axis.title.x=element_text(size=16,face="plain",colour="black"),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  scale_x_discrete(limits = c("System", "Level","Type", "Driver.Response","Aim"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") 

alluvial


## Section 2: Meta-analysis
data <- read.csv("data_incl.csv", sep=";", dec=".")
summary(data)

#delete NA
data<-data[!is.na(data$newresp),]

#a quick glance
Fig1 <-ggplot(data, aes(x=log10cellsize,y=newresp, col=caseID, shape =resptype))+
  geom_point(alpha=0.5)+
  theme_bw()+
  xlab(expression("Cell Size (LN µm³)"))+
  ylab(expression("Function"))+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=16,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=16,face="bold",colour="black"))+
  theme(legend.position="none")+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))+
  facet_wrap(~resp.subset, scales="free_y")
Fig1


#create effect sizes for the overall meta-analysis 
#slope and its standard error as effect size and sampling variance per unique caseID
results = data %>% 
  group_by(caseID, resp.categ, resp.scale,resp.subset) %>% # selects the cases to keep separate
  do(tidy(lm(newresp~log10cellsize, data = .))) %>% #does the linear regression
  filter(term == "log10cellsize") #extracts the slopes
summary(results)

#quick look at the effect size
bwplot(estimate~resp.categ|resp.scale,results)

#unique identifier per subset
results$USI<-do.call(paste, c(results[c("resp.categ", "resp.scale")], sep = "_"))
USI<-unique(results$USI)

#create an empty data frame
rma.df<-data.frame()

# the following loop cycles through all unique combinations of response category and response scale
#and provides mean effec size, significnce level (against 0) and confidence intervals
for(i in 1:length(USI)){
  temp<-results[results$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>1){#does the next step only if at least 2 data points are present
    lm1<-rma(estimate,std.error, data=temp)#makes the MA
    mean.ES <- coef(summary(lm1))[1, 1]#selects the mean effect size
    se.ES<- coef(summary(lm1))[1, 2]#selects its standard error
    p.ES<-coef(summary(lm1))[1, 4]#gives the p-value
    CIl.ES<-coef(summary(lm1))[1, 5]#gives lower CI
    CIu.ES<-coef(summary(lm1))[1, 6]#gives higher CI
    N<-dim(temp)[1]#gives the number of effect sizes
    rma.df<-rbind(rma.df,data.frame(temp[1,c("resp.categ", "resp.scale")],
                                    mean.ES,se.ES,p.ES,
                                    CIl.ES,CIu.ES,N))
    rm(temp)
  }
}

rma.df

#plot the effect sizes
# create a logicaö sequence for categories to plot and add a jitter
rma.df$resp.categ <- factor(rma.df$resp.categ,levels = c("growth", "element content", "resource uptake", "C fixation"))
pd<-position_dodge(.2)

MAresults<-ggplot(rma.df, aes(x=resp.categ,y=mean.ES,shape=resp.scale))+
  geom_hline(yintercept=0)+
  geom_point(size=5,position=pd)+
  scale_shape_manual(values=c(1,16))+
  geom_errorbar(aes(ymin=CIl.ES,ymax=CIu.ES),width=0,size=1,position=pd)+
  theme_bw()+
  xlab(expression("Response"))+
  ylab(expression("Mean effect size"))+#ylim(-2,4)+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=16,colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=16,colour="black",angle=0))+
  #theme(legend.position="none")+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))+coord_flip()

MAresults



# redo with separate slopes for phyla
names(data)
data$USI2<-do.call(paste, c(data[c("caseID","resp.categ", "resp.scale","resp.subset","Phylum")], sep = "_"))
USI2<-unique(data$USI2)

# here we opbtain a slope and its standard error as effect size per unique caseID AND phyla
# we only do so for phyla covering at least 3 data points and 1 order of magnitude in size

slope<-data.frame()
for(i in 1:length(USI2)){
  temp<-data[data$USI2==USI2[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does next step only if at least 3 datapoints are available
    if((max(temp$log10cellsize)-min(temp$log10cellsize))>1){#does next step only if at least 3 datapoints are available
      lm1<-lm(newresp~log10cellsize, temp)#makes a linear regreassion
      icpt <- coef(summary(lm1))[1, 1]#selects the slope
      slp <- coef(summary(lm1))[2, 1]#selects the slope
      se.slp<- coef(summary(lm1))[2, 2]#selects its standard error
      p<-anova(lm1)$'Pr(>F)'[1]#gives the p-value
      N<-dim(temp)[1] # gives the number of data points
      slope<-rbind(slope,data.frame(temp[1,c("caseID","resp.categ", "resp.scale","resp.subset","Phylum")],
                                              icpt,slp,se.slp,p,N))
      rm(temp)
    }
  }
}
summary(slope)


#quick look at the effect size
bwplot(slp~resp.categ|resp.scale,slope)

#create phylum specific effect sizes for plotting
#unique identifier per subset
slope$USI<-do.call(paste, c(slope[c("resp.categ", "resp.scale","Phylum")], sep = "_"))
USI<-unique(slope$USI)

#create an empty data frame

rma.df2<-data.frame()

# the following loop cycles through all unique cases

for(i in 1:length(USI)){
  temp<-slope[slope$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>1){#does the next step only if at least 2 data points are present
    lm1<-rma(slp,se.slp, data=temp)#makes the MA
    mean.ES <- coef(summary(lm1))[1, 1]#selects the mean effect size
    se.ES<- coef(summary(lm1))[1, 2]#selects its standard error
    p.ES<-coef(summary(lm1))[1, 4]#gives the p-value
    CIl.ES<-coef(summary(lm1))[1, 5]#gives lower CI
    CIu.ES<-coef(summary(lm1))[1, 6]#gives higher CI
    N<-dim(temp)[1]#gives the number of effect sizes
    rma.df2<-rbind(rma.df2,data.frame(temp[1,c("resp.categ", "resp.scale","Phylum")],
                                    mean.ES,se.ES,p.ES,
                                    CIl.ES,CIu.ES,N))
    rm(temp)
  }
}

rma.df2



# a way to plot the effect sizes
rma.df2$resp.categ <- factor(rma.df2$resp.categ,levels = c("growth", "element content", "resource uptake", "C fixation"))
pd2<-position_dodge(.5)

MAresults2<-ggplot(rma.df2, aes(x=resp.categ,y=mean.ES,col=Phylum, shape=resp.scale))+
  geom_hline(yintercept=0)+
  geom_point(size=5,position=pd2)+
  scale_shape_manual(values=c(1,16))+
  geom_errorbar(aes(ymin=CIl.ES,ymax=CIu.ES),width=0,size=1,position=pd2)+
  theme_bw()+
  xlab(expression("Response"))+
  ylab(expression("Mean effect size"))+#ylim(-2,4)+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=16,colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=16,colour="black",angle=0))+
  #theme(legend.position="none")+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))+coord_flip()


MAresults2


# formal meta-analysis testing for differences by phylogeny within each response category and scale
names(slope)

cfix.phylum<-rma(slp,se.slp**2,
                  mod=~resp.scale/Phylum,
                  data=slope[slope$resp.categ=="C fixation",])
summary(cfix.phylum)

uptake.phylum<-rma(slp,se.slp**2,
                 mod=~resp.scale/Phylum,
                 data=slope[slope$resp.categ=="resource uptake",])
summary(uptake.phylum)

elementcontent.phylum<-rma(slp,se.slp**2,
                 mod=~resp.scale/Phylum,
                 data=slope[slope$resp.categ=="element content",])
summary(elementcontent.phylum)

growth.phylum<-rma(slp,se.slp**2,
                 mod=~Phylum,
                 data=slope[slope$resp.categ=="growth",])
summary(growth.phylum)

# redo with separate slopes for system

names(data)
data$USI3<-do.call(paste, c(data[c("caseID","resp.categ", "resp.scale","resp.subset","system")], sep = "_"))
USI3<-unique(data$USI3)


# here we opbtain a slope and its standard error as effect size per unique caseID AND system

slope3<-data.frame()

for(i in 1:length(USI3)){
  temp<-data[data$USI3==USI3[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does next step only if at least 3 datapoints are available
    if((max(temp$log10cellsize)-min(temp$log10cellsize))>1){#does next step only if at least 3 datapoints are available
      lm1<-lm(newresp~log10cellsize, temp)#makes a linear regreassion
      icpt <- coef(summary(lm1))[1, 1]#selects the slope
      slp <- coef(summary(lm1))[2, 1]#selects the slope
      se.slp<- coef(summary(lm1))[2, 2]#selects its standard error
      p<-anova(lm1)$'Pr(>F)'[1]#gives the p-value
      N<-dim(temp)[1] # gives the number of data points
      slope3<-rbind(slope3,data.frame(temp[1,c("caseID","resp.categ", "resp.scale","resp.subset","system")],
                                    icpt,slp,se.slp,p,N))
      rm(temp)
    }
  }
}

summary(slope3)


#quick look at the effect size
bwplot(slp~resp.categ|resp.scale,slope3)


#doing the metaanalysis with metafor

#unique identifier per subset

slope3$USI<-do.call(paste, c(slope3[c("resp.categ", "resp.scale","system")], sep = "_"))
USI<-unique(slope3$USI)


#create an empty data frame

rma.df3<-data.frame()

# the following loop cycles through all unique cases

for(i in 1:length(USI)){
  temp<-slope3[slope3$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>0){#does the next step only if at least 2 data points are present
    lm1<-rma(slp,se.slp, data=temp)#makes the MA
    mean.ES <- coef(summary(lm1))[1, 1]#selects the mean effect size
    se.ES<- coef(summary(lm1))[1, 2]#selects its standard error
    p.ES<-coef(summary(lm1))[1, 4]#gives the p-value
    CIl.ES<-coef(summary(lm1))[1, 5]#gives lower CI
    CIu.ES<-coef(summary(lm1))[1, 6]#gives higher CI
    N<-dim(temp)[1]#gives the number of effect sizes
    rma.df3<-rbind(rma.df3,data.frame(temp[1,c("resp.categ", "resp.scale","system")],
                                      mean.ES,se.ES,p.ES,
                                      CIl.ES,CIu.ES,N))
    rm(temp)
  }
}

rma.df3

names(rma.df3)

# a way to plot the effect sizes
rma.df3$resp.categ <- factor(rma.df3$resp.categ,levels = c("growth", "element content", "resource uptake", "C fixation"))
pd2<-position_dodge(.5)

MAresults3<-ggplot(rma.df3, aes(x=resp.categ,y=mean.ES,col=system, shape=resp.scale))+
  geom_hline(yintercept=0)+
  geom_point(size=5,position=pd2)+
  scale_shape_manual(values=c(1,16))+
  geom_errorbar(aes(ymin=CIl.ES,ymax=CIu.ES),width=0,size=1,position=pd2)+
  theme_bw()+
  xlab(expression("Response"))+
  ylab(expression("Mean effect size"))+#ylim(-2,4)+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=16,colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=16,colour="black",angle=0))+
  #theme(legend.position="none")+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))+coord_flip()


MAresults3



# formal meta-analysis testing for differences by system
names(slope3)

cfix.system<-rma(slp,se.slp**2,
                 mod=~resp.scale/system,
                 data=slope3[slope3$resp.categ=="C fixation",])
summary(cfix.system)

uptake.system<-rma(slp,se.slp**2,
                   mod=~resp.scale/system,
                   data=slope3[slope3$resp.categ=="resource uptake",])
summary(uptake.system)

elementcontent.system<-rma(slp,se.slp**2,
                           mod=~resp.scale/system,
                           data=slope3[slope3$resp.categ=="element content",])
summary(elementcontent.system)

growth.system<-rma(slp,se.slp**2,
                   mod=~system,
                   data=slope3[slope3$resp.categ=="growth",])
summary(growth.system)




## Section 3: Raw data plots

#bring the p-value into the data for plotting
names(data)
names(results)
data<-merge(data,results[,c(1,2,3,4,9)],by=c("caseID","resp.categ","resp.scale","resp.subset"))
data$sign<-"yes"
data$sign[data$p.value<0.05]<-"no"
data$sign<-as.factor(data$sign)
summary(data)


#in order to arrange plots we create subset data

Cfixdata<-data[data$resp.categ=="C fixation",]
summary(Cfixdata)
unique(Cfixdata$resp.subset)
#two variable have very few data, in order to keep plot visible these are not plotted
Cfixdata<-Cfixdata[Cfixdata$resp.subset!="C affinity per Cavail",]
Cfixdata<-Cfixdata[Cfixdata$resp.subset!="C half saturation",]

#create logical arrangement of response variables
Cfixdata$facet<-factor(Cfixdata$resp.subset,levels=c('C fixation rate','respiration rate','exudation rate',
         'C affinity per light','Irradiance [compensation]','exudation fraction',
         'C fixation rate per C','specific respiration rate','ratio resp to ps'))
#replace unclear subset name (as of reviewer 2)
levels(Cfixdata$facet) <- sub("ratio resp to ps", "respir/photosyn", levels(Cfixdata$facet))
unique(Cfixdata$facet)

#plot
photsyn1 <-ggplot(Cfixdata, aes(x=log10cellsize,y=newresp, col=MA_ID,shape=resp.scale,linetype=sign))+
  geom_point(alpha=0.8)+
  theme_bw()+
  geom_smooth(method="lm", se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  xlab(expression("Cell size (log µm³)"))+
  ylab(expression("Response (log transformed)"))+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=16,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=16,face="bold",colour="black"))+
  theme(legend.position="none")+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))+
  facet_wrap(~facet,scales="free_y")
photsyn1
photsyn1<-tag_facet(photsyn1)
photsyn1<-photsyn1 + theme(strip.text = element_text(size=12,face="plain",colour="black"))




#same for resource uptake

Uptdata <-data[data$resp.categ=="resource uptake",]
unique(Uptdata$resp.subset)
Uptdata$facet<-factor(Uptdata$resp.subset,levels=c('N-spec N uptake','N uptake','N half saturation',
                                                   'N affinity','P-spec P uptake','P uptake',
                                                   'P half saturation'))


Upt1 <-ggplot(Uptdata, aes(x=log10cellsize,y=newresp, col=MA_ID,shape=resp.scale,linetype=sign))+
  geom_point(alpha=0.8)+
  theme_bw()+
  geom_smooth(method="lm", se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  xlab(expression("Cell size (Log µm³)"))+
  ylab(expression("Response (log transformed)"))+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=16,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=16,face="bold",colour="black"))+
  theme(legend.position="none")+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))+
  facet_wrap(~facet, scales="free_y", nrow=2)
Upt1
Upt1<-tag_facet(Upt1)
Upt1<-Upt1 + theme(strip.text = element_text(size=12,face="plain",colour="black"))
Upt1




# same for the two remaining categories

ECdata<-data[data$resp.categ=="element content"|
               data$resp.categ=="growth",]

unique(ECdata$resp.subset)
ECdata$facet<-factor(ECdata$resp.subset,levels=c('C content','N content','P content',
                                                   'Chla content','N storage','N:C ratio',
                                                   'Chla per biovolume ','growth rate',
                                                   'sedimentation rate'))

EC1 <-ggplot(ECdata, aes(x=log10cellsize,y=newresp, col=MA_ID,shape=resp.scale, linetype=sign))+
  geom_point(alpha=0.8)+
  theme_bw()+
  geom_smooth(method="lm", se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  xlab(expression("Cell size (log µm³)"))+
  ylab(expression("Response (log transformed)"))+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=16,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=16,face="bold",colour="black"))+
  theme(legend.position="none")+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))+
  facet_wrap(~facet, scales="free_y")
EC1
EC1<-tag_facet(EC1)
EC1<-EC1 + theme(strip.text = element_text(size=12,face="plain",colour="black"))
EC1


