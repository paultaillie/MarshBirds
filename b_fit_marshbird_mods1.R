# fit models to marshbird data

# SUBJECT

# Paul J. Taillie
# North Carolina State University
# 12/12/17

# load packages
library(tidyverse)
library(unmarked)

# #clear environment
remove(list=ls())

# load packages
library(tidyverse)
library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(gridExtra)


# Read in data observation data
time.188<-read.csv("processed_data/time.csv")
date.188<-read.csv("processed_data/date.csv")
type.188<-read.csv("processed_data/type.csv")
BLRA<-read.csv("processed_data/BLRA.csv")
CLRA<-read.csv("processed_data/CLRA.csv")
LEBI<-read.csv("processed_data/LEBI.csv")
VIRA<-read.csv("processed_data/VIRA.csv")
SESP<-read.csv("processed_data/SESP.csv")
covs.lidar<-read.csv("processed_data/siteCovs_LiDAR.csv")
site.covs<-read.csv("processed_data/siteCovs.csv")


#Adjust LiDAR covs
lidar<-covs.lidar%>%
  mutate(PointYr=paste(as.character(covs.lidar$Point),as.character(covs.lidar$year),sep="."))%>%
  mutate(in.ydat=PointYr%in%site.covs$PointYr)%>%
  filter(in.ydat==TRUE)%>%
  mutate(PointYr=as.numeric(PointYr))%>%
  arrange(PointYr)

# Get rid of Cedar Island Points from field collected covariates
field<-site.covs[which(site.covs$PointYr%in%lidar$PointYr==T),]


#add field covariate to lidar
covs<-field%>%
  mutate(marsh.area=lidar$Marsh)%>%
  mutate(shrub.area=lidar$Shrub)%>%
  mutate(forest.area=lidar$Forest)%>%
  mutate(num.fires=lidar$total)%>%
  mutate(property=lidar$property)

# Get rid of Cedar Island surveys
VIRA<-VIRA[which(VIRA$X%in%covs$PointYr==T),]
BLRA<-BLRA[which(BLRA$X%in%covs$PointYr==T),]
CLRA<-CLRA[which(CLRA$X%in%covs$PointYr==T),]
LEBI<-LEBI[which(LEBI$X%in%covs$PointYr==T),]
SESP<-SESP[which(SESP$X%in%covs$PointYr==T),]
date.in<-date.188[which(date.188$X%in%covs$PointYr==T),]
time.in<-time.188[which(time.188$X%in%covs$PointYr==T),]
type.188$X<-date.188$X
type.in<-type.188[which(type.188$X%in%covs$PointYr==T),]


#standarize function
standardize<-function(xx){  
  xx.mean=mean(xx,na.rm=T)
  xx.sd=sd(xx,na.rm=T)
  xx.standardized<-(xx-xx.mean)/xx.sd
  return(xx.standardized)
}

#Standardize Site Covariates
siteCovs<- covs%>%
  select(Herb_hits:Shrub.Woody,marsh.area,shrub.area,forest.area)%>%  #select covariate columns
  mutate_all(funs(standardize))
siteCovs$year=factor(covs$Year) #add year back on as factor
siteCovs$fire=covs$num.fires  #add fire
siteCovs$property=factor(covs$property)

#make woody_hits categorical
siteCovs2<-siteCovs%>%
  mutate(Woody_cat=cut(Woody_hits,breaks=c(-1,-.4,8),labels=c("No Shrub","Shrub")))
# Standardize Observation Covariates and store in list
obs.covs=list(
  time=standardize(as.matrix(time.in[,2:9])),
  date=standardize(as.matrix(date.in[,2:9])),
  #type<-as.matrix(type.in%>%select(V1:V8)%>%mutate_all(funs(factor)))
  type=as.matrix(type.in[,2:9])
)

#Setup Unmarked Frame w/ categorical woody variable
umf.VIRA <- unmarkedFrameOccu(y=VIRA[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs2)
umf.CLRA <- unmarkedFrameOccu(y=CLRA[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs2)
umf.LEBI <- unmarkedFrameOccu(y=LEBI[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs2)
umf.SESP <- unmarkedFrameOccu(y=SESP[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs2)
umf.BLRA <- unmarkedFrameOccu(y=BLRA[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs)

#fit models
summary(mod1.vira<-occu(~type+time~Woody_cat+Cladium+Juncus,umf.VIRA))
summary(mod1.clra<-occu(~type+time~Woody_cat+Juncus+Cladium,umf.CLRA))
summary(mod1.lebi<-occu(~type+time~Woody_cat+Juncus+Cladium,umf.LEBI))
summary(mod1.sesp<-occu(~type+time~Woody_cat+Juncus+Cladium,umf.SESP))
summary(mod1.blra<-occu(~1~Cladium,umf.BLRA))



# set up prediction covariates
min.juncus=min(siteCovs2$Juncus)
max.juncus=max(siteCovs2$Juncus)
min.cladium=min(siteCovs2$Cladium)
max.cladium=max(siteCovs2$Cladium)
nd.1<-data.frame(Woody_cat=rep("No Shrub",10),
                 Juncus=seq(max.juncus,min.juncus,,10),
                 Cladium=seq(min.cladium,max.cladium,,10))
nd.2<-data.frame(Woody_cat=rep("Shrub",10),
                 Juncus=seq(max.juncus,min.juncus,,10),
                 Cladium=seq(min.cladium,max.cladium,,10))
nd<-rbind(nd.1,nd.2)

#do predictions
dat.vira<-predict(mod1.vira,type="state",newdata=nd,appendData=T)
dat.clra<-predict(mod1.clra,type="state",newdata=nd,appendData=T)
dat.lebi<-predict(mod1.lebi,type="state",newdata=nd,appendData=T)
dat.sesp<-predict(mod1.sesp,type="state",newdata=nd,appendData=T)

#make plots
clra.plot<-ggplot(dat.clra)+
  geom_line(aes(x=Cladium,y=Predicted,ymin=lower,ymax=upper,color=Woody_cat),size=1.3)+
  ggtitle("Clapper/King Rail")+
  ylim(0,1)+theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+xlab("Herbaceous Composition")
sesp.plot<-ggplot(dat.sesp)+
  geom_line(aes(x=Cladium,y=Predicted,ymin=lower,ymax=upper,color=Woody_cat),size=1.3)+
  ggtitle("Seaside Sparrow")+
  ylim(0,1)+theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+xlab("Herbaceous Composition")
vira.plot<-ggplot(dat.vira)+
  geom_line(aes(x=Cladium,y=Predicted,ymin=lower,ymax=upper,color=Woody_cat),size=1.3)+
  ggtitle("Virgina Rail")+
  ylim(0,1)+theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+xlab("Herbaceous Composition")
lebi.plot<-ggplot(dat.lebi)+
  geom_line(aes(x=Cladium,y=Predicted,ymin=lower,ymax=upper,color=Woody_cat),size=1.3)+
  ggtitle("Least Bittern")+
  ylim(0,1)+theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+xlab("Herbaceous Composition")

grid.arrange(clra.plot,sesp.plot,vira.plot,lebi.plot,nrow=2,ncol=2)



# ---------  fire plot    --------------------------------
siteCovs3<-siteCovs
siteCovs3$fire.cat<- case_when(
  siteCovs3$fire<0~"unburned",
  siteCovs3$fire>1~"burned")
umf.VIRA <- unmarkedFrameOccu(y=VIRA[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs3)
umf.CLRA <- unmarkedFrameOccu(y=CLRA[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs3)
umf.LEBI <- unmarkedFrameOccu(y=LEBI[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs3)
umf.SESP <- unmarkedFrameOccu(y=SESP[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs3)
summary(mod3.vira<-occu(~type~Juncus+Cladium+fire,umf.VIRA))
summary(mod3.clra<-occu(~type~Juncus+Cladium+fire,umf.CLRA))
summary(mod3.lebi<-occu(~type~Juncus+Cladium+fire,umf.LEBI))
summary(mod3.sesp<-occu(~type~Juncus+Cladium+fire,umf.SESP))

nd3<-data.frame(Woody_hits=rep(0,2),
                Juncus=rep(0,2),
                Cladium=rep(0,2),
                fire=factor(c(0,1)))
dat3.vira<-predict(mod3.vira,type="state",newdata=nd3,appendData=T)
dat3.clra<-predict(mod3.clra,type="state",newdata=nd3,appendData=T)
dat3.lebi<-predict(mod3.lebi,type="state",newdata=nd3,appendData=T)
dat3.sesp<-predict(mod3.sesp,type="state",newdata=nd3,appendData=T)

vira.fire<-ggplot(dat3.vira)+
  geom_pointrange(aes(x=fire,y=Predicted,ymin=lower,ymax=upper))+
  ggtitle("Virginia Rail")+
  ylim(0,1)+ylab("Occupancy Probability")+xlab(" ")
clra.fire<-ggplot(dat3.clra)+
  geom_pointrange(aes(x=fire,y=Predicted,ymin=lower,ymax=upper))+
  ggtitle("Clapper/King Rail")+
  ylim(0,1)+ylab(" ")+xlab(" ")
lebi.fire<-ggplot(dat3.lebi)+
  geom_pointrange(aes(x=fire,y=Predicted,ymin=lower,ymax=upper))+
  ggtitle("Least Bittern")+
  ylim(0,1)+ylab(" ")+xlab(" ")
sesp.fire<-ggplot(dat3.sesp)+
  geom_pointrange(aes(x=fire,y=Predicted,ymin=lower,ymax=upper))+
  ggtitle("Seaside Sparrow")+
  ylim(0,1)+ylab(" ")+xlab(" ")

grid.arrange(vira.fire,clra.fire,sesp.fire,lebi.fire,nrow=1,ncol=4)

#  ----------------- ------------- -----------------------------




summary(mod2.vira<-occu(~type+time~Phragmites,umf.VIRA))
summary(mod2.clra<-occu(~type+time~Phragmites,umf.CLRA))
summary(mod2.lebi<-occu(~type+time~Phragmites,umf.LEBI))
summary(mod2.sesp<-occu(~type+time~Phragmites,umf.SESP))





















