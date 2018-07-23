# cov relationships


# COMPARE COVARIATE RELATIONSHIPS

# 

# Paul J. Taillie
# North Carolina State University
# 12/12/17



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
covs<-read.csv("processed_data/siteCovs.csv")
canopy.change<-read.csv("processed_data/dist_to_forest.csv")
covs.arcgis<-read.csv("processed_data/covs_ArcGIS.csv")


#add mean vegetation height change
#hist(canopy.change$Can_ht_change,breaks=15)
covs$CanopyChange<-rep(NA,nrow(covs))
covs$Dist_to_For<-rep(NA,nrow(covs))

for (i in 1:nrow(covs)){
  temp.index=which(canopy.change$PointID==covs$PointID[i])
  covs$CanopyChange[i]<-canopy.change$Can_ht_change[temp.index]
  covs$Dist_to_For[i]<-(canopy.change$Dist_to_Forest_ft[temp.index]*.3048)
  
}

#add fires
covs$num.fires<-rep(NA,nrow(covs))
covs.arcgis$PointYr=as.numeric(
  paste(as.character(covs.arcgis$Point),
        as.character(covs.arcgis$year),
        sep="."))
for (i in 1:nrow(covs)){
  temp.index=which(covs.arcgis$PointYr==covs$PointYr[i])
  covs$num.fires[i]<-covs.arcgis$total[temp.index]
}

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
  select(Herb_hits:Shrub.Woody,Dist_to_For)%>%  #select covariate columns
  mutate_all(funs(standardize))
siteCovs$year=factor(covs$Year) #add year back on as factor
siteCovs$fire=covs$num.fires  #add fire 
siteCovs$CanopyChange=covs$CanopyChange  #add canopy change as factor



ggplot(siteCovs)+
  geom_smooth(aes(x=Dist_to_For,y=Woody_hits),method="lm")


