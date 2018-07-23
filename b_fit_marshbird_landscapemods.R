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
                              siteCovs=siteCovs)
umf.CLRA <- unmarkedFrameOccu(y=CLRA[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs)
umf.LEBI <- unmarkedFrameOccu(y=LEBI[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs)
umf.SESP <- unmarkedFrameOccu(y=SESP[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs)
umf.BLRA <- unmarkedFrameOccu(y=BLRA[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs)



