# marshbird covs

# SUBJECT

# Paul J. Taillie
# North Carolina State University
# DATE STARTED

# load packages
library(tidyverse)
library(unmarked)

# #clear environment
remove(list=ls())

# load packages
library(tidyverse)

#load data
covs.raw_2016<-read.csv(file="raw_data/marshbird_2016_veg_accessed 1_31_17.csv")
covs.raw_2017<-read.csv(file="raw_data/marshbird_2017_veg _9_1_17.csv")
covs.all<-rbind(covs.raw_2016,covs.raw_2017) #combine years

# extract relevant metrics
for (i in 7:25){
  covs.all[which(is.na(covs.all[,i])==TRUE),i]<-0
}


covs<-covs.all %>%
  group_by(PointID,Year)%>%
  summarize(Herb_hits=mean(Herb_hits),
            Herb_height=mean(Herb_height),
            Woody_hits=mean(Woody_hits),
            Woody_height=mean(Woody_height),
            Bare=mean(Bare),
            Wrack=mean(Wrack),
            Tree=mean(Tree),
            Snag=mean(Snag),
            Juncus=mean(Juncus),
            Cladium=mean(Cladium),
            Typha=mean(Typha),
            SpartinaPatens=mean(SpartinaPatens),
            Phragmites=mean(Phragmites),
            Scirpus=mean(Scirpus),
            Distichlis=mean(Distichlis),
            Andropogon=mean(Andropogon),
            Shrub.Woody=mean(Shrub.Woody))

covs$PointID<-as.character(covs$PointID)
for (i in 1:length(covs$PointID)){
  if(nchar(covs$PointID[i])==1){covs$PointID[i]<-paste("00",as.character(covs$PointID[i]),sep="")}
  if(nchar(covs$PointID[i])==2){covs$PointID[i]<-paste("0",as.character(covs$PointID[i]),sep="")}
  if(nchar(covs$PointID[i])==3){covs$PointID[i]<-as.character(covs$PointID[i])}
}
covs$PointYr<-paste(covs$PointID,as.character(covs$Year),sep=".")

temp<-c( "126.2016", "127.2016", "128.2016" ,"129.2016","134.2016" ,"135.2016")
temp1<-which(covs$PointYr%in%temp)
covs<-covs[-temp1,]
#write.csv(covs,file="siteCovs.csv")









