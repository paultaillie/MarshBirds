# Format Kaleidoscope output for analysis in unmarked

# Paul J. Taillie
# North Carolina State University
# 1/11/18

# load packages
library(tidyverse)
library(stringr)
library(reshape2)

# #clear environment
remove(list=ls())


#load data
pass3.in<-read.csv("raw_data/pass3_results.csv")

# fix weird folder names
pass3<-pass3.in%>%
  mutate(SITE=as.character(FOLDER))%>%
  filter(SITE!="132")%>%  # remove accidently included 2016 data
  filter(SITE!="133")
for (i in 1:nrow(pass3)){
  if (nchar(pass3$SITE[i])>3){pass3$SITE[i]<-substr(pass3$SITE[i],1,3)}
}
unique(pass3$SITE)
         
# separate filename into date and time
pass3.date.time<-as.data.frame(str_split_fixed(pass3$IN.FILE,pattern="_",n=3))
colnames(pass3.date.time)<-c("ARU","date","time")
#convert date to julian
date1<-as.POSIXlt(as.character(pass3.date.time$date),format="%Y%m%d")
date2<-as.numeric(date1$yday)
pass3$julian_date<-date2
# convert time to min after midnight
 time<-pass3.date.time%>%
  mutate(time.hr=substr(as.character(time),1,2))%>%
  mutate(time.min=substr(as.character(time),3,4))%>%
  mutate(min_after_midnight=((as.numeric(time.hr)*60)+as.numeric(time.min)))%>%
  select(min_after_midnight)
pass3$time<-as.numeric(time[,1])

#  ----------  build observation covariate matrices  ----------------------
# date
obsCovs.date<-pass3%>%
  group_by(SITE)%>%
  mutate(visits=dense_rank(IN.FILE))%>%
  select(SITE,visits,julian_date)%>%
  dcast(SITE~visits,fun.aggregate=mean)
#time
obsCovs.time<-pass3%>%
  group_by(SITE)%>%
  mutate(visits=dense_rank(IN.FILE))%>%
  select(SITE,visits,time)%>%
  dcast(SITE~visits,fun.aggregate=mean)


# ----------  Set up y data (observations) for BLRA  ------------------------
# surveyed = 0 and no survey = NA
date.matrix<-as.matrix(obsCovs.date[,2:41])
index.temp=which(is.na(date.matrix)==F)
y.blra.aru=matrix(NA,23,40)
y.blra.aru[index.temp]<-0
# lump "BlRA" and "BLRA growl"
pass3$MANUAL.ID[which(pass3$MANUAL.ID=="BLRA growl")]<-"BLRA"
# BLRA detections
blra.det.hist<-pass3%>%
  mutate(blra=factor(MANUAL.ID))
levels(blra.det.hist$blra)<-list(
    blra="BLRA",
    other=levels(pass3$MANUAL.ID)[which(levels(pass3$MANUAL.ID)!="BLRA")])

blra.det.hist2<-blra.det.hist%>%
  mutate(blra2=as.numeric(case_when(blra=="blra"~1,blra=="other"~0)))%>%
  group_by(SITE)%>%
  mutate(visits=dense_rank(IN.FILE))%>%
  select(SITE,visits,blra2)%>%
  dcast(SITE~visits,fun.aggregate=max)

for (i in 1:nrow(blra.det.hist2)){
  for (j in 1:ncol(blra.det.hist2)){
    if(blra.det.hist2[i,j]<0){blra.det.hist2[i,j]<-NA}
  }
}


# ------       LEBI              ----------------
# ----------  Set up y data (observations) for LEBI  ------------------------
# surveyed = 0 and no survey = NA
date.matrix<-as.matrix(obsCovs.date[,2:41])
index.temp=which(is.na(date.matrix)==F)
y.lebi.aru=matrix(NA,23,40)
y.lebi.aru[index.temp]<-0

# LEBI detections
lebi.det.hist<-pass3%>%
  mutate(lebi=factor(MANUAL.ID))
levels(lebi.det.hist$lebi)<-list(
  lebi="LEBI",
  other=levels(pass3$MANUAL.ID)[which(levels(pass3$MANUAL.ID)!="LEBI")])

lebi.det.hist2<-lebi.det.hist%>%
  mutate(lebi2=as.numeric(case_when(lebi=="lebi"~1,lebi=="other"~0)))%>%
  group_by(SITE)%>%
  mutate(visits=dense_rank(IN.FILE))%>%
  select(SITE,visits,lebi2)%>%
  dcast(SITE~visits,fun.aggregate=max)

for (i in 1:nrow(lebi.det.hist2)){
  for (j in 1:ncol(lebi.det.hist2)){
    if(lebi.det.hist2[i,j]<0){lebi.det.hist2[i,j]<-NA}
  }
}

 # ------  write data to "processed_data" folder  ----------------
write.csv(obsCovs.date,"processed_data/date.aru.csv")
write.csv(obsCovs.time,"processed_data/time.aru.csv")
write.csv(blra.det.hist2,"processed_data/blra.aru.csv")
write.csv(lebi.det.hist2,"processed_data/lebi.aru.csv")









