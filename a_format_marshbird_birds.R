# marshbird birds and obsCovs

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

#load data
birds_2016<-read.csv(file="raw_data/marshbird_2016_birds_accessed 11_14_2016.csv") #point 54 deleted (P3)
birds_2017<-read.csv(file="raw_data/marshbird_2017_birds_9_1_17.csv") 

# add "year"
birds_2016$year=rep(2016,nrow(birds_2016))
birds_2017$year=rep(2017,nrow(birds_2017))
birds.all<-rbind(birds_2016,birds_2017)
for (i in 13:22){
  birds.all[which(is.na(birds.all[,i])==TRUE),i]<-0
}

# sort out species
levels(birds.all$Species)
birds.all$Species<-as.character(birds.all$Species)
birds.all$Species[which(birds.all$Species=="KI/CLRA")]<-"CLRA"
birds.all$Species[which(birds.all$Species=="KCLR")]<-"CLRA"
birds.all$Species[which(birds.all$Species=="KCRA")]<-"CLRA"
birds.all$Species[which(birds.all$Species=="KIRA")]<-"CLRA"
birds.all$Species[which(birds.all$Species=="CL/KIRA")]<-"CLRA"
birds.all$Species[which(birds.all$Species=="KCLRA")]<-"CLRA"

# convert date to julian date
date<-as.POSIXlt(as.character(birds.all$Date),format="%m/%d/%Y")
date2<-as.numeric(date$yday)
birds.all$Date<-date2
# convert time to minutes after sunrise
time<-as.character(birds.all$Time)
time.char=nchar(time)
index.time=which(time.char<5)
time[index.time]<-paste("0",time[index.time],sep="")
time2<-as.POSIXlt(time,format="%I:%M")
which(is.na(time2)==TRUE)
time3<-(60*as.numeric(time2$hour))+as.numeric(time2$min)
birds.all$Time<-time3

# new column for PointID*Year
birds.all$PointID<-as.character(birds.all$PointID)
for (i in 1:length(birds.all$PointID)){
  if(nchar(birds.all$PointID[i])==1){birds.all$PointID[i]<-paste("00",as.character(birds.all$PointID[i]),sep="")}
  if(nchar(birds.all$PointID[i])==2){birds.all$PointID[i]<-paste("0",as.character(birds.all$PointID[i]),sep="")}
  if(nchar(birds.all$PointID[i])==3){birds.all$PointID[i]<-as.character(birds.all$PointID[i])}
}
birds.all$PointYr<-paste(birds.all$PointID,as.character(birds.all$year),sep=".")
birds.all$PointYr<-as.character(birds.all$PointYr)
birds.all$Visit<-as.character(birds.all$Visit)

# summarize counts by passive vs. active
birds.all$passive=rep(0,nrow(birds.all))
birds.all$active=rep(0,nrow(birds.all))
for (i in 1:nrow(birds.all)){
  birds.all$passive[i]<-max(birds.all[i,13:17])
  birds.all$active[i]<-max(birds.all[i,18:22])
}

 # - ----------- function ----------------------------------------------
det.hist<-function(birds.all,spec,det.type){
  temp1<-as.data.frame.matrix(table(birds.all$PointYr,birds.all$Visit)) #determine visits that were surveyed
  
  det.hist<-data.frame(visit1=rep(NA,length(unique(birds.all$PointYr))),
                       visit2=rep(NA,length(unique(birds.all$PointYr))),
                       visit3=rep(NA,length(unique(birds.all$PointYr))),
                       visit4=rep(NA,length(unique(birds.all$PointYr))))
  rownames(det.hist)<-rownames(temp1)
  if (det.type=="passive"){index.type=26}
  if (det.type=="active"){index.type=27}
  
  for(i in 1:nrow(temp1)){
    for (j in 1:ncol(temp1)){
      if(temp1[i,j]>0){det.hist[i,j]<-0} #replace all surveyed visits with 0
    }} #end loops
  
  temp2<-birds.all[birds.all$Species==spec,] #pull out species of interest
  for (i in 1:nrow(det.hist)){
    temp3<-temp2[temp2$PointYr==rownames(det.hist)[i],] #pull out point of interest
    if(length(which(temp3$Visit==1))>0){det.hist[i,1]<-length(which(temp3[temp3$Visit==1,det.type]>0))}
    if(length(which(temp3$Visit==2))>0){det.hist[i,2]<-length(which(temp3[temp3$Visit==2,det.type]>0))}
    if(length(which(temp3$Visit==3))>0){det.hist[i,3]<-length(which(temp3[temp3$Visit==3,det.type]>0))}
    if(length(which(temp3$Visit==4))>0){det.hist[i,4]<-length(which(temp3[temp3$Visit==4,det.type]>0))}
  }
  if (spec=="SESP")
    {index.2016=which(substr(rownames(det.hist),5,8)=="2016")
    det.hist[index.2016,]<-NA}
    
  
  
  return(det.hist)
}# ---------------   end function    ---------------------------------




# ----------- Observation Covariate Matrices   -------------------------
matrix.blank<-matrix(NA,length(unique(birds.all$PointYr)),length(unique(birds.all$Visit)))
colnames(matrix.blank)<-unique(birds.all$Visit)
row.names(matrix.blank)<-unique(birds.all$PointYr)
date.blank<-matrix.blank #create blank matrices
time.blank<-matrix.blank
for(i in 1:nrow(birds.all)){
  time.blank[birds.all$PointYr[i], birds.all$Visit[i]]<-birds.all$Time[i]  #populate with date and time
  date.blank[birds.all$PointYr[i], birds.all$Visit[i]]<-birds.all$Date[i]}
time.1=cbind(time.blank,time.blank)  #stack side by side for passive and active
date.1=cbind(date.blank,date.blank)
time<-time.1[order(rownames(time.1)),] #re-order to match detection history
date<-date.1[order(rownames(date.1)),]

# Type (passive vs. active observation covariate)
matrix.type0<-matrix("passive",length(unique(birds.all$PointYr)),length(unique(birds.all$Visit)))
matrix.type1<-matrix("active",length(unique(birds.all$PointYr)),length(unique(birds.all$Visit)))
type=cbind(matrix.type0,matrix.type1)

# write obsCovs to csv's
#write.csv(date,file="processed_data/date.csv")
#write.csv(time,file="processed_data/time.csv")
#write.csv(type,file="processed_data/type.csv")

####     -----------    det hist  ----------------------------  ######
# 
bird="SESP"
y.passive<- det.hist(birds.all,bird,"passive")
y.active<- det.hist(birds.all,bird,"active")
y<-data.frame(y.passive,y.active)
colnames(y)<-c("1.passive",
                       "2.passive",
                       "3.passive",
                       "4.passive",
                       "1.active",
                       "2.active",
                       "3.active",
                       "4.active")
write.csv(y,file="processed_data/SESP.csv")






























