# 

# MARSHBIRD MODEL AVERAGING

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

# Standardize Observation Covariates and store in list
obs.covs=list(
  time=standardize(as.matrix(time.in[,2:9])),
  date=standardize(as.matrix(date.in[,2:9])),
  #type<-as.matrix(type.in%>%select(V1:V8)%>%mutate_all(funs(factor)))
  type=as.matrix(type.in[,2:9])
)

# table 1  - covariate ranges
stats<-covs%>%
  select(Dist_to_For, Juncus, Herb_hits, Woody_height, num.fires, Cladium,Year)%>%
  filter(Year==2016)%>%
  summarize_all(funs(min, mean, max))
table.1<-data.frame(min=t(stats[1,1:7]),
                    mean=t(stats[1,8:14]),
                    max=t(stats[1,15:21]))
#  write.csv(table.1,file="tables/cov_ranges.csv")


# correlation
cor(data.frame(
  Herb_hits=  covs$Herb_hits,
  Woody_height=  covs$Woody_height,
  Juncus=  covs$Juncus,
  Cladium=  covs$Cladium,
  Dist_to_For=covs$Dist_to_For
))

# -------------- function to use model averaged coefficients to predict response (occupancy)  -----------
plot.occ<-function(SPEC,obs.covs,covs){
  
#Standardize Site Covariates
siteCovs<- covs%>%
  select(Herb_hits:Shrub.Woody,Dist_to_For,num.fires)%>%  #select covariate columns
  mutate_all(funs(standardize))
siteCovs$year=factor(covs$Year) #add year back on as factor

#set up unmarked frame
umf <- unmarkedFrameOccu(y=SPEC[,2:9],
                              obsCovs=obs.covs,
                              siteCovs=siteCovs)

#fit global model and dredge
global.mod<-occu(~type+time+date+I(time^2)+I(date^2)~Herb_hits+Woody_height+Cladium+Juncus+year+Dist_to_For+num.fires,umf)
dredge.lebi<-dredge(global.mod,rank=AIC, subset=`psi(year)`)
#get model averaged coefficients
top.mods<-get.models(dredge.lebi,subset=delta<2)
avgm<-model.avg(top.mods)
results.out<-as.data.frame(avgm$coefArray[1,1:2,])
#set up prediction data
newdata.blank<-data.frame(
  Herb_hits=rep(0,20),
  Woody_height=rep(0,20),  
  Dist_to_For=rep(0,20),
  Cladium=rep(0,20),
  Juncus=rep(0,20),
  num.fires=rep(0,20),
  year=factor(rep("2017",20),levels=c("2016","2017")))
  
#woody
newdata.woody<-newdata.blank
Woody_height.X<-seq(min(covs$Woody_height),max(covs$Woody_height),length.out=20)
newdata.woody$Woody_height<-standardize(Woody_height.X)  
woody.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.woody,se.fit=TRUE))
woody.pred$Woody_height=Woody_height.X

woody.plot<-ggplot(woody.pred)+
  geom_ribbon(aes(x=Woody_height,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
  geom_line(aes(x=Woody_height,y=fit),size=1.3)+
  ylim(0,1)+
  theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+
  xlab("Woody vegetation height")

#Distance to forest
newdata.Dist_to_For<-newdata.blank
Dist_to_For.X<-seq(min(covs$Dist_to_For),max(covs$Dist_to_For),length.out=20)
newdata.Dist_to_For$Dist_to_For<-standardize(Dist_to_For.X)  
Dist_to_For.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.Dist_to_For,se.fit=TRUE))
Dist_to_For.pred$Dist_to_For=Dist_to_For.X

Dist_to_For.plot<-ggplot(Dist_to_For.pred)+
  geom_ribbon(aes(x=Dist_to_For,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
  geom_line(aes(x=Dist_to_For,y=fit),size=1.3)+
  ylim(0,1)+
  theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+
  xlab("Distance to Forest (m)")

#Herb hits
newdata.Herb<-newdata.blank
Herb_hits.X<-seq(min(covs$Herb_hits),max(covs$Herb_hits),length.out=20)
newdata.Herb$Herb_hits<-standardize(Herb_hits.X)  
Herb.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.Herb,se.fit=TRUE))
Herb.pred$Herb_hits=Herb_hits.X

Herb.plot<-ggplot(Herb.pred)+
  geom_ribbon(aes(x=Herb_hits,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
  geom_line(aes(x=Herb_hits,y=fit),size=1.3)+
  ylim(0,1)+
  theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+
  xlab("Herbaceous Vegetation (hits)")
#Cladium
newdata.cladium<-newdata.blank
Cladium.X<-seq(min(covs$Cladium),max(covs$Cladium),length.out=20)
newdata.cladium$Cladium<-standardize(Cladium.X)  
cladium.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.cladium,se.fit=TRUE))
cladium.pred$Cladium=Cladium.X

cladium.plot<-ggplot(cladium.pred)+
  geom_ribbon(aes(x=Cladium,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
  geom_line(aes(x=Cladium,y=fit),size=1.3)+
  ylim(0,1)+
  theme(axis.ticks=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+
  xlab("Cladium Cover")
#Juncus
newdata.juncus<-newdata.blank
Juncus.X<-seq(min(covs$Juncus),max(covs$Juncus),length.out=20)
newdata.juncus$Juncus<-standardize(Juncus.X)  
juncus.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.juncus,se.fit=TRUE))
juncus.pred$Juncus=Juncus.X

juncus.plot<-ggplot(juncus.pred)+
  geom_ribbon(aes(x=Juncus,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
  geom_line(aes(x=Juncus,y=fit),size=1.3)+
  ylim(0,1)+
  theme(axis.ticks=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+
  xlab("Juncus Cover")

#Fire
newdata.num.fires<-newdata.blank[1:5,]
num.fires.X<-seq(1,5)
newdata.num.fires$num.fires<-standardize(num.fires.X)  
num.fires.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.num.fires,se.fit=TRUE))
num.fires.pred$num.fires=num.fires.X

num.fires.plot<-ggplot(num.fires.pred)+
  geom_ribbon(aes(x=num.fires,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
  geom_line(aes(x=num.fires,y=fit),size=1.3)+
  ylim(0,1)+
  theme(axis.ticks=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+
  xlab("Number of Fires")

#Combined Gradient plot
newdata.combo<-newdata.blank
newdata.combo$Juncus      <-rev(standardize(Juncus.X))  
newdata.combo$Dist_to_For <-rev(standardize(Dist_to_For.X))
newdata.combo$Woody_height<-standardize(Woody_height.X) 
newdata.combo$Cladium    <-standardize(Cladium.X)
combo.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.combo,se.fit=TRUE))
combo.pred$gradient=seq(1:20)




return(list(woody.plot,
cladium.plot,
juncus.plot,
Herb.plot,
Dist_to_For.plot,
num.fires.plot,
combo.pred,
results.out))
}############  --------  end function   -------------------------------

#run function on each species
SESP.results<-plot.occ(SESP,obs.covs,covs)
CLRA.results<-plot.occ(CLRA,obs.covs,covs)
LEBI.results<-plot.occ(LEBI,obs.covs,covs)
VIRA.results<-plot.occ(VIRA,obs.covs,covs)


# ------------------------------- TAble 2 data  (parameter estimates)-------------------------------
params<-bind_rows(
  SESP.results[[8]],
  CLRA.results[[8]],
  LEBI.results[[8]],
  VIRA.results[[8]])
params.out<-params%>%
  mutate(spec=c("SESP","SESP","CLRA","CLRA","LEBI","LEBI","VIRA","VIRA"))%>%
  mutate(quantity=c("mean","se","mean","se","mean","se","mean","se"))


#   write.csv(params.out,file="processed_data/param_estimates2.csv")


#  ------------------------  Figure - Shore to Forest Gradient  ---------------------------------

# add species
SESP.gradient<-mutate(SESP.results[[7]],species="Seaside Sparrow")
CLRA.gradient<-mutate(CLRA.results[[7]],species="Clapper/King Rail")
LEBI.gradient<-mutate(LEBI.results[[7]],species="Least Bittern")
VIRA.gradient<-mutate(VIRA.results[[7]],species="Virgina Rail")

# combine data
gradient.dat.all<-bind_rows(
  SESP.gradient,
  CLRA.gradient,
  LEBI.gradient,
  VIRA.gradient)


combo.plot<-ggplot(gradient.dat.all)+
  geom_ribbon(aes(x=gradient,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
  geom_line(aes(x=gradient,y=fit),size=1.3)+
  facet_grid(species~.)+
  ylim(0,1)+
  theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
  ylab("Occupancy Probability")+
  xlab("Shoreline-to-forest gradient")


