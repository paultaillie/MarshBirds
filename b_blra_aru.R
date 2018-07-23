# Fit models for BLRA and LEBI with in person surveys AND data from 
#     kaleidoscope analysis of ARU recordings

# Paul J. Taillie
# North Carolina State University
# 1/16/2018

# load packages
library(tidyverse)
library(unmarked)
library(MuMIn)
library(gridExtra)
# #clear environment
remove(list=ls())


# Read in "aru" data
time.aru<-read.csv("processed_data/time.aru.csv")
date.aru<-read.csv("processed_data/date.aru.csv")
blra.aru<-read.csv("processed_data/blra.aru.csv")
lebi.aru<-read.csv("processed_data/lebi.aru.csv")
type.aru<-data.frame(
  X=as.numeric(paste0(time.aru$SITE,".","2017")),
     matrix("ARU",23,40))


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


# ------   format covariates (copied from "b_fit_marshbird_modaverages.R") -----

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


#  combine in person and ARU
# blra
blra.aru2<-blra.aru%>%
  select(-X)%>%
  mutate(X=as.numeric(paste0(SITE,".","2017")))
blra.all=full_join(BLRA,blra.aru2,by="X")
blra.all=select(blra.all,-SITE)
# lebi
lebi.aru2<-lebi.aru%>%
  select(-X)%>%
  mutate(X=as.numeric(paste0(SITE,".","2017")))
lebi.all=full_join(LEBI,lebi.aru2,by="X")
lebi.all=select(lebi.all,-SITE)
#time
time.aru2<-time.aru%>%
  select(-X)%>%
  mutate(X=as.numeric(paste0(SITE,".","2017")))
time.all=full_join(time.in,time.aru2,by="X")
time.all=select(time.all,-SITE)
#date
date.aru2<-date.aru%>%
  select(-X)%>%
  mutate(X=as.numeric(paste0(SITE,".","2017")))
date.all=full_join(date.in,date.aru2,by="X")
date.all=select(date.all,-SITE)
#type
type.all=full_join(type.in,type.aru,by="X")

# Standardize Observation Covariates and store in list
obs.covs=list(
  time=standardize(as.matrix(time.all[,2:49])),
  date=standardize(as.matrix(date.all[,2:49])),
  type=as.matrix(type.all[,2:49])
)
obs.covs.in.person=list(
  time=standardize(as.matrix(time.all[,2:9])),
  date=standardize(as.matrix(date.all[,2:9])),
  type=as.matrix(type.all[,2:9])
)
#Standardize Site Covariates
siteCovs<- covs%>%
  select(Herb_hits:Shrub.Woody,Dist_to_For,num.fires)%>%  #select covariate columns
  mutate_all(funs(standardize))
siteCovs$year=factor(covs$Year) #add year back on as factor



#  -----------   Initial Modelling  -------------------

# Build unmarkedFrames
# blra
umf.blra.all <- unmarkedFrameOccu(y=blra.all[,2:49],
                         obsCovs=obs.covs,
                         siteCovs=siteCovs)
umf.blra.in.person <- unmarkedFrameOccu(y=blra.all[,2:9],
                                  obsCovs=obs.covs.in.person,
                                  siteCovs=siteCovs)
# lebi
umf.lebi.all <- unmarkedFrameOccu(y=lebi.all[,2:49],
                                  obsCovs=obs.covs,
                                  siteCovs=siteCovs)
umf.lebi.in.person <- unmarkedFrameOccu(y=lebi.all[,2:9],
                                        obsCovs=obs.covs.in.person,
                                        siteCovs=siteCovs)
#Fit null models to compare occupancy estimates
#blra
summary(blra.mod.all<-occu(~type+time+date+
                     I(time^2)+I(date^2)~year,
                   umf.blra.all))
summary(blra.mod.in.person<-occu(~type+time~year, umf.blra.in.person))
summary(blra.mod.in.person2<-occu(~type+date~year, umf.blra.in.person))

#lebi
summary(lebi.mod.all<-occu(~type+time+date+
                             I(time^2)+I(date^2)~year,
                           umf.lebi.all))
summary(lebi.mod.in.person<-occu(~type+time+date+
                             I(time^2)+I(date^2)~year, umf.lebi.in.person))


# -------------   ARU vs. in-person Ocupancy Probability plot     ------------

# collect results
# blra
   #all
    blra.aru.results=c(blra.mod.all@estimates@estimates$state@estimates[1],
                       confint(blra.mod.all,type="state")[1,])
    
   #in person
    blra.ip.results=c(blra.mod.in.person@estimates@estimates$state@estimates[1],
                       confint(blra.mod.in.person,type="state")[1,])
blra.results<-bind_rows(blra.aru.results,blra.ip.results)%>%
  rename(c("(Intercept)"="Mean","0.025"="Upper","0.975"="Lower"))%>%
  mutate(model=c("ARU","In-person"))%>%
  mutate(species=rep("BLRA",2))
# lebi
#all
lebi.aru.results=c(lebi.mod.all@estimates@estimates$state@estimates[1],
                   confint(lebi.mod.all,type="state")[1,])

#in person
lebi.ip.results=c(lebi.mod.in.person@estimates@estimates$state@estimates[1],
                  confint(lebi.mod.in.person,type="state")[1,])
lebi.results<-bind_rows(lebi.aru.results,lebi.ip.results)%>%
  rename(c("(Intercept)"="Mean","0.025"="Upper","0.975"="Lower"))%>%
  mutate(model=c("ARU","In-person"))%>%
  mutate(species=rep("LEBI",2))
plot1.dat=bind_rows(blra.results,lebi.results)
  
# build plot
ggplot(plot1.dat)+
  geom_pointrange(aes(x=species,
                      y=plogis(Mean),
                      ymin=plogis(Lower),
                      ymax=plogis(Upper),
                      shape=model),size=1,
                  position=position_dodge(width=.4))+
  ylim(c(0,1))+
  ylab("Occupancy Probability")

#  ----------   Detection covariate plots    --------------------------------

# calculate covariate ranges
time.min.ip<-min(time.all[,2:9],na.rm=T)
time.max.ip<-max(time.all[,2:9],na.rm=T)
time.min.aru<-min(time.all[,2:49],na.rm=T)
time.max.aru<-max(time.all[,2:49],na.rm=T)

date.min.ip<-min(date.all[,2:9],na.rm=T)
date.max.ip<-max(date.all[,2:9],na.rm=T)
date.min.aru<-min(date.all[,2:49],na.rm=T)
date.max.aru<-max(date.all[,2:49],na.rm=T)
time.pred.ip<-seq(time.min.ip,time.max.ip,length.out=20)
time.pred.aru<-seq(time.min.aru,time.max.aru,length.out=20)
date.pred.ip<-seq(date.min.ip,date.max.ip,length.out=20)
date.pred.aru<-seq(date.min.aru,date.max.aru,length.out=20)

# build "newdata"
newdata.ip.time<-data.frame(
  time=standardize(time.pred.ip),
  date=rep(0,20),
  type=factor(rep("active",20),levels=c("active","passive")))

newdata.aru.time<-data.frame(
  time=standardize(time.pred.aru),
  date=rep(0,20),
  type=factor(rep("active",20),levels=c("active","passive","ARU")))

newdata.ip.date<-data.frame(
  time=rep(0,20),
  date=standardize(date.pred.ip),
  type=factor(rep("active",20),levels=c("active","passive")))

newdata.aru.date<-data.frame(
  time=rep(0,20),
  date=standardize(date.pred.aru),
  type=factor(rep("active",20),levels=c("active","passive","ARU")))
# do predictions
blra.ip.pred.time<-as.data.frame(predict(blra.mod.in.person,
                       type="det",
                       newdata=data.frame(type=factor(rep("active",20),
                                                      levels=c("active","passive")),
                                          time=standardize(time.pred.ip)),
                                    
                                    se.fit=TRUE))
blra.ip.pred.date<-as.data.frame(predict(blra.mod.in.person2,
                                         type="det",
                                         newdata=data.frame(type=factor(rep("active",20),
                                                                        levels=c("active","passive")),
                                                            date=standardize(date.pred.ip)),
                                         
                                         se.fit=TRUE))
blra.aru.pred.time<-as.data.frame(predict(blra.mod.all,
                                    type="det",
                                    newdata=newdata.aru.time,
                                    se.fit=TRUE))
blra.aru.pred.date<-as.data.frame(predict(blra.mod.all,
                                          type="det",
                                          newdata=newdata.aru.date,
                                          se.fit=TRUE))

#put data together
blra.aru.plot.time1<-blra.aru.pred.time%>%
  mutate(time.x=time.pred.aru)%>%
  mutate(type=rep("ARU",20))
blra.aru.plot.time2<-blra.ip.pred.time%>%
  mutate(time.x=time.pred.ip)%>%
  mutate(type=rep("In-person",20))
blra.time.det=rbind(blra.aru.plot.time1,blra.aru.plot.time2)

blra.aru.plot.date1<-blra.aru.pred.date%>%
  mutate(date.x=date.pred.aru)%>%
  mutate(type=rep("ARU",20))
blra.aru.plot.date2<-blra.ip.pred.date%>%
  mutate(date.x=date.pred.ip)%>%
  mutate(type=rep("In-person",20))
blra.date.det=rbind(blra.aru.plot.date1,blra.aru.plot.date2)
# build plots
blra.time.plot<-ggplot(blra.time.det)+
  geom_ribbon(aes(x=(time.x/60),ymin=lower,ymax=upper,fill=type),alpha=.5)+
  scale_fill_grey()+
  geom_line(aes(x=time.x/60,y=Predicted,linetype=type),size=1.3)+
  #geom_line(aes(x=time.x/60,y=.474),linetype="dashed")+
  ylim(0,1)+xlim(0,24)+
  theme(legend.position="none")+
  ylab("Detection Probability")+
  xlab("Time of day")


blra.date.plot<-ggplot(blra.date.det)+
  geom_ribbon(aes(x=(date.x),ymin=lower,ymax=upper,fill=type),alpha=.5)+
  scale_fill_grey()+
  geom_line(aes(x=date.x,y=Predicted,linetype=type),size=1.3)+
  #geom_line(aes(x=time.x/60,y=.474),linetype="dashed")+
  ylim(0,1)+
  theme(legend.position="none")+
  ylab("Detection Probability")+
  xlab("Julian Date")
grid.arrange(blra.time.plot,blra.date.plot,nrow=2)
#original
blra.aru.plot<-ggplot(blra.aru.plot.dat)+
  geom_ribbon(aes(x=(time.x/60),ymin=lower,ymax=upper),fill="grey69")+
  geom_line(aes(x=time.x/60,y=Predicted),size=1.3)+
  geom_line(aes(x=time.x/60,y=.474),linetype="dashed")+
  ylim(0,1)+xlim(0,24)+
  theme(legend.position="none")+
  ylab("Detection Probability")+
  xlab("Time of day")



#  ------  Table 1  -  model selection  ------- ---------------

#fit aru models
blra.covs.all1<-occu(~type~year,umf.blra.all)
blra.covs.all2<-occu(~type+time~year,umf.blra.all)
blra.covs.all3<-occu(~type+time+date~year,umf.blra.all)
blra.covs.all4<-occu(~type+time+date+I(date^2)~year,umf.blra.all)
blra.covs.all5<-occu(~type+time+I(time^2)+date+I(date^2)~year,umf.blra.all)

blra.aru.mods<-fitList(
  blra.covs.all1,
  blra.covs.all2,
  blra.covs.all3,
  blra.covs.all4,
  blra.covs.all5)

blra.aru.results<-modSel(blra.aru.mods)@Full

#fit in person models
blra.covs.ip1<-occu(~type~year,umf.blra.in.person)
blra.covs.ip2<-occu(~type+time~year,umf.blra.in.person)
blra.covs.ip3<-occu(~type+time+date~year,umf.blra.in.person)
blra.covs.ip4<-occu(~type+time+date+I(date^2)~year,umf.blra.in.person)
blra.covs.ip5<-occu(~type+time+I(time^2)+date+I(date^2)~year,umf.blra.in.person)

blra.ip.mods<-fitList(
  blra.covs.ip1,
  blra.covs.ip2,
  blra.covs.ip3,
  blra.covs.ip4,
  blra.covs.ip5)

blra.ip.results<-modSel(blra.ip.mods)@Full
blra.ip.results2<-blra.ip.results%>%
  add_column("SEp(typeARU)" = rep(NA,5), .after = 8)%>%
  add_column("p(typeARU)" = rep(NA,5), .after = 8)
table1<-rbind(blra.ip.results2,blra.aru.results)

# write table 1 to .csv
#         write.csv(t(table1),"tables/param_estimates.csv")
#         write.csv(table1,"tables/modsel.csv")


#  -----------------    Table 2  - summary of sites with focal species  ----

blra.occ.sites<-data.frame(X=blra.all[,1],
                           occ=rowSums(blra.all[2:49],na.rm=T))
blra.occ.sites%>%
  mutate(X.character=as.character(X))%>%
  mutate(year=unlist(strsplit(X.character,".",fixed=T))[seq(2,length(X.character)*2,2)])%>%
  #group_by(year)%>%
  summarize(num.occ.sites=length(which(occ>0)))

blra.occ.sites.ip<-data.frame(X=blra.all[,1],
                           occ=rowSums(blra.all[2:9],na.rm=T))
blra.occ.sites.ip%>%
  mutate(X.character=as.character(X))%>%
  mutate(year=unlist(strsplit(X.character,".",fixed=T))[seq(2,length(X.character)*2,2)])%>%
  group_by(year)%>%
  summarize(num.occ.sites=length(which(occ>0)))



lebi.occ.sites<-data.frame(X=lebi.all[,1],
                           occ=rowSums(lebi.all[2:49],na.rm=T))
lebi.occ.sites%>%
  mutate(X.character=as.character(X))%>%
  mutate(year=unlist(strsplit(X.character,".",fixed=T))[seq(2,length(X.character)*2,2)])%>%
  group_by(year)%>%
  summarize(num.occ.sites=length(which(occ>0)))

###############################################################################
#  ---    Black Rail Model averaging   ---------------------------------------
###############################################################################

summary(Dist_to_For.mod   <-occu(~type~ year+Dist_to_For,umf.blra.all))
summary(Juncus.mod          <-occu(~type~ year+Juncus,umf.blra.all))
summary(Herb_hits.mod          <-occu(~type~ year+Herb_hits,umf.blra.all))
summary(year.mod          <-occu(~type~ year,umf.blra.all))
summary(Woody_height.mod          <-occu(~type~ year+Woody_height,umf.blra.all))
summary(num.fires.mod          <-occu(~type~ year+num.fires,umf.blra.all))
summary(cladium.mod          <-occu(~type~ year+Cladium,umf.blra.all))





global.mod<-occu(~type~Cladium+year+Dist_to_For+num.fires,umf.blra.all)

dredge.blra<-dredge(global.mod,rank=AIC, fixed=c("year"))
#get model averaged coefficients
top.mods<-get.models(dredge.blra,subset=delta<2)
avgm<-model.avg(top.mods)
results.out<-as.data.frame(avgm$coefArray[1,1:2,])

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

blra.results<-list(woody.plot,
            cladium.plot,
            juncus.plot,
            Herb.plot,
            Dist_to_For.plot,
            num.fires.plot,
            results.out))
















#  scrap   ---------------------------------------------------


#  other blra.full models
summary(blra.covs.all<-occu(~type+time+date+
                             I(time^2)+I(date^2)~year+Woody_hits+Cladium,
                           umf.blra.all))

summary(blra.covs.ip<-occu(~type+time+date+
                             I(time^2)+I(date^2)~year+Woody_hits+Cladium,
                                     umf.blra.in.person))



blra.occ.sites<-data.frame(X=blra.all[,1],
                           occ=rowSums(blra.all[2:49],na.rm=T))
blra.occ.sites%>%
  mutate(X.character=as.character(X))%>%
  mutate(year=unlist(strsplit(X.character,".",fixed=T))[seq(2,length(X.character)*2,2)])%>%
  group_by(year)%>%
  summarize(num.occ.sites=length(which(occ>0)))


