# Preliminary analysis of Kaleidoscope results

# Paul J. Taillie
# North Carolina State University
# 1/9/18

# load packages
library(tidyverse)

# #clear environment
remove(list=ls())


#load data
pass1<-read.csv("raw_data/pass1_results.csv")
pass2<-read.csv("raw_data/pass2_results.csv")
pass3.in<-read.csv("raw_data/pass3_results.csv")



# -------------    Pass 1       -----------------------
# change other rails to "non-focal"
pass1.vocs.by.spec<-pass1%>%
  mutate(spec=as.character(MANUAL.ID))%>%
  mutate(spec=case_when(
    pass1$MANUAL.ID == "CLRA grunt" ~ "non focal",
    pass1$MANUAL.ID == "CLRA kek" ~ "non focal",
    pass1$MANUAL.ID == "focal mix" ~ "non focal",
    pass1$MANUAL.ID == "LEBI" ~ "LEBI",
    pass1$MANUAL.ID == "BLRA" ~ "BLRA",
    pass1$MANUAL.ID == "BLRA growl" ~ "BLRA growl",
    pass1$MANUAL.ID == "non focal" ~ "non focal",
    pass1$MANUAL.ID == "SESP" ~ "non focal",
    pass1$MANUAL.ID == "SORA" ~ "non focal",
    pass1$MANUAL.ID == "VIRA call" ~ "non focal",
    pass1$MANUAL.ID == "VIRA grunt" ~ "non focal",
    pass1$MANUAL.ID == "" ~ ""
  ))%>%
  group_by(spec)%>%
  summarize(count=n())%>%
  mutate(percent=count/66.84)%>%
  mutate(percent=round(percent,2))
pass1.vocs.by.spec # table summarizing labeled vocalizations by species
write.csv(pass1.vocs.by.spec,file="tables/pass1/pass1_vocs_by_spec.csv")

#extract rows from pass1 that were labeled
pass1.labeled<-pass1[pass1$MANUAL.ID!="",]
#number of labeled clusters (all of them - 77) 
length(unique(pass1.labeled$TOP3MATCH)) 
#table of number of labeled vocalizations per cluster
labels.per.cluster<-as.data.frame(table(pass1.labeled$TOP3MATCH))
write.csv(labels.per.cluster,file="tables/pass1/labels_per_cluster.csv")


###   -------------  Pass 2   ----------------------------

# number of vocalizations
nrow(pass2)
# number of BLRA vocalizations
nrow(pass2[pass2$TOP1MATCH=="BLRA",])
# number of LEBI vocalizations
nrow(pass2[pass2$TOP1MATCH=="LEBI",])

#extract BLRA from pass 2
pass2.blra<-pass2[pass2$TOP1MATCH=="BLRA",]
# summarize revisions to those vocalizations called BLRA by kaleidoscope 
#    in pass 2
pass2.blra.revisions<-as.data.frame(table(pass2.blra$MANUAL.ID))
write.csv(pass2.blra.revisions,file="tables/pass2/blra_revisions.csv")

#extract LEBI from pass 2
pass2.lebi<-pass2[pass2$TOP1MATCH=="LEBI",]
# summarize revisions to those vocalizations called BLRA by kaleidoscope 
#    in pass 2
pass2.lebi.revisions<-as.data.frame(table(pass2.lebi$MANUAL.ID))
write.csv(pass2.lebi.revisions,file="tables/pass2/lebi_revisions.csv")



#  -------   Pass 3   -----------------------------------------------------
# remove data accidentally added from 2016
# remove data  accidentally included from 2016
pass3<-pass3.in%>%
  mutate(SITE=as.character(FOLDER))%>%
  filter(SITE!="132")%>%
  filter(SITE!="133")


# Number of BLRA detected
pass3.blra<-pass3[pass3$TOP1MATCH.=="BLRA",]
nrow(pass3.blra)
pass3.blra.g<-pass3[pass3$TOP1MATCH.=="BLRA growl",]
nrow(pass3.blra.g)
# Number of LEBI detected
pass3.lebi<-pass3[pass3$TOP1MATCH.=="LEBI",]
nrow(pass3.lebi)

# summary of BLRA
pass3.blra.confirm<-as.data.frame(table(pass3.blra$MANUAL.ID))
pass3.blra.confirm
pass3.blra.g.confirm<-as.data.frame(table(pass3.blra.g$MANUAL.ID))
pass3.blra.g.confirm

# summary of LEBI
pass3.lebi.confirm<-as.data.frame(table(pass3.lebi$MANUAL.ID))
pass3.lebi.confirm

#sites with true BLRA
pass3.true.blra<-pass3.blra[pass3.blra$MANUAL.ID=="BLRA"|
                              pass3.blra$MANUAL.ID=="BLRA growl",]
blra.sites<-as.data.frame(table(pass3.true.blra$FOLDER))
write.csv(blra.sites,file="tables/pass3/blra_sites.csv")

#sites with true LEBI
pass3.true.lebi<-pass3.lebi[pass3.lebi$MANUAL.ID=="LEBI",]
lebi.sites<-as.data.frame(table(pass3.true.lebi$FOLDER))
write.csv(lebi.sites,file="tables/pass3/lebi_sites.csv")
















