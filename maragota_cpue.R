##########################################################
# Analysis of Labrus bergylta catch rates in Galicia#
# using data sampled onboard fishing vessels by the UTPB #
##########################################################

#Date start: April 2014 by Jaime Otero (bulk of code from pollachius analysis)
#Last update: 19-02-2015 by Alex Alonso
#many code string were carried out by Jaime Otero during POollachius pollachius analysis
#input data comes from "preparation_data_base_casgass.R" code that generate faneca.txt file
#General objective: to analyse variation of catch rates from monitoring program data of UTPB
#
# Project: CASGASS
#
# Authors: Jaime Otero (jotero@iim.csic.es) and Alex Alonso (alex@iim.csic.es)
# Further modifications by: Alexandre Alonso (alex@iim.csic.es)
# Institute of Marine Research (Vigo)

#libraries

library(pscl)
library(mgcv)
library(MASS)
library(splines)
library(effects)
library(ggplot2)
library(boot)
library(snow)
library(lmtest)
library(sp)
library(gstat)
library(glmmADMB) 
library(car)
library(bbmle)
#install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R",type="source")
library(coefplot2)
#Joshep M. Hilbe libraries for count data analysis
library(msme)
library(COUNT)

source(file="C:\\Users\\alex\\Dropbox\\casgass\\HighStatLib.r")
source(file="C:\\Users\\alex\\Dropbox\\casgass\\multiple_ggplot.r")
source(file="C:\\Users\\alex\\Dropbox\\casgass\\CI_mean.r")

#set working directory for input and output
setwd("D:\\iim-csic\\proyectos\\ICES_CASGASS\\analysis\\wrasse")

#to start since the last point
load("wrasse.RData")

#1# Load species data
#el código está hecho en base al análisis de pollachius
#adaptar el código a cada especies
wrasse<-read.table(file="C:\\Users\\alex\\Dropbox\\casgass\\pinto.txt",header=T,dec=".",sep=",")

head(wrasse)
dim(wrasse)
length(unique(wrasse$Idflota)) # Number of vessels surveyed


# ----------------------------- #
#2# Evaluate fishing per gear, distribution of response and other elementary information

pinto.tot<-sum(wrasse$Ntot,na.rm=T) # Total number of caught pouting
pintos<-as.numeric(tapply(wrasse$Ntot,wrasse$Gear,sum,na.rm=T)) # Total number per gear
gear.hauls<-tapply(wrasse$Ntot,wrasse$Gear,length) # Hauls per gear
fish<-ifelse(wrasse$Ntot==0,0,1) # Variable to classify 0/1 hauls
zero.fish<-gear.hauls-tapply(fish,wrasse$Gear,sum,na.rm=T) # Zeros per gear
gear.zeros<-round((zero.fish*100)/gear.hauls,2) # Percentage of zeros per gear

basic.info<-data.frame(cbind(gear.hauls,gear.zeros,(pintos*100)/pinto.tot))
colnames(basic.info)<-c("hauls","zeros","catch")
basic.info$gear<-c("Boliche","Bou de Man","Bou de Vara","Cerco Piobardeira","Liña cordel",
                   "Miños","Nasa choco","Nasa nécora",
                   "Nasa peixe","Nasa pulpo","Voitron","Palangrillo","Racu","Raeira",
                   "Rastro camaron","Trasmallo","Vetas")

pdf(file="gearsWrasse.pdf",width=15,height=8)

ggplot(data=basic.info,aes(x=gear,y=catch))+
  geom_bar(stat="identity")+
  scale_y_continuous(limits=c(0,60),"Frequency")+
  scale_x_discrete("Gear")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=14),
        axis.text.y=element_text(size=14))+
  geom_text(aes(label=paste("[",hauls,",",zeros,"]"),vjust=-0.5),size=3)

dev.off()

par(mar=c(5,5,3,3))
plot(table(wrasse$Ntot),ylab="Frequency",xlab="Number of wrasse caught per haul")
par(mar=c(5,5,3,3))
plot(table(wrasse$Ntot),ylab="Frequency",xlab="Number of wrasse caught per haul",xlim=c(0,10))

basic.info

basic.info[order(-basic.info$catch), ]
basic.info[order(basic.info$zeros), ]
#seleccionamos artes de enmalle: vetas y miños y trasmallos en base
#al porcentage de capturas muestreadas
#entre las dos representan casi el 70% y tienen porcentages de zeros <50%

# ----------------------------- #
#3# Select data set that only include the 2 main gears
wrasse1<-wrasse[wrasse$Gear=="TRASMALLOS" | wrasse$Gear=="VETAS"| wrasse$Gear=="MINOS",]

head(wrasse1)
dim(wrasse1)
length(unique(wrasse1$Idflota))
length(unique(wrasse1$Idjornada))
length(unique(wrasse1$Idlance))

# ----------------------------- #
#4# Add variables and recode as factors when needed

#4.1# Variables already in the data

summary(wrasse1)

wrasse1<-wrasse1[-which(is.na(wrasse1$Lat)),]
wrasse1<-wrasse1[-which(is.na(wrasse1$Seafloor)),]
wrasse1<-wrasse1[-which(is.na(wrasse1$Depth)),]

dim(wrasse1)
str(wrasse1)
summary(wrasse1)

wrasse1$fyear<-factor(wrasse1$Year)
wrasse1$fgear<-factor(wrasse1$Gear)
wrasse1$fcrew<-factor(ifelse(wrasse1$Crew <= 3,1,2)) # Crew as factor
wrasse1$fZoneA<-factor(wrasse1$ZoneA)
wrasse1$fZoneO<-factor(wrasse1$ZoneO)
wrasse1$Observer<-factor(wrasse1$Observer)

#4.2# Gear size

#4.2.1# Recorded gear size in each trip

gearSize<-read.csv2(file="C:\\Users\\alex\\Dropbox\\casgass\\dimensiones_artes_pesca.csv",header=T,dec=".",sep=";")

head(gearSize)

gearSize1<-gearSize[gearSize$ARTE=="TRASMALLOS",] # Only nasa-peixe
gearSize2<-gearSize[gearSize$ARTE=="VETAS",] # Only Vetas
gearSize3<-gearSize[gearSize$ARTE=="MINOS",] # Only Vetas
gearSize1<-gearSize1[,c(1:3,5)]
gearSize2<-gearSize2[,c(1:3,5)]
gearSize3<-gearSize3[,c(1:3,5)]
head(gearSize1)
head(gearSize2)
head(gearSize3)

dim(gearSize1)
dim(gearSize2)
dim(gearSize3)

#trasmallo
par(mfrow=c(1,2))
summary(gearSize1$longitud);hist(gearSize1$longitud) # Looks OK, units m (I guess)
summary(gearSize1$altura);hist(gearSize1$altura) # Strange values
summary(gearSize1$altura);hist(gearSize1$altura,breaks=100,xlim=c(0,100))
# Remove mistakes and scale in m
gearSize1$altura<-ifelse(gearSize1$altura<5,gearSize1$altura,
                         ifelse(gearSize1$altura>5&gearSize1$altura<40,gearSize1$altura/10,
                                ifelse(gearSize1$altura>1000,gearSize1$altura/1000,gearSize1$altura/100)))
summary(gearSize1$altura);hist(gearSize1$altura) #it looks better, but check altura less than 1m

#veta
par(mfrow=c(2,2))
summary(gearSize2$longitud);hist(gearSize2$longitud,breaks=50) # Looks OK
summary(gearSize2$altura);hist(gearSize2$altura,breaks=50) # Strange values
summary(gearSize2$altura);hist(gearSize2$altura,breaks=50,xlim=c(0,350)) # Strange values
which(gearSize2$altura<100)
gearSize2[which(gearSize2$altura<100),]
summary(gearSize2[gearSize2$altura<100,])
# Remove mistakes and scale in m
gearSize2$altura<-ifelse(gearSize2$altura<5,gearSize2$altura,
                         ifelse(gearSize2$altura>5&gearSize2$altura<100,gearSize2$altura/10,gearSize2$altura/100)) 
summary(gearSize2$altura);hist(gearSize2$altura,breaks=50) #it looks much better

#miño
par(mfrow=c(1,2))
summary(gearSize3$longitud);hist(gearSize3$longitud,breaks=50) # Looks OK
summary(gearSize3$altura);hist(gearSize3$altura,breaks=50) # Strange values
hist(gearSize3$altura,breaks=50,xlim=c(1000,3000))
hist(gearSize3$altura,breaks=50,xlim=c(0,100))
which(gearSize3$altura<100)
summary(gearSize3[gearSize3$altura<125,])
# Remove mistakes and scale in m
gearSize3$altura<-ifelse(gearSize3$altura==2500,2.5,ifelse(gearSize3$altura < 5,gearSize3$altura,gearSize3$altura/100))
summary(gearSize3$altura);hist(gearSize3$altura,breaks=50) #it looks much better

#join them 
gearSize4<-data.frame(rbind(gearSize1,gearSize2,gearSize3))
head(gearSize4)
tail(gearSize4)
dim(gearSize4)

#adding gear size information to main data set
wrasse1<-merge(wrasse1,gearSize4,by.x=c("Idlance","Gear"),by.y=c("Idlance","ARTE")) # Merge datasets
head(wrasse1)

#treatment of gearsize by Gear independently
wrasse2<-wrasse1[wrasse1$Gear=="TRASMALLOS",]
wrasse3<-wrasse1[wrasse1$Gear=="VETAS",]
wrasse4<-wrasse1[wrasse1$Gear=="MINOS",]

dim(wrasse2)
head(wrasse2)
dim(wrasse3)
head(wrasse3)
dim(wrasse4)
head(wrasse4)

# Fill in NAs with vessel-average data
#trasmallo
ids<-unique(wrasse2$Idflota) 
long<-length(ids)
lista<-as.list(1:long)
names(lista)<-ids

for (i in 1:long) {lista[[i]]<-wrasse2[wrasse2$Idflota==ids[i],]}
for (i in 1:long) {lista[[i]]$longitud2<-ifelse(sum(lista[[i]]$longitud,na.rm=T) > 0,mean(lista[[i]]$longitud,na.rm=T),NA)}
for (i in 1:long) {lista[[i]]$altura2<-ifelse(sum(lista[[i]]$altura,na.rm=T) > 0,mean(lista[[i]]$altura,na.rm=T),NA)}

wrasse2<-do.call(rbind,lista)
head(wrasse2)
summary(wrasse2)

wrasse2$longitud3<-ifelse(is.na(wrasse2$longitud),wrasse2$longitud2,wrasse2$longitud)
wrasse2$altura3<-ifelse(is.na(wrasse2$altura),wrasse2$altura2,wrasse2$altura)
head(wrasse2)
tail(wrasse2)
summary(wrasse2)

wrasse2$longitud4<-ifelse(is.na(wrasse2$longitud3),
                           mean(wrasse2$longitud3,na.rm=T),wrasse2$longitud3) # Fill in with fleet average data for those vessels without any information
wrasse2$altura4<-ifelse(is.na(wrasse2$altura3),
                         mean(wrasse2$altura3,na.rm=T),wrasse2$altura3)
head(wrasse2)
tail(wrasse2)
summary(wrasse2)

#vetas
ids<-unique(wrasse3$Idflota) 
long<-length(ids)
lista<-as.list(1:long)
names(lista)<-ids

for (i in 1:long) {lista[[i]]<-wrasse3[wrasse3$Idflota==ids[i],]}
for (i in 1:long) {lista[[i]]$longitud2<-ifelse(sum(lista[[i]]$longitud,na.rm=T) > 0,mean(lista[[i]]$longitud,na.rm=T),NA)}
for (i in 1:long) {lista[[i]]$altura2<-ifelse(sum(lista[[i]]$altura,na.rm=T) > 0,mean(lista[[i]]$altura,na.rm=T),NA)}

wrasse3<-do.call(rbind,lista)
head(wrasse3)
summary(wrasse3)

wrasse3$longitud3<-ifelse(is.na(wrasse3$longitud),wrasse3$longitud2,wrasse3$longitud)
wrasse3$altura3<-ifelse(is.na(wrasse3$altura),wrasse3$altura2,wrasse3$altura)
head(wrasse3)
tail(wrasse3)
summary(wrasse3)

wrasse3$longitud4<-ifelse(is.na(wrasse3$longitud3),
                           mean(wrasse3$longitud3,na.rm=T),wrasse3$longitud3) # Fill in with fleet average data for those vessels without any information
wrasse3$altura4<-ifelse(is.na(wrasse3$altura3),
                         mean(wrasse3$altura3,na.rm=T),wrasse3$altura3)
head(wrasse3)
tail(wrasse3)
summary(wrasse3)

#miños
ids<-unique(wrasse4$Idflota) 
long<-length(ids)
lista<-as.list(1:long)
names(lista)<-ids

for (i in 1:long) {lista[[i]]<-wrasse4[wrasse4$Idflota==ids[i],]}
for (i in 1:long) {lista[[i]]$longitud2<-ifelse(sum(lista[[i]]$longitud,na.rm=T) > 0,mean(lista[[i]]$longitud,na.rm=T),NA)}
for (i in 1:long) {lista[[i]]$altura2<-ifelse(sum(lista[[i]]$altura,na.rm=T) > 0,mean(lista[[i]]$altura,na.rm=T),NA)}

wrasse4<-do.call(rbind,lista)
head(wrasse4)
summary(wrasse4)

wrasse4$longitud3<-ifelse(is.na(wrasse4$longitud),wrasse4$longitud2,wrasse4$longitud)
wrasse4$altura3<-ifelse(is.na(wrasse4$altura),wrasse4$altura2,wrasse4$altura)
head(wrasse4)
tail(wrasse4)
summary(wrasse4)

wrasse4$longitud4<-ifelse(is.na(wrasse4$longitud3),
                           mean(wrasse4$longitud3,na.rm=T),wrasse4$longitud3) # Fill in with fleet average data for those vessels without any information
wrasse4$altura4<-ifelse(is.na(wrasse4$altura3),
                         mean(wrasse4$altura3,na.rm=T),wrasse4$altura3)
summary(wrasse4)

dim(wrasse2)
dim(wrasse3)
dim(wrasse4)
wrasse5<-data.frame(rbind(wrasse2,wrasse3,wrasse4))
dim(wrasse5)
summary(wrasse5)

wrasse5$Area<-wrasse5$Pieces*wrasse5$longitud4*wrasse5$altura4 # TRASMALLOS, MIÑ‘OS & VETAS Area (m2) 
summary(wrasse5$Area)
hist(wrasse5$Area)

#4.2.3# Differences between Nasas, Vetas and Miños

head(wrasse5)
wrasse5$Gear<-factor(wrasse5$Gear)
wrasse5$Gear<-relevel(wrasse5$Gear, ref = "TRASMALLOS")
levels(wrasse5$Gear)

table(wrasse5$Year,wrasse5$Month,wrasse5$Gear)
table(wrasse5$Year,wrasse5$ZoneO,wrasse5$Gear)

gear1<-lm(log(Area)~Gear,data=wrasse5)
summary(gear1)
plot(allEffects(gear1))

gear2<-lm(log(Depth)~Gear,data=wrasse5)
summary(gear2)
plot(allEffects(gear2))

gear3<-lm(log(Soak)~Gear,data=wrasse5)
summary(gear3)
plot(allEffects(gear3))

gear4<-lm(Month~Gear,data=wrasse5)
summary(gear4)
plot(allEffects(gear4))

gear5<-lm(GRT~Gear,data=wrasse5)
summary(gear5)
plot(allEffects(gear5))

gear6<-lm(Crew~Gear,data=wrasse5)
summary(gear6)
plot(allEffects(gear6))

gear7<-lm(ZoneO~Gear,data=wrasse5)
summary(gear7)
plot(allEffects(gear7))

# ----------------------------- #
#5# Add environmental data (At this time I run models without AR structure)

#5.1# Load upwelling and remove seasonal cycle

qx<-read.csv2(file="C:\\Users\\alex\\Dropbox\\casgass\\Upwelling.csv",header=T,dec=",",sep=",")

head(qx)
qx$QX<-qx$QX/1000 # Change scale

summary(qx$QX)
#Based on Álvarez-Salgado et al. (2002) extreme values where considered as follows in the unpwelling time series
IntQua<-3*(summary(qx$QX)[[5]]-summary(qx$QX)[[2]])

hlow<-summary(qx$QX)[[2]]-IntQua
hupp<-summary(qx$QX)[[5]]+IntQua

qx$QX<-ifelse(qx$QX < hlow,hlow,qx$QX)
qx$QX<-ifelse(qx$QX > hupp,hupp,qx$QX)

g1<-gam(QX~s(DoY,k=6,bs="cc")+s(Timer,k=4,bs="cr"),data=qx,na.action=na.exclude)
summary(g1)
#plot(g1,pages=1)

qx$qxAno<-residuals(g1)

head(qx)
qx$QY<-qx$QY/1000 # Change scale
g2<-gam(QY~s(DoY,k=6,bs="cc")+s(Timer,k=4,bs="cr"),data=qx,na.action=na.exclude)
summary(g2)
#plot(g2,pages=1)
qx$qyAno<-residuals(g2)

#5.2# Load sst and remove seasonal cycle for each zone

sst<-read.table(file="C:\\Users\\alex\\Dropbox\\casgass\\sstGal.selection.txt",header=T,dec=".",sep=" ")

head(sst)

z1<-205 # Zone 1
z2<-229 # Zone 2
z3<-228 # Zone 3
z4<-252 # Zone 4
z5<-276 # Zone 5
z6<-324 # Zone 6
z7<-327 # Zone 7
z8<-352 # Zone 8
z9<-379 # Zone 9

#5.2.1# SST Zone 1

g3<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z1,],na.action=na.exclude)
summary(g3)
#plot(g3,pages=1)
sstAnoZ1<-residuals(g3)

#5.2.2# SST zone 2

g4<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z2,],na.action=na.exclude)
summary(g4)
#plot(g4,pages=1)
sstAnoZ2<-residuals(g4)

#5.2.3# SST zone 3

g5<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z3,],na.action=na.exclude)
summary(g5)
#plot(g5,pages=1)
sstAnoZ3<-residuals(g5)

#5.2.4# SST zone 4

g6<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z4,],na.action=na.exclude)
summary(g6)
#plot(g6,pages=1)
sstAnoZ4<-residuals(g6)

#5.2.5# SST zone 5

g7<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z5,],na.action=na.exclude)
summary(g7)
#plot(g7,pages=1)
sstAnoZ5<-residuals(g7)

#5.2.6# SST zone 6

g8<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z6,],na.action=na.exclude)
summary(g8)
#plot(g8,pages=1)
sstAnoZ6<-residuals(g8)

#5.2.7# SST zone 7

g9<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z7,],na.action=na.exclude)
summary(g9)
#plot(g9,pages=1)
sstAnoZ7<-residuals(g9)

#5.2.8# SST zone 8

g10<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z8,],na.action=na.exclude)
summary(g10)
#plot(g10,pages=1)
sstAnoZ8<-residuals(g10)

#5.2.9# SST zone 9

g11<-gam(sst~s(doy,k=6,bs="cc")+s(timer,k=4,bs="cr"),data=sst[sst$box==z9,],na.action=na.exclude)
summary(g11)
#plot(g11,pages=1)
sstAnoZ9<-residuals(g11)

#5.3# Merge data sets

oceano<-cbind(qx,sstAnoZ1,sstAnoZ2,sstAnoZ3,sstAnoZ4,sstAnoZ5,
              sstAnoZ6,sstAnoZ7,sstAnoZ8,sstAnoZ9)

head(oceano)

#5.4# Calculate mean qx and sst for each pollack haul

#5.4.1# Set up the 'windows' defined as an 'habitat preference indicator'

win.Up<-18 # 18 days for upwelling (the three previous weeks)
win.Te<-6 # Days for sst: 3 (the previous 3 days), 6 (the previous week), 12 (the two previous weeks)

#5.4.2# Calculate environmental covariate (from the window day to the preceeding day of the catch)

wrasse5$QxM<-as.vector(rep(NA,dim(wrasse5)[1]))

for (i in 1:dim(wrasse5)[1]) {wrasse5$QxM[i]<-mean(oceano$qxAno[(which(oceano$Year==wrasse5$Year[i] &
                                                                           oceano$DoY==wrasse5$Julian[i])-win.Up):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                            oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)}

wrasse5$QyM<-as.vector(rep(NA,dim(wrasse5)[1]))

for (i in 1:dim(wrasse5)[1]) {wrasse5$QyM[i]<-mean(oceano$qyAno[(which(oceano$Year==wrasse5$Year[i] &
                                                                           oceano$DoY==wrasse5$Julian[i])-win.Up):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                            oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)}

wrasse5$sstM<-as.vector(rep(NA,dim(wrasse5)[1]))

for (i in 1:dim(wrasse5)[1]) {wrasse5$sstM[i]<-
                                 
                                 if (wrasse5$ZoneA[i]=="1") {mean(oceano$sstAnoZ1[(which(oceano$Year==wrasse5$Year[i] &
                                                                                            oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                             oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)} 
                               
                               else {if (wrasse5$ZoneA[i]=="2") {mean(oceano$sstAnoZ2[(which(oceano$Year==wrasse5$Year[i] &
                                                                                                oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                                 oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)} 
                                     
                                     else {if (wrasse5$ZoneA[i]=="3") {mean(oceano$sstAnoZ3[(which(oceano$Year==wrasse5$Year[i] &
                                                                                                      oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                                       oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)} 
                                           
                                           else {if (wrasse5$ZoneA[i]=="4") {mean(oceano$sstAnoZ4[(which(oceano$Year==wrasse5$Year[i] &
                                                                                                            oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                                             oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)} 
                                                 
                                                 else {if (wrasse5$ZoneA[i]=="5") {mean(oceano$sstAnoZ5[(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                  oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                                                   oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)} 
                                                       
                                                       else {if (wrasse5$ZoneA[i]=="6") {mean(oceano$sstAnoZ6[(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                        oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                                                         oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)} 
                                                             
                                                             else {if (wrasse5$ZoneA[i]=="7") {mean(oceano$sstAnoZ7[(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                              oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                                                               oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)} 
                                                                   
                                                                   else {if (wrasse5$ZoneA[i]=="8") {mean(oceano$sstAnoZ8[(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                    oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                                                                     oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)} 
                                                                         
                                                                         else {mean(oceano$sstAnoZ9[(which(oceano$Year==wrasse5$Year[i] &
                                                                                                             oceano$DoY==wrasse5$Julian[i])-win.Te):(which(oceano$Year==wrasse5$Year[i] &
                                                                                                                                                              oceano$DoY==wrasse5$Julian[i])-1)],na.rm=T)}
                                                                         
                                                                   }
                                                             }
                                                       }
                                                 }
                                           }
                                     }
                               }
}


# ----------------------------- #
#6# Preliminary steps

head(wrasse5)
summary(wrasse5)

names(wrasse5)
wrasse6<-wrasse5[,-c(30:37)]
dim(wrasse6)
names(wrasse6)

#6.1# Response variable

#6.1.1# Percentage of zeroes and other stuff

(length(which(wrasse6$Ntot==0))*100)/nrow(wrasse6) # Percentage of zeros
dim(wrasse6)[1] # Number of hauls
length(unique(wrasse6$Idflota)) # Number of vessels

wrasse6$Gear<-factor(wrasse6$Gear)
poll.dist<-as.data.frame(table(wrasse6$Ntot,wrasse6$Gear))
colnames(poll.dist)<-c("Count","Gear","Freq")
head(poll.dist)

pdf(file="respWrasse.pdf",width=10,height=10)

ggplot(data=poll.dist,aes(x=Count,y=Freq))+facet_wrap(~Gear,nrow=3,scales="free_y")+
  geom_bar(stat="identity")+
  scale_y_continuous("Frequency")+
  scale_x_discrete("Number of pouting caught per haul")+
  theme(axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=8))

dev.off()

par(mar=c(5,5,3,3))
hist(wrasse6$Ntot,prob=T,ylab="Probability",xlab="Nº of fished pouting",main="")

hist(log(wrasse6$Wtot))

plot(log(Ntot)~log(Wtot),data=wrasse6,ylab="Total number",xlab="Total biomass")

#6.1.2# Plots of spatial distribution of effort (m2 of net and fishing hour) 

p1<-ggplot(data=wrasse6,aes(x=Lon,y=Lat))+facet_wrap(~Gear,nrow=3,scales="free_y")

pdf(file="effortWrasse.pdf",width=8,height=12)

p1+geom_point(aes(color=Area))+
  scale_color_gradient(low="blue",high="red")

dev.off()

#6.1.3# Distribution of cases across explanatory variables

table(wrasse6$Year,wrasse6$Gear) # Few data in 1999
table(wrasse6$Year,wrasse6$Month,wrasse6$Gear)
table(wrasse6$fcrew,wrasse6$Gear)
table(wrasse6$fgear)
table(wrasse6$ZoneA,wrasse6$Gear)
table(wrasse6$ZoneO,wrasse6$Gear)

#problema de muestreo con NASAS, años submuestreados
#imposible obtener una serie temporal completa
#de todos modos haremos los modelos para nasas y enmalle

#6.2# Set up the offset

#6.2.1# Exploration of effort variables

head(wrasse6)

par(mfrow=c(2,3))
plot(GRT~Crew,data=wrasse6)
abline(lm(GRT~Crew,data=wrasse6),col="red",lwd=2)
plot(GRT~Area,data=wrasse6)
abline(lm(GRT~Area,data=wrasse6),col="red",lwd=2)
plot(GRT~Soak,data=wrasse6)
abline(lm(GRT~Soak,data=wrasse6),col="red",lwd=2)
plot(Crew~Area,data=wrasse6)
abline(lm(Crew~Area,data=wrasse6),col="red",lwd=2)
plot(Crew~Soak,data=wrasse6)
abline(lm(Crew~Soak,data=wrasse6),col="red",lwd=2)
plot(Area~Soak,data=wrasse6)
abline(lm(Area~Soak,data=wrasse6),col="red",lwd=2)

par(mfrow=c(1,2))
boxplot(Area~fcrew,data=wrasse6,notch=T,xlab="Crew",ylab="Area")
abline(lm(Area~fcrew,data=wrasse6),col="red",lwd=2)
boxplot(Soak~fcrew,data=wrasse6,notch=T,xlab="Crew",ylab="Soak time")
abline(lm(Soak~fcrew,data=wrasse6),col="red",lwd=2)

#6.3# Calculation of nighttime fishing

library(lubridate) # Necesaria para obtener horas y minutos de objetos POSIXct

clock.M <- function (t) {hour(t)*60 + minute(t)} # Calcula el minuto entre 1 y 1440 de un dia

pho<-read.csv2(file="C:\\Users\\alex\\Dropbox\\casgass\\Sunrise_set.csv",header=T,dec=".",sep=";")

head(pho)
pho<-pho[,c(1:3,7:12)] # Nos quedamos con las columnas que interesan

#6.3.1# Convertir a POSIXct y obtener informacion sobre las fechas y horas

head(wrasse6)

wrasse6$Deployment<-as.POSIXct(wrasse6$Deployment,format="%Y-%m-%d %H:%M:%S") 
wrasse6$Retrieval<-as.POSIXct(wrasse6$Retrieval,format="%Y-%m-%d %H:%M:%S")
pho$Rise<-as.POSIXct(pho$Rise,format="%d/%m/%Y %H:%M:%S")
pho$Set<-as.POSIXct(pho$Set,format="%d/%m/%Y %H:%M:%S")

wrasse6$yLarg<-year(wrasse6$Deployment) # Obtener año de largada
wrasse6$mLarg<-month(wrasse6$Deployment) # Obtener mes de largada
wrasse6$dLarg<-day(wrasse6$Deployment) # Obtener dia de largada
wrasse6$yVir<-year(wrasse6$Retrieval) # Obtener año de virada
wrasse6$mVir<-month(wrasse6$Retrieval) # Obtener mes de virada
wrasse6$dVir<-day(wrasse6$Retrieval) # Obtener dia de virada

wrasse6$hLarg<-clock.M(wrasse6$Deployment) # Calcula minuto largada
wrasse6$hVir<-clock.M(wrasse6$Retrieval) # Calcula minuto virada
pho$hRise<-clock.M(pho$Rise) # Calcula minuto Sunrise
pho$hSet<-clock.M(pho$Set) # Calcula minuto Sunset

#6.3.2# Obtener informacion cruzando los datos de pesca con fotoperiodo

# Obtener dia unico para la largada de entre todos los dias posibles desde
# el 1/1/1999 (dia 1) hasta el 31/12/2013 (dia 5479) 

wrasse6$tLarg<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$tLarg[i]<-pho$Timer[which(wrasse6$yLarg[i]==pho$Year & wrasse6$mLarg[i]==pho$Month & wrasse6$dLarg[i]==pho$Day)]} 

# Obtener dia unico para la virada de entre todos los dias posibles desde
# el 1/1/1999 (dia 1) hasta el 31/12/2013 (dia 5479) 

wrasse6$tVir<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$tVir[i]<-pho$Timer[which(wrasse6$yVir[i]==pho$Year & wrasse6$mVir[i]==pho$Month & wrasse6$dVir[i]==pho$Day)]}  

# Obtener minuto Sunrise el dia de largada

wrasse6$hRiseL<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$hRiseL[i]<-pho$hRise[which(wrasse6$tLarg[i]==pho$Timer)]} 

# Obtener minuto Sunset el dia de largada

wrasse6$hSetL<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$hSetL[i]<-pho$hSet[which(wrasse6$tLarg[i]==pho$Timer)]} 

# Obtener minuto Sunrise el dia de virada

wrasse6$hRiseV<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$hRiseV[i]<-pho$hRise[which(wrasse6$tVir[i]==pho$Timer)]} 

# Obtener minuto Sunset el dia de virada

wrasse6$hSetV<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$hSetV[i]<-pho$hSet[which(wrasse6$tVir[i]==pho$Timer)]} 

# Obtener minutos nocturnos transcurridos entre dias en los casos en que
# pasÃ³ mas de un dia entre largada y virada

wrasse6$minNigTot<-as.vector(rep(NA,dim(wrasse6)[1]))

for (i in 1:dim(wrasse6)[1]) {wrasse6$minNigTot[i]<-if (wrasse6$tVir[i]-wrasse6$tLarg[i]<=1) {0}
                               else {sum(pho$Night[which(pho$Timer==wrasse6$tLarg[i]+1) : which(pho$Timer== wrasse6$tVir[i])])}
} 

#6.3.3# Calcular minutos nocturnos

#6.3.3.1# Minutos nocturnos si la largada y la virada tienen lugar el mismo dia (6 combinaciones posibles) 

wrasse6$minNig1<-as.vector(rep(NA,dim(wrasse6)[1]))

for (i in 1:dim(wrasse6)[1]) {wrasse6$minNig1[i]<-if (wrasse6$tVir[i]-wrasse6$tLarg[i]==0) {
  
  # Largada y Virada ocurren antes del Sunrise
  ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] <= wrasse6$hRiseL[i],wrasse6$hVir[i]-wrasse6$hLarg[i],
         
         # Largada ocurre antes del Sunrise y Virada entre Sunrise y Sunset 
         ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] > wrasse6$hRiseL[i] & wrasse6$hVir[i] <= wrasse6$hSetL[i],wrasse6$hRiseL[i]-wrasse6$hLarg[i], 
                
                # Largada ocurre antes del Sunrise y Virada despues del Sunset
                ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] > wrasse6$hSetL[i],(wrasse6$hRiseL[i]-wrasse6$hLarg[i])+(wrasse6$hVir[i]-wrasse6$hSetL[i]),
                       
                       # Largada ocurre entre Sunrise y Sunset y Virada ocurre entre Sunrise y Sunset 
                       ifelse(wrasse6$hLarg[i] >= wrasse6$hRiseL[i] & wrasse6$hLarg[i] <= wrasse6$hSetL & wrasse6$hVir[i] >= wrasse6$hRiseL & wrasse6$hVir[i] <= wrasse6$hSetL[i],0,
                              
                              # Largada ocurre entre Sunrise y Sunset y Virada ocurre despues del Sunset
                              ifelse(wrasse6$hLarg[i] >= wrasse6$hRiseL[i] & wrasse6$hLarg[i] <= wrasse6$hSetL & wrasse6$hVir[i] > wrasse6$hSetL[i],wrasse6$hVir[i]-wrasse6$hSetL[i],
                                     
                                     # Largada y Virada ocurren despues del Sunset
                                     wrasse6$hVir[i]-wrasse6$hLarg[i])))))
  
} else {0}

} 

#6.3.3.2# Minutos nocturnos si la virada tiene lugar al dia siguiente de la largada (9 combinaciones posibles)

wrasse6$minNig2<-as.vector(rep(NA,dim(wrasse6)[1]))

for (i in 1:dim(wrasse6)[1]) {wrasse6$minNig2[i]<-if (wrasse6$tVir[i]-wrasse6$tLarg[i]==1) {
  
  # Largada ocurre antes del Sunrise y Virada ocurre antes del Sunrise del dia siguiente 
  ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] <= wrasse6$hRiseV[i],(wrasse6$hRiseL[i]-wrasse6$hLarg[i])+((1440-wrasse6$hSetL[i])+wrasse6$hVir[i]),
         
         # Largada ocurre antes del Sunrise y Virada ocurre entre el Sunrise y el Sunset del dia siguiente
         ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] > wrasse6$hRiseV[i] & wrasse6$hVir[i] <= wrasse6$hSetV[i],(wrasse6$hRiseL[i]-wrasse6$hLarg[i])+((1440-wrasse6$hSetL[i])+wrasse6$hRiseV[i]),
                
                # Largada ocurre antes del Sunrise y Virada ocurre despues del Sunset del dia siguiente
                ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] > wrasse6$hSetV[i],(wrasse6$hRiseL[i]-wrasse6$hLarg[i])+(1440-wrasse6$hSetL[i])+wrasse6$hRiseV[i]+(wrasse6$hVir[i]-wrasse6$hSetV[i]),
                       
                       # Largada ocurre entre el Sunrise y el Sunset y Virada ocurre antes del Sunrise del dia siguiente 
                       ifelse(wrasse6$hLarg[i] >= wrasse6$hRiseL[i] & wrasse6$hLarg[i] <= wrasse6$hSetL[i] & wrasse6$hVir[i] <= wrasse6$hRiseV[i],(1440-wrasse6$hSetL[i])+wrasse6$hVir[i],
                              
                              # Largada ocurre entre el Sunrise y el Sunset y Virada ocurre entre el Sunrise y el Sunset del dia siguiente 
                              ifelse(wrasse6$hLarg[i] >= wrasse6$hRiseL[i] & wrasse6$hLarg[i] <= wrasse6$hSetL[i] & wrasse6$hVir[i] >= wrasse6$hRiseV[i] & wrasse6$hVir[i] <= wrasse6$hSetV[i],(1440-wrasse6$hSetL[i])+wrasse6$hRiseV[i],
                                     
                                     # Largada ocurre entre el Sunrise y el Sunset y Virada ocurre despues del Sunset del dia siguiente
                                     ifelse(wrasse6$hLarg[i] >= wrasse6$hRiseL[i] & wrasse6$hLarg[i] <= wrasse6$hSetL[i] & wrasse6$hVir[i] > wrasse6$hSetV[i],(1440-wrasse6$hSetL[i])+wrasse6$hRiseV[i]+(wrasse6$hVir[i]-wrasse6$hSetV[i]),
                                            
                                            # Largada ocurre despues del Sunset y virada ocurre antes del Sunrise del dia siguiente
                                            ifelse(wrasse6$hLarg[i] > wrasse6$hSetL[i] & wrasse6$hVir[i] <= wrasse6$hRiseV[i],(1440-wrasse6$hLarg[i])+wrasse6$hVir[i],
                                                   
                                                   # Largada ocurre despues del Sunset y virada ocurre entre el Sunrise y el Sunset del dia siguiente
                                                   ifelse(wrasse6$hLarg[i] > wrasse6$hSetL[i] & wrasse6$hVir[i] >= wrasse6$hRiseV[i] & wrasse6$hVir[i] <= wrasse6$hSetV[i],(1440-wrasse6$hLarg[i])+wrasse6$hRiseV[i],
                                                          
                                                          # Largada ocurre despues del Sunset y virada ocurre despues del Sunset del dia siguiente
                                                          (1440-wrasse6$hLarg[i])+wrasse6$hRiseV[i]+(wrasse6$hVir[i]-wrasse6$hSetV[i])))))))))
  
} else {0}

}

#6.3.3.3# Minutos nocturnos si entre la largada y la virada pasa mas de un dia (las mismas 9 combinaciones posibles)

wrasse6$minNig3<-as.vector(rep(NA,dim(wrasse6)[1]))

for (i in 1:dim(wrasse6)[1]) {wrasse6$minNig3[i]<-if (wrasse6$tVir[i]-wrasse6$tLarg[i] > 1) {
  
  # Largada ocurre antes del Sunrise y Virada ocurre antes del Sunrise varios dias despues
  ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] <= wrasse6$hRiseV[i],(wrasse6$hRiseL[i]-wrasse6$hLarg[i])+(wrasse6$minNigTot[i]-(wrasse6$hRiseV[i]-wrasse6$hVir[i])),
         
         # Largada ocurre antes del Sunrise y Virada ocurre entre el Sunrise y el Sunset varios dias despues
         ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] > wrasse6$hRiseV[i] & wrasse6$hVir[i] <= wrasse6$hSetV[i],(wrasse6$hRiseL[i]-wrasse6$hLarg[i])+wrasse6$minNigTot[i],
                
                # Largada ocurre antes del Sunrise y Virada ocurre despues del Sunset varios dias despues
                ifelse(wrasse6$hLarg[i] < wrasse6$hRiseL[i] & wrasse6$hVir[i] > wrasse6$hSetV[i],(wrasse6$hRiseL[i]-wrasse6$hLarg[i])+wrasse6$minNigTot[i]+(wrasse6$hVir[i]-wrasse6$hSetV[i]),
                       
                       # Largada ocurre entre el Sunrise y el Sunset y Virada ocurre antes del Sunrise varios dias despues
                       ifelse(wrasse6$hLarg[i] >= wrasse6$hRiseL[i] & wrasse6$hLarg[i] <= wrasse6$hSetL[i] & wrasse6$hVir[i] <= wrasse6$hRiseV[i],wrasse6$minNigTot[i]-(wrasse6$hRiseV[i]-wrasse6$hVir[i]),
                              
                              # Largada ocurre entre el Sunrise y el Sunset y Virada ocurre entre el Sunrise y el Sunset varios dias despues
                              ifelse(wrasse6$hLarg[i] >= wrasse6$hRiseL[i] & wrasse6$hLarg[i] <= wrasse6$hSetL[i] & wrasse6$hVir[i] >= wrasse6$hRiseV[i] & wrasse6$hVir[i] <= wrasse6$hSetV[i],wrasse6$minNigTot[i],
                                     
                                     # Largada ocurre entre el Sunrise y el Sunset y Virada ocurre despues del Sunset varios dias despues
                                     ifelse(wrasse6$hLarg[i] >= wrasse6$hRiseL[i] & wrasse6$hLarg[i] <= wrasse6$hSetL[i] & wrasse6$hVir[i] > wrasse6$hSetV[i],wrasse6$minNigTot[i]+(wrasse6$hVir[i]-wrasse6$hSetV[i]),
                                            
                                            # Largada ocurre despues del Sunset y virada ocurre antes del Sunrise varios dias despues
                                            ifelse(wrasse6$hLarg[i] > wrasse6$hSetL[i] & wrasse6$hVir[i] <= wrasse6$hRiseV[i],wrasse6$minNigTot[i]-((wrasse6$hLarg[i]-wrasse6$hSetL[i])-(wrasse6$hRiseV[i]-wrasse6$hVir[i])),
                                                   
                                                   # Largada ocurre despues del Sunset y virada ocurre entre el Sunrise y el Sunset varios dias despues
                                                   ifelse(wrasse6$hLarg[i] > wrasse6$hSetL[i] & wrasse6$hVir[i] >= wrasse6$hRiseV[i] & wrasse6$hVir[i] <= wrasse6$hSetV[i],wrasse6$minNigTot[i]-(wrasse6$hLarg[i]-wrasse6$hSetL[i]),
                                                          
                                                          # Largada ocurre despues del Sunset y virada ocurre despues del Sunset varios dias despues
                                                          (wrasse6$minNigTot[i]-(wrasse6$hLarg[i]-wrasse6$hSetL[i]))+(wrasse6$hVir[i]-wrasse6$hSetV[i])))))))))
  
} else {0}

}

#6.3.4# Nueva variable '% minutos nocturnos'

#6.3.4.1# Check that soak from UTPB is correct

# wrasse6$Soaktime<-as.numeric(difftime(wrasse6$Retrieval,wrasse6$Deployment,units="mins"))

# plot(wrasse6$Soak~Soaktime)
# soakTest<-ifelse(wrasse6$Soak==Soaktime,0,1)
# sum(soakTest) # Los dos valores son iguales asi que usamos el suyo
# nn<-which(soakTest==1,)
# wrasse6[nn,]

#6.3.4.2# Calculate proportion of night time

wrasse6$caladoNight<-as.vector(rep(NA,dim(wrasse6)[1])) 

for (i in 1:dim(wrasse6)[1]) {wrasse6$caladoNight[i]<-if (wrasse6$tVir[i]==wrasse6$tLarg[i]) {
  
  # % minutos nocturnos si largada y virada ocurren el mismo dia
  round(wrasse6$minNig1[i]/wrasse6$Soak[i],2)
  
} else {if (wrasse6$tVir[i]-wrasse6$tLarg[i]==1) {
  
  # % minutos nocturnos si virada ocurre al dia siguiente de largada
  round(wrasse6$minNig2[i]/wrasse6$Soak[i],2)
  
} else {
  
  # % minutos nocturnos si virada ocurre varios dias despues de largada
  round(wrasse6$minNig3[i]/wrasse6$Soak[i],2)
}
}
}

summary(wrasse6$caladoNight)
hist(wrasse6$caladoNight)

#6.3.5# Standarizar hora de largada segun fotoperiodo

#6.3.5.1# Obtener minuto Sunrise el dia siguiente de largada

names(wrasse6)

wrasse6$hRiseL1<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$hRiseL1[i]<-pho$hRise[which(wrasse6$tLarg[i]+1==pho$Timer)]} 

#6.3.5.2# Corregir hora de largada segun fotoperiodo: si la hora de largada ocurre antes del
# sunset restamos hora de largada-hora sunrise (con lo que largadas antes de sunrise
# tendran valores negativos, y largadas despues de sunrise hasta sunset tedran valores
# positivos. Si la largada ocurre despues del sunset los valores seran muy negativos
# por que los contaremos respecto del sunrise del dia siguiente)

wrasse6$hLargC<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$hLargC[i]<-if (wrasse6$hLarg[i] <= wrasse6$hSetL[i]) {
  
  wrasse6$hLarg[i] - wrasse6$hRiseL[i]
  
} else {
  
  ((1440 - wrasse6$hLarg[i]) + wrasse6$hRiseL1[i])*-1
  
}
}

#6.3.6# Standarizar hora de virada segun fotoperiodo

#6.3.6.1# Obtener minuto Sunrise el dia siguiente de largada

names(wrasse6)

wrasse6$hSetV1<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$hSetV1[i]<-pho$hSet[which(wrasse6$tVir[i]-1==pho$Timer)]} 

#6.3.6.2# Corregir hora de virada segun fotoperiodo: si la hora de virada ocurre antes del
# sunrise restamos hora de virada-hora sunset del dia anterior (con lo que viradas mas
# lejos del sunset tienen valores muy negativos). Si la virada ocurre despues del sunrise
# los valores seran positivos

wrasse6$hVirC<-as.vector(rep(NA,dim(wrasse6)[1]))
for (i in 1:dim(wrasse6)[1]) {wrasse6$hVirC[i]<-if (wrasse6$hVir[i] <= wrasse6$hRiseV[i]) {
  
  ((1440 - wrasse6$hSetV1[i]) + wrasse6$hVir[i])*-1
  
} else {
  
  wrasse6$hVir[i] - wrasse6$hRiseV[i]
  
}
}

#6.3.7# Uso de factores (vale 0 si la proporcion de calado nocturno es
# mayor que 75%, si esta entre el 25 y 75% vale 1, y si es menor o igual que
# el 25%, esto es, el lance es practicamente todo diurno, vale 2)

wrasse6$Period<-factor(ifelse(wrasse6$caladoNight <= 0.25,2,
                               ifelse(wrasse6$caladoNight > 0.25 & wrasse6$caladoNight <= 0.75,1,0)))

#6.4 added by Alex Alonso

#we are not going to work with it

#añadimos una nueva variable indicativa del tamaño promedio de los individuos capturados
#Wtot/Ntot puede ser utilizada como indicativo de cambios de distribución ontogénicos???
#summary(wrasse6$Ntot)
#summary(wrasse6$Wtot)
#plot(wrasse6$Ntot,wrasse6$Wtot)
#necesitamos eliminar NA de los datos Wtot para poder calcular esta nueva variable
#wrasse6$Wtot<-ifelse(wrasse6$Ntot==0,0,wrasse6$Wtot)
#par(mfrow=c(1,2))
#dotchart(wrasse6$Wtot)
#hist(wrasse6$Wtot,breaks=100)
#summary(wrasse6$Wtot)
#plot(wrasse6$Ntot,wrasse6$Wtot)
#calculamos la relación
#wrasse6$AvgSize<-ifelse(wrasse6$Ntot==0,0,wrasse6$Wtot/wrasse6$Ntot)
#hist(wrasse6$AvgSize,breaks=100)
#summary(wrasse6$AvgSize)
#does it make sense 0 values for this variable?
#plot(wrasse6$Ntot,wrasse6$AvgSize)


# ----------------------------- #
#7# Exploratory steps before modelling

#7.1# Distribution of all potential explanatory variables

#likely outlier in enmalle >500
which(wrasse6$Depth>500) # Solo un lance por encima de 200 m (utpb dice que es real; pero no representativo de la pesquería)
#wrasse7<-wrasse6[wrasse6$Depth<500,]
#no lo eliminamos
wrasse7<-wrasse6
levels(wrasse7$Gear)
levels(wrasse7$fgear)
wrasse7$Gear <- relevel(wrasse7$Gear, ref="TRASMALLOS") 
wrasse7$fgear <- relevel(wrasse7$fgear, ref="TRASMALLOS") 

d1<-ggplot(data=wrasse7,aes(x=Crew))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Crew")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))
d2<-ggplot(data=wrasse7,aes(x=GRT))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("GRT")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d3<-ggplot(data=wrasse7,aes(x=Area))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous(expression(paste("Gear area"," ","(",m^2,")")))+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d4<-ggplot(data=wrasse7,aes(x=Soak/1440))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Soak time (days)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d5<-ggplot(data=wrasse7,aes(x=hLarg/60))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Deployment (Local time)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d6<-ggplot(data=wrasse7,aes(x=hVir/60))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Retrieval (Local time)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d7<-ggplot(data=wrasse7,aes(x=Julian))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Day of the Year")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d8<-ggplot(data=wrasse7,aes(x=Lat))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Latitude (ºN)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d9<-ggplot(data=wrasse7,aes(x=Lon))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Longitude (ºW)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d10<-ggplot(data=wrasse7,aes(x=Depth))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Depth (m)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d11<-ggplot(data=wrasse7,aes(x=QxM))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous(expression(paste("Qx"," ","(",m^3,s^-1,km^-1,") Ã—",10^3)))+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d12<-ggplot(data=wrasse7,aes(x=sstM))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("SST (ºC)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d13<-ggplot(data=wrasse7,aes(x=AvgSize))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Individual size (kg)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  theme(legend.position="none")
d14<-ggplot(data=wrasse7,aes(x=caladoNight))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Night-time soak (%)")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position="none")

pdf(file="explWrasse.pdf",width=20,height=10)

multiplot(d1,d7,d8,d2,d5,d9,d3,d6,d10,d4,d14,d12,cols=4)

dev.off()

ggplot(data=wrasse7,aes(x=hLargC))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Corrected deployment time")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))

ggplot(data=wrasse7,aes(x=hVirC))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Corrected retrieval time")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Trasmallos","Miños","Vetas"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))

summary(wrasse7[,c("GRT","Area","hLarg","hVir","Depth","Soak")])

#7.2 Simple relationships

names(wrasse7)
write.table(wrasse7,"wrasse_dat.txt",col.names = TRUE,sep=",")

wrasse.dat<-read.table("wrasse_dat.txt",header=T,dec=".",sep=",")

names(wrasse.dat)
dim(wrasse.dat)
summary(wrasse.dat)


#data base
wrasse.dat<-wrasse.dat[,-c(34:39,42:51,53,55)]
names(wrasse.dat)


par(mfrow=c(4,4))
plot(log(Ntot+1)~log(GRT),data=wrasse.dat)
boxplot(log(Ntot+1)~fcrew,data=wrasse.dat)
boxplot(log(Ntot+1)~fyear,data=wrasse.dat)
plot(log(Ntot+1)~Julian,data=wrasse.dat)
plot(log(Ntot+1)~Depth,data=wrasse.dat)
plot(log(Ntot+1)~QxM,data=wrasse.dat)
plot(log(Ntot+1)~QyM,data=wrasse.dat)
plot(log(Ntot+1)~sstM,data=wrasse.dat)
boxplot(log(Ntot+1)~ZoneA,data=wrasse.dat)
boxplot(log(Ntot+1)~ZoneO,data=wrasse.dat)
plot(log(Ntot+1)~Soak,data=wrasse.dat)
plot(log(Ntot+1)~caladoNight,data=wrasse.dat)
hist(wrasse.dat$GRT)
hist(wrasse.dat$Depth)
hist(wrasse.dat$QxM)
hist(wrasse.dat$sstM)

#7.3# Adding Offset

wrasse.dat$OffSet<-wrasse.dat$Area*(wrasse.dat$Soak/60)
summary(wrasse.dat$OffSet);hist(log(wrasse.dat$OffSet),xlab="OffSet",main="")

#7.4# Collinearity

head(wrasse.dat)
dim(wrasse.dat)

#transform depth and GRT for enmalles
wrasse.dat$lGRT<-log(wrasse.dat$GRT)
wrasse.dat$lDepth<-log(wrasse.dat$Depth)

#VIF
vifs<-c("lGRT","Crew","Year","Julian","Depth","QxM","sstM","caladoNight")
corvif(wrasse.dat[,vifs])
#collianlity crew~GRT

# ----------------------------- #
#8# Modelling of catches using non-mixed models

#factores
names(wrasse.dat)
str(wrasse.dat)
wrasse.dat$fyear<-as.factor(wrasse.dat$fyear)
wrasse.dat$fcrew<-as.factor(wrasse.dat$fcrew)
wrasse.dat$fZoneA<-as.factor(wrasse.dat$fZoneA)
wrasse.dat$fZoneO<-as.factor(wrasse.dat$fZoneO)

summary(wrasse.dat$GRT);hist(wrasse.dat$GRT)
#so we are going to do 5 main classes to avoid extrange patterns in linear response
# "C5", "C5-10", "C10-15","C15-20","C>20"
wrasse.dat$fGRT<-as.factor(ifelse(wrasse.dat$GRT<5,"C5",
                                       ifelse(wrasse.dat$GRT<10,"C5-10",
                                              ifelse(wrasse.dat$GRT<15,"C10-15",
                                                     ifelse(wrasse.dat$GRT<20,"C15-20","C>20")))))
wrasse.dat$fGRT<-factor(wrasse.dat$fGRT,levels = c("C5", "C5-10", "C10-15","C15-20","C>20"))
str(wrasse.dat$fGRT)
boxplot(Ntot~fGRT,data=wrasse.dat)
table(wrasse.dat$fGRT)

wrasse.dat$Gear<-relevel(wrasse.dat$Gear,ref="TRASMALLOS")
levels(wrasse.dat$Gear)

#8.1# GAM Poisson modelling

poiss0<-gam(Ntot~offset(log(OffSet))+
              lGRT*Gear+
              fyear*fZoneO+s(Julian,k=6,bs="cc")+
              s(lDepth,by=Gear,k=3)+
              s(sstM,k=3)+
              s(caladoNight,by=Gear,k=3)+
              Seafloor,
            family=poisson,data=wrasse.dat)

summary(poiss0)
anova(poiss0)
plot(poiss0,pages=1,all.terms=T,scale=0,shade="T")
sum(residuals(poiss0,type="pearson")^2)/poiss0$df.resid # Overdispersion

#8.2# GAM NegBin modelling

negbin0<-gam(Ntot~offset(log(OffSet))+
                    lGRT*Gear+
                    fyear*fZoneO+s(Julian,k=6,bs="cc")+
                    s(lDepth,by=Gear,k=3)+
                    s(sstM,k=3)+
                    s(caladoNight,by=Gear,k=3)+
                    Seafloor,
                  family=nb(),data=wrasse.dat)

summary(negbin0)
anova(negbin0)
plot(negbin0,pages=1,all.terms=T,scale=0,shade="T")
sum(residuals(negbin0,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(negbin0))+1)) 

#8.3# GLM Poisson modelling

poiss1<-glm(Ntot~offset(log(OffSet))+
              Gear*lGRT+
              fyear*fZoneO+
              poly(Julian,2)+
              lDepth*Gear+
              sstM+
              caladoNight*Gear+
              Seafloor,
            family=poisson,data=wrasse.dat)

summary(poiss1)
Anova(poiss1)
plot(allEffects(poiss1))
sum(residuals(poiss1,type="pearson")^2)/poiss1$df.resid # Overdispersion

#8.4# Negative Binomial modelling

negbin1<-glm.nb(Ntot~offset(log(OffSet))+
                  Gear*lGRT+
                  fyear*fZoneO+
                  poly(Julian,2)+
                  lDepth*Gear+
                  sstM+
                  caladoNight*Gear+
                  Seafloor,
                data=wrasse.dat)

summary(negbin1)
Anova(negbin1)
plot(allEffects(negbin1))
sum(residuals(negbin1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(negbin1))+1))

par(mfrow=c(1,3))
hist(residuals(negbin1,type="deviance"))
plot(residuals(negbin1,type="deviance")~log(fitted(negbin1)))
boxplot(residuals(negbin1,type="deviance")~wrasse.dat$Idflota)
abline(0,0,col="red")
#even negative binomial does not fix overdispersion

#we are not going to try to modeled with Hilbe package (msme, COUNT)
#neither NB mixed
#and go directly to ZI models

odTest(negbin1)

#8.5# ZERO INFLATED

#8.5.1# exploring binomial 
wrasse.dat$binNtot<-ifelse(wrasse.dat$Ntot==0,0,1)

bin0<-gam(binNtot~lDepth+sstM+caladoNight+Seafloor,
          family=binomial,data=wrasse.dat)
summary(bin0)
anova(bin0)
plot(bin0,pages=1,all.terms=T,scale=0,shade="T")
#depth and calado night result in positive effects on pouting catch

#8.5.2# ZI Models
zipoiss1<-zeroinfl(Ntot~offset(log(OffSet))+
                     Gear*lGRT+
                     fyear*fZoneO+
                     poly(Julian,2)+
                     lDepth*Gear+
                     sstM+
                     caladoNight*Gear+
                     Seafloor|1,
                   dist="poisson",link="logit",data=wrasse.dat)

summary(zipoiss1)
sum(residuals(zipoiss1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zipoiss1))))

zipoiss2<-zeroinfl(Ntot~offset(log(OffSet))+
                     Gear*lGRT+
                     fyear*fZoneO+
                     poly(Julian,2)+
                     lDepth*Gear+
                     sstM+
                     caladoNight*Gear+
                     Seafloor|lDepth+Seafloor,
                   dist="poisson",link="logit",data=wrasse.dat)

summary(zipoiss2)
sum(residuals(zipoiss2,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zipoiss2))))

zinb1<-zeroinfl(Ntot~offset(log(OffSet))+
                  Gear*lGRT+
                  fyear*fZoneO+
                  poly(Julian,2)+
                  lDepth*Gear+
                  sstM+
                  caladoNight*Gear+
                  Seafloor|1,
                dist="negbin",link="logit",data=wrasse.dat)

summary(zinb1)
sum(residuals(zinb1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zinb1))+1))

zinb2<-zeroinfl(Ntot~offset(log(OffSet))+
                  Gear*lGRT+
                  fyear*fZoneO+
                  poly(Julian,2)+
                  lDepth*Gear+
                  sstM+
                  caladoNight*Gear+
                  Seafloor|lDepth+Seafloor,
                dist="negbin",link="logit",data=wrasse.dat)

summary(zinb2)
sum(residuals(zinb2,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zinb2))+1))

#8.5.3# Mixed model NB
#http://glmmadmb.r-forge.r-project.org/glmmADMB.html

#glmmadmb does not allow to include explanatory variables in ZI part
wrasse.dat$Idflota<-factor(wrasse.dat$Idflota)

zinbm1 <- glmmadmb(Ntot~Gear*lGRT+
                     fyear*fZoneO+
                     poly(Julian,2)+
                     lDepth*Gear+
                     sstM+
                     caladoNight*Gear+
                     Seafloor+
                     offset(log(OffSet))+
                     (1|Idflota),
                   data=wrasse.dat,
                   zeroInflation=TRUE,
                   family="nbinom")

sum(residuals(zinbm1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zinbm1))+1))
#it does not work :(
#Parameters were estimated, but not standard errors were not

#8.6# full model comparison
# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

#8.6.1# AIC/BIC

#AIC
AICtab(poiss0,poiss1,negbin0,negbin1,zipoiss1,zipoiss2,zinb1,zinb2)

#BIC
aic=AIC(poiss0,poiss1,negbin0,negbin1,zipoiss1,zipoiss2,zinb1,zinb2,k=2)
aic[order(aic$AIC), ]
bic=AIC(poiss0,poiss1,negbin0,negbin1,zipoiss1,zipoiss2,zinb1,zinb2,k=log(dim(wrasse.dat)[1]))
bic[order(bic$AIC), ]

#8.6.2# overdispersion parameter

GAMpoisson=sum(residuals(poiss0,type="pearson")^2)/poiss0$df.resid
GLMpoisson=sum(residuals(poiss1,type="pearson")^2)/poiss1$df.resid
GLMNB=sum(residuals(negbin1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(negbin1))+1)) 
GAMNB=sum(residuals(negbin0,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(negbin0))+1)) 
ZIpoisson1=sum(residuals(zipoiss1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zipoiss1))))
ZIpoisson2=sum(residuals(zipoiss2,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zipoiss2))))
ZINB1=sum(residuals(zinb1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zinb1))+1))
ZINB2=sum(residuals(zinb2,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(zinb2))+1))

model=c("GAMpoisson","GLMpoisson","GLMNB","GAMNB","ZIpoisson1","ZIpoisson2","ZINB1","ZINB2")
theta=c(GAMpoisson,GLMpoisson,GLMNB,GAMNB,ZIpoisson1,ZIpoisson2,ZINB1,ZINB2)
M=data.frame(cbind(model,theta))
M=M[order(M$theta),]
M

#If the dispersion parameter is significantly greater than one, indicating overdispersion
#(variance greater than the mean),
#then the scale parameter should be used to adjust the variance.  
#Failing to account for the overdispersion can result in inflated test statistics. 
#However, when the dispersion parameter is less than one, then the test statistics become more conservative,
#which is not considered as much of a problem.

#8.6.3# loglikelihood

fm<-list("GAMPois"=poiss0,"GLMPois"=poiss1,
         "GAMNB"=negbin0,"GLMNB"=negbin1,
         "ZIPOISS1"=zipoiss1,"ZIPOISS2"=zipoiss2,
         "ZINB1"=zinb1,"ZINB2"=zinb2)

logliks<-rbind(logLik=sapply(fm,function (x) round(logLik(x),digits=0)),
                        Df=sapply(fm,function (x) attr(logLik(x),"df")))
logliks

#8.6.4# Vuong test to compare NB vs ZINB
# do not run with random models
#On the Misuse of The Vuong Test 
#http://cybermetrics.wlv.ac.uk/paperdata/misusevuong.pdf
#It is beyond doubt that the widespread practice of using Vuong's test
#for non-nested models as a test of zero-inflation is erroneous

vuong(poiss1,negbin1)
vuong(negbin1,zinb1)
vuong(zinb1,zinb2)

#8.6.5# Predicting probabilities

phat.pois<-predprob(poiss1)
phat.pois.mn<-apply(phat.pois,2,mean)
phat.nb<-predprob(negbin1)
phat.nb.mn<-apply(phat.nb,2,mean)
phat.zipoiss1<-predprob(zipoiss1)
phat.zipoiss.mn1<-apply(phat.zipoiss1,2,mean)
phat.zipoiss2<-predprob(zipoiss2)
phat.zipoiss.mn2<-apply(phat.zipoiss2,2,mean)
phat.zinb1<-predprob(zinb1)
phat.zinb.mn1<-apply(phat.zinb1,2,mean)
phat.zinb2<-predprob(zinb2)
phat.zinb.mn2<-apply(phat.zinb2,2,mean)

with(wrasse.dat,{
  hist(Ntot,prob=TRUE,col="grey",breaks=seq(-0.5,1450,1),xlab="",
       main="",ylim=c(0,1),xlim=c(-0.5,30))
  lines(x=seq(0,344,1),y=phat.pois.mn,type="l",lwd=2,col=1) #negro
  lines(x=seq(0,344,1),y=phat.nb.mn,type="l",lwd=2,col=2) #rojo
  lines(x=seq(0,344,1),y=phat.zipoiss.mn1,type="l",lwd=2,col=3) #verde
  lines(x=seq(0,344,1),y=phat.zipoiss.mn2,type="l",lwd=2,col=4) #azul
  lines(x=seq(0,344,1),y=phat.zinb.mn1,type="l",lwd=2,col=5) #turquesa
  lines(x=seq(0,344,1),y=phat.zinb.mn2,type="l",lwd=2,col=6) #rosa
})

#8.6.6# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

aics<-AIC(poiss0,poiss1,negbin0,negbin1,zipoiss1,zipoiss2,zinb1,zinb2,k=2) # AIC
bics<-AIC(poiss0,poiss1,negbin0,negbin1,zipoiss1,zipoiss2,zinb1,zinb2,k=log(dim(wrasse.dat)[1])) # BIC

logliks

zeroes<-round(c("Obs"=sum(wrasse.dat$Ntot<1),
                "GLMPois"=sum(dpois(0,fitted(poiss1))),
                "GLMNB"=sum(dnbinom(0,mu=fitted(negbin1),size=negbin1$theta)),
                "ZIPOISS1"=sum(predict(zipoiss1,type="prob")[,1]),
                "ZIPOISS2"=sum(predict(zipoiss2,type="prob")[,1]),
                "ZINB1"=sum(predict(zinb1,type="prob")[,1]),
                "ZINB2"=sum(predict(zinb2,type="prob")[,1])))

modelsComp<-data.frame(cbind(logliks[c(4,8,10,12,14,16)],bics[c(2,4:8),2],aics[c(2,4:8),2],logliks[c(3,7,9,11,13,15)],zeroes[2:7]))
colnames(modelsComp)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes[[1]],")"))
modelsComp$deltaAIC<-round(modelsComp$AIC-min(modelsComp$AIC),2)
modelsComp$deltaBIC<-round(modelsComp$BIC-min(modelsComp$BIC),2)
modelsComp$ModLikelihood<-round(exp(-modelsComp$deltaAIC/2),2)
modelsComp$AICweight<-round(modelsComp$ModLikelihood/sum(modelsComp$ModLikelihood),2)
modelsComp<-modelsComp[order(modelsComp$deltaAIC),]
modelsComp

#final candidate models

#NO random effects
summary(zinb2)

#Question about the over-dispersion parameter is in general a tricky one.
#A large over-dispersion parameter could be due to a miss-specified model or
#could be due to a real process with over-dispersion.
#Adding an over-dispersion problem does not necessarily improve a miss-specified model.
#http://www.ats.ucla.edu/stat/r/dae/zinbreg.htm

#8.7# Selection of predictors

#8.7.1# nasa model selection
ZINB1<-zinb2
summary(ZINB1)

#year:zone
ZINB2<-zeroinfl(Ntot~offset(log(OffSet))+
                  Gear*lGRT+
                  fyear+fZoneO+
                  poly(Julian,2)+
                  lDepth*Gear+
                  sstM+
                  caladoNight*Gear+
                  Seafloor|lDepth+Seafloor,
                dist="negbin",link="logit",data=wrasse.dat)
summary(ZINB2)
lrtest(ZINB1,ZINB2)
vuong(ZINB1,ZINB2)
AICtab(ZINB1,ZINB2)

#Gear:caladoNight (reintroducing interaction)
ZINB3<-zeroinfl(Ntot~offset(log(OffSet))+
                  Gear*lGRT+
                  fyear*fZoneO+
                  poly(Julian,2)+
                  lDepth*Gear+
                  sstM+
                  caladoNight+
                  Seafloor|lDepth+Seafloor,
                dist="negbin",link="logit",data=wrasse.dat)
summary(ZINB3)
lrtest(ZINB2,ZINB3)
lrtest(ZINB1,ZINB3)
vuong(ZINB2,ZINB3)
vuong(ZINB1,ZINB3)
AICtab(ZINB2,ZINB3)
AICtab(ZINB1,ZINB3)

#sstM (reintroducing interaction)
ZINB4<-zeroinfl(Ntot~offset(log(OffSet))+
                  Gear*lGRT+
                  fyear*fZoneO+
                  poly(Julian,2)+
                  lDepth*Gear+
                  caladoNight*Gear+
                  Seafloor|lDepth+Seafloor,
                dist="negbin",link="logit",data=wrasse.dat)
summary(ZINB4)
lrtest(ZINB3,ZINB4)
lrtest(ZINB1,ZINB4)
vuong(ZINB3,ZINB4)
vuong(ZINB1,ZINB4)
AICtab(ZINB3,ZINB4)
AICtab(ZINB1,ZINB4)

sum(residuals(ZINB1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(ZINB1))+1))
#slightly overdispersed, so, this can lead to inflated significance

#8.7.2# comparing models
AICzinb<-AIC(ZINB1,ZINB2,ZINB3,ZINB4,k=2) # AIC
BICzinb<-AIC(ZINB1,ZINB2,ZINB3,ZINB4,k=log(dim(wrasse.dat)[1])) # BIC

fm.zinb<-list("M1"=ZINB1,"M2"=ZINB2,"M3"=ZINB3,"M4"=ZINB4)
Logliks.zinb<-rbind(logLik=sapply(fm.zinb,function (x) round(logLik(x),digits=0)),
                   Df=sapply(fm.zinb,function (x) attr(logLik(x),"df")))

ZERO.zinb<-round(c("Obs"=sum(wrasse.dat$Ntot<1),
                  "M1"=sum(predict(ZINB1,type="prob")[,1]),
                  "M2"=sum(predict(ZINB2,type="prob")[,1]),
                  "M3"=sum(predict(ZINB3,type="prob")[,1]),
                  "M4"=sum(predict(ZINB4,type="prob")[,1])))

MComp.zinb<-data.frame(cbind(Logliks.zinb[c(2,4,6,8)],BICzinb[[2]],AICzinb[[2]],Logliks.zinb[c(1,3,5,7)],ZERO.zinb[2:5]))
colnames(MComp.zinb)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",ZERO.zinb[[1]],")"))
MComp.zinb$deltaAIC<-round(MComp.zinb$AIC-min(MComp.zinb$AIC),2)
MComp.zinb$deltaBIC<-round(MComp.zinb$BIC-min(MComp.zinb$BIC),2)
MComp.zinb$ModLikelihood<-round(exp(-MComp.zinb$deltaAIC/2),2)
MComp.zinb$AICweight<-round(MComp.zinb$ModLikelihood/sum(MComp.zinb$ModLikelihood),2)
MComp.zinb<-MComp.zinb[order(MComp.zinb$deltaAIC),]
MComp.zinb

#8.8# Basic model checking

wrasse.dat$residZINB<-residuals(ZINB1,type="pearson")
wrasse.dat$fitZINB<-fitted(ZINB1,type="response")

#8.8.1# Normality (not really relevant) and heterogeneity

hn<-ggplot(wrasse.dat, aes(x=residZINB))+
  geom_histogram(binwidth = 0.2,aes(fill = ..count..))+
  scale_fill_gradient("Count", low = "green", high = "red") +
  geom_vline(aes(xintercept=0),col="black")+
  scale_x_continuous("Pearson residuals",limits=c(0,20)) #limiting x range... is much more wider
frn<-ggplot(data=wrasse.dat,aes(x=fitZINB,y=residZINB,colour=residZINB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red")+
  scale_y_continuous("Pearson residuals")+
  scale_x_continuous("Fitted values")+
  geom_hline(aes(yintercept=0),col="black")+
  ggtitle("Residual vs. Fitted") # some lack of fit at low and high values

pdf(file="ZINBassumption_wrasse.pdf",width=10,height=10)

multiplot(hn,frn)

dev.off()

#8.8.2# Residuals vs predictors

rpn1<-ggplot(wrasse.dat, aes(x=factor(Idflota), y=residZINB))+
  geom_boxplot(col="gray50",aes(fill = GRT))+
  scale_fill_gradient(low = "green", high = "red") +
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Boat")
rpn2<-ggplot(wrasse.dat, aes(x=fyear, y=residZINB))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Year")+
  facet_wrap(~fZoneO,scales="free")
rpn3<-ggplot(wrasse.dat, aes(x=fyear, y=residZINB))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Year")+
  facet_wrap(~Seafloor,scales="free")
rpn4<-ggplot(wrasse.dat, aes(x=Julian, y=residZINB,colour=fitZINB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Julian")
rpn5<-ggplot(wrasse.dat, aes(x=Depth, y=residZINB,colour=fitZINB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Depth(m)")
rpn6<-ggplot(wrasse.dat, aes(x=caladoNight, y=residZINB,colour=fitZINB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("% Night time")

pdf(file="ZINBresid-expl_wrasse.pdf",width=10,height=16)

multiplot(rpn2,rpn4,rpn1,rpn3,rpn5,rpn6,cols=2)

dev.off()

#8.8.3# Spatial patterns

ggplot(data=wrasse.dat,aes(x=Lon,y=Lat))+
  geom_point(aes(colour=residZINB))+
  scale_colour_gradient(low="blue",high="red")

wrasse.dat.Spat<-data.frame(wrasse.dat$residZINB,wrasse.dat$Lon,wrasse.dat$Lat)
coordinates(wrasse.dat.Spat)<- ~wrasse.dat.Lon+wrasse.dat.Lat
vario1<-variogram(wrasse.dat.residZINB~1,data=wrasse.dat.Spat)
plot(vario1,pch=16,col=1,cex=1.5) # No spatial autocorrelation using the default
# classical method of moments estimate of the variogram (this assumes normality)
vario2<-variogram(wrasse.dat.residZINB~1,data=wrasse.dat.Spat,cressie=TRUE)
plot(vario2,pch=16,col=1,cex=1.5) # No clear spatial autocorrelation using the
# robust variogram estimate

var1zinb<-ggplot(data=vario1,aes(x=dist,y=gamma))+
  geom_point(aes(size=np),col="#3182bd")+
  scale_y_continuous("Sample variogram",limits=c(0,3.5))+
  scale_x_continuous("Distance")+
  scale_size_continuous("Number of \npoint pairs")+
  ggtitle("Classical variogram")+
  theme(legend.position=c(0,0),legend.justification=c(-0.3,-0.15))
var1zinb

var2zinb<-ggplot(data=vario2,aes(x=dist,y=gamma))+
  geom_point(aes(size=np),col="#3182bd")+
  scale_y_continuous("Sample variogram",limits=c(0,0.2))+
  scale_x_continuous("Distance")+
  scale_size_continuous("Number of \npoint pairs")+
  ggtitle("Robust variogram")+
  theme(legend.position="none")
var2zinb

pdf(file="ZINBresid-wrasse.pdf",width=13,height=8)

multiplot(resPlot1,var1zinb,var2zinb,cols=3)

dev.off()

#8.9# Model coefficients interpretation (care with ln-transformation of Depth)??

exp(coef(ZINB1)[c(1:61)]) #count part
exp(coef(ZINB1)[c(62:65)]) #zero part

#8.10# different data frame for the analysis
#awful model checking
#we are no able to catch zero counts

#try model only with MINOS that present the lowest number of zeros
dat<-wrasse.dat[wrasse.dat$Gear=="MINOS",]
trial<-zeroinfl(Ntot~offset(log(OffSet))+
                  lGRT+
                  fyear*fZoneO+
                  poly(Julian,2)+
                  lDepth+
                  sstM+
                  caladoNight+
                  Seafloor|lDepth+Seafloor,
                dist="negbin",link="logit",data=dat)
summary(trial)
sum(residuals(trial,type="pearson")^2)/(nrow(dat)-(length(coef(trial))+1))


#ploting comparison
dat$residZINB<-residuals(trial,type="pearson")
dat$fitZINB<-fitted(trial,type="response")

hn2<-ggplot(dat, aes(x=residZINB))+
  geom_histogram(binwidth = 0.2,aes(fill = ..count..))+
  scale_fill_gradient("Count", low = "green", high = "red") +
  geom_vline(aes(xintercept=0),col="black")+
  scale_x_continuous("Pearson residuals") #limiting x range... is much more wider
frn2<-ggplot(data=dat,aes(x=fitZINB,y=residZINB,colour=residZINB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red")+
  scale_y_continuous("Pearson residuals")+
  scale_x_continuous("Fitted values")+
  geom_hline(aes(yintercept=0),col="black")+
  ggtitle("Residual vs. Fitted") # some lack of fit at low and high values

multiplot(hn,hn2,frn,frn2,cols=2)

rpn1b<-ggplot(dat, aes(x=factor(Idflota), y=residZINB))+
  geom_boxplot(col="gray50",aes(fill = GRT))+
  scale_fill_gradient(low = "green", high = "red") +
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Boat")
rpn2b<-ggplot(dat, aes(x=fyear, y=residZINB))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Year")+
  facet_wrap(~fZoneO,scales="free")
rpn3b<-ggplot(dat, aes(x=fyear, y=residZINB))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Year")+
  facet_wrap(~Seafloor,scales="free")
rpn4b<-ggplot(dat, aes(x=Julian, y=residZINB,colour=fitZINB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Julian")
rpn5b<-ggplot(dat, aes(x=Depth, y=residZINB,colour=fitZINB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Depth(m)")
rpn6b<-ggplot(dat, aes(x=caladoNight, y=residZINB,colour=fitZINB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("% Night time")

multiplot(rpn2,rpn2b,rpn4,rpn4b,rpn1,rpn1b,rpn3,rpn3b,rpn5,rpn5b,rpn6,rpn6b,cols=3)

ggplot(data=dat,aes(x=Lon,y=Lat))+
  geom_point(aes(colour=residZINB))+
  scale_colour_gradient(low="blue",high="red")

dat.Spat<-data.frame(dat$residZINB,dat$Lon,dat$Lat)
coordinates(dat.Spat)<- ~dat.Lon+dat.Lat
vario1b<-variogram(dat.residZINB~1,data=dat.Spat)
plot(vario1b,pch=16,col=1,cex=1.5) # No spatial autocorrelation using the default
# classical method of moments estimate of the variogram (this assumes normality)
vario2b<-variogram(wrasse.dat.residZINB~1,data=wrasse.dat.Spat,cressie=TRUE)
plot(vario2b,pch=16,col=1,cex=1.5) # No clear spatial autocorrelation using the
# robust variogram estimate

var1zinb2<-ggplot(data=vario1b,aes(x=dist,y=gamma))+
  geom_point(aes(size=np),col="#3182bd")+
  scale_y_continuous("Sample variogram",limits=c(0,3))+
  scale_x_continuous("Distance")+
  scale_size_continuous("Number of \npoint pairs")+
  ggtitle("Classical variogram")+
  theme(legend.position=c(0,0),legend.justification=c(-0.3,-0.15))

var2zinb2<-ggplot(data=vario2b,aes(x=dist,y=gamma))+
  geom_point(aes(size=np),col="#3182bd")+
  scale_y_continuous("Sample variogram",limits=c(0,0.2))+
  scale_x_continuous("Distance")+
  scale_size_continuous("Number of \npoint pairs")+
  ggtitle("Robust variogram")+
  theme(legend.position="none")

multiplot(hn,hn2,var1zinb,var1zinb2,frn,frn2,var2zinb,var2zinb2,cols=4)

#we have lost number of samples but it seems that the model improve a bit
#so, start over again selection of predictors

#8.11# alternative model only with Miños
#poisson
minos.poiss1<-glm(Ntot~offset(log(OffSet))+
                    lGRT+
                    fyear*fZoneO+
                    poly(Julian,2)+
                    lDepth+
                    sstM+
                    caladoNight+
                    Seafloor,
                  family=poisson,data=dat)

summary(minos.poiss1)
Anova(minos.poiss1)
plot(allEffects(minos.poiss1))
sum(residuals(minos.poiss1,type="pearson")^2)/minos.poiss1$df.resid # Overdispersion

#Negative Binomial
minos.negbin1<-glm.nb(Ntot~offset(log(OffSet))+
                        lGRT+
                        fyear*fZoneO+
                        poly(Julian,2)+
                        lDepth+
                        sstM+
                        caladoNight+
                        Seafloor,
                      data=dat)

summary(minos.negbin1)
Anova(minos.negbin1)
plot(allEffects(minos.negbin1))
sum(residuals(minos.negbin1,type="pearson")^2)/(nrow(dat)-(length(coef(minos.negbin1))+1))

odTest(minos.negbin1)

#ZERO INFLATED
#exploring binomial 
dat$binNtot<-ifelse(dat$Ntot==0,0,1)

minos.bin0<-glm(binNtot~lDepth+sstM+caladoNight+Seafloor,
                family=binomial,data=dat)
summary(minos.bin0)
anova(minos.bin0)
plot(allEffects(minos.bin0))

#ZI Models
minos.zipoiss1<-zeroinfl(Ntot~offset(log(OffSet))+
                           lGRT+
                           fyear*fZoneO+
                           poly(Julian,2)+
                           lDepth+
                           sstM+
                           caladoNight+
                           Seafloor|1,
                         dist="poisson",link="logit",data=dat)

summary(minos.zipoiss1)
sum(residuals(minos.zipoiss1,type="pearson")^2)/(nrow(dat)-(length(coef(minos.zipoiss1))))

minos.zipoiss2<-zeroinfl(Ntot~offset(log(OffSet))+
                           lGRT+
                           fyear*fZoneO+
                           poly(Julian,2)+
                           lDepth+
                           sstM+
                           caladoNight+
                           Seafloor|lDepth+Seafloor,
                         dist="poisson",link="logit",data=dat)

summary(minos.zipoiss2)
sum(residuals(minos.zipoiss2,type="pearson")^2)/(nrow(dat)-(length(coef(minos.zipoiss2))))

minos.zinb1<-zeroinfl(Ntot~offset(log(OffSet))+
                        lGRT+
                        fyear*fZoneO+
                        poly(Julian,2)+
                        lDepth+
                        sstM+
                        caladoNight+
                        Seafloor|1,
                      dist="negbin",link="logit",data=dat)

summary(minos.zinb1)
sum(residuals(minos.zinb1,type="pearson")^2)/(nrow(dat)-(length(coef(minos.zinb1))+1))

minos.zinb2<-zeroinfl(Ntot~offset(log(OffSet))+
                        lGRT+
                        fyear*fZoneO+
                        poly(Julian,2)+
                        lDepth+
                        sstM+
                        caladoNight+
                        Seafloor|lDepth+Seafloor,
                      dist="negbin",link="logit",data=dat)

summary(minos.zinb2)
sum(residuals(minos.zinb2,type="pearson")^2)/(nrow(dat)-(length(coef(minos.zinb2))+1))

#Mixed model NB it does not estimate standar errors

#model comparison
# AIC/BIC
aic.minos=AIC(minos.poiss1,minos.negbin1,minos.zipoiss1,minos.zipoiss2,minos.zinb1,minos.zinb2,k=2)
aic.minos[order(aic.minos$AIC), ]
bic.minos=AIC(minos.poiss1,minos.negbin1,minos.zipoiss1,minos.zipoiss2,minos.zinb1,minos.zinb2,k=log(dim(wrasse.dat)[1]))
bic.minos[order(bic.minos$AIC), ]

#overdispersion parameter
minos.GLMpoisson=sum(residuals(minos.poiss1,type="pearson")^2)/minos.poiss1$df.resid
minos.GLMNB=sum(residuals(negbin1,type="pearson")^2)/(nrow(wrasse.dat)-(length(coef(negbin1))+1)) 
minos.ZIpoisson1=sum(residuals(minos.zipoiss1,type="pearson")^2)/(nrow(dat)-(length(coef(minos.zipoiss1))))
minos.ZIpoisson2=sum(residuals(minos.zipoiss2,type="pearson")^2)/(nrow(dat)-(length(coef(minos.zipoiss2))))
minos.ZINB1=sum(residuals(minos.zinb1,type="pearson")^2)/(nrow(dat)-(length(coef(minos.zinb1))+1))
minos.ZINB2=sum(residuals(minos.zinb2,type="pearson")^2)/(nrow(dat)-(length(coef(minos.zinb2))+1))

minos.model=c("minos.GLMpoisson","minos.GLMNB","minos.ZIpoisson1","minos.ZIpoisson2","minos.ZINB1","minos.ZINB2")
minos.theta=c(minos.GLMpoisson,minos.GLMNB,minos.ZIpoisson1,minos.ZIpoisson2,minos.ZINB1,minos.ZINB2)
minos.M=data.frame(cbind(minos.model,minos.theta))
minos.M$minos.theta<-round(as.numeric(levels(minos.M$minos.theta))[minos.M$minos.theta],digits = 2)
minos.M=minos.M[order(minos.M$minos.theta),]
minos.M

#loglikelihood
fm.minos<-list("GLMPois"=minos.poiss1,
               "GLMNB"=minos.negbin1,
               "ZIPOISS1"=minos.zipoiss1,"ZIPOISS2"=minos.zipoiss2,
               "ZINB1"=minos.zinb1,"ZINB2"=minos.zinb2)

logliks.minos<-rbind(logLik=sapply(fm.minos,function (x) round(logLik(x),digits=0)),
                     Df=sapply(fm.minos,function (x) attr(logLik(x),"df")))
logliks.minos

#vuong test
vuong(minos.poiss1,minos.negbin1)
vuong(minos.negbin1,minos.zinb1)
vuong(minos.zinb1,minos.zinb2)

#8.6.5# Predicting probabilities

minos.phat.pois<-predprob(minos.poiss1)
minos.phat.pois.mn<-apply(minos.phat.pois,2,mean)
minos.phat.nb<-predprob(minos.negbin1)
minos.phat.nb.mn<-apply(minos.phat.nb,2,mean)
minos.phat.zipoiss1<-predprob(minos.zipoiss1)
minos.phat.zipoiss.mn1<-apply(minos.phat.zipoiss1,2,mean)
minos.phat.zipoiss2<-predprob(minos.zipoiss2)
minos.phat.zipoiss.mn2<-apply(minos.phat.zipoiss2,2,mean)
minos.phat.zinb1<-predprob(minos.zinb1)
minos.phat.zinb.mn1<-apply(minos.phat.zinb1,2,mean)
minos.phat.zinb2<-predprob(minos.zinb2)
minos.phat.zinb.mn2<-apply(minos.phat.zinb2,2,mean)

with(wrasse.dat,{
  hist(Ntot,prob=TRUE,col="grey",breaks=seq(-0.5,1450,1),xlab="",
       main="",ylim=c(0,1),xlim=c(-0.5,30))
  lines(x=seq(0,173,1),y=minos.phat.pois.mn,type="l",lwd=2,col=1) #negro
  lines(x=seq(0,173,1),y=minos.phat.nb.mn,type="l",lwd=2,col=2) #rojo
  lines(x=seq(0,173,1),y=minos.phat.zipoiss.mn1,type="l",lwd=2,col=3) #verde
  lines(x=seq(0,173,1),y=minos.phat.zipoiss.mn2,type="l",lwd=2,col=4) #azul
  lines(x=seq(0,173,1),y=minos.phat.zinb.mn1,type="l",lwd=2,col=5) #turquesa
  lines(x=seq(0,173,1),y=minos.phat.zinb.mn2,type="l",lwd=2,col=6) #rosa
})

#Comparison of models
zeroes.minos<-round(c("Obs"=sum(dat$Ntot<1),
                      "GLMPois"=sum(dpois(0,fitted(minos.poiss1))),
                      "GLMNB"=sum(dnbinom(0,mu=fitted(minos.negbin1),size=minos.negbin1$theta)),
                      "ZIPOISS1"=sum(predict(minos.zipoiss1,type="prob")[,1]),
                      "ZIPOISS2"=sum(predict(minos.zipoiss2,type="prob")[,1]),
                      "ZINB1"=sum(predict(minos.zinb1,type="prob")[,1]),
                      "ZINB2"=sum(predict(minos.zinb2,type="prob")[,1])))

modelsComp.minos<-data.frame(cbind(logliks.minos[c(2,4,6,8,10,12)],bic.minos[[2]],aic.minos[[2]],logliks.minos[c(1,3,5,7,9,11)],zeroes.minos[2:7]))
colnames(modelsComp.minos)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes[[1]],")"))
modelsComp.minos$deltaAIC<-round(modelsComp.minos$AIC-min(modelsComp.minos$AIC),2)
modelsComp.minos$deltaBIC<-round(modelsComp.minos$BIC-min(modelsComp.minos$BIC),2)
modelsComp.minos$ModLikelihood<-round(exp(-modelsComp.minos$deltaAIC/2),2)
modelsComp.minos$AICweight<-round(modelsComp.minos$ModLikelihood/sum(modelsComp.minos$ModLikelihood),2)
modelsComp.minos<-modelsComp.minos[order(modelsComp.minos$deltaAIC),]
modelsComp.minos

#select model predictors
minos.ZINB1<-minos.zinb2
summary(minos.ZINB1)

#year:zone
minos.ZINB2<-zeroinfl(Ntot~offset(log(OffSet))+
                        lGRT+
                        fyear+fZoneO+
                        poly(Julian,2)+
                        lDepth+
                        sstM+
                        caladoNight+
                        Seafloor|lDepth+Seafloor,
                      dist="negbin",link="logit",data=dat)
summary(minos.ZINB2)
lrtest(minos.ZINB1,minos.ZINB2)
vuong(minos.ZINB1,minos.ZINB2)
AICtab(minos.ZINB1,minos.ZINB2)

#sst (reintroducing interaction)
minos.ZINB3<-zeroinfl(Ntot~offset(log(OffSet))+
                        lGRT+
                        fyear*fZoneO+
                        poly(Julian,2)+
                        lDepth+
                        caladoNight+
                        Seafloor|lDepth+Seafloor,
                      dist="negbin",link="logit",data=dat)
summary(minos.ZINB3)
lrtest(minos.ZINB2,minos.ZINB3)
lrtest(minos.ZINB1,minos.ZINB3)
vuong(minos.ZINB2,minos.ZINB3)
vuong(minos.ZINB1,minos.ZINB3)
AICtab(minos.ZINB2,minos.ZINB3)
AICtab(minos.ZINB1,minos.ZINB3)

#sstM (reintroducing interaction)
minos.ZINB4<-zeroinfl(Ntot~offset(log(OffSet))+
                        lGRT+
                        fyear*fZoneO+
                        poly(Julian,2)+
                        lDepth+
                        Seafloor|lDepth+Seafloor,
                      dist="negbin",link="logit",data=dat)
summary(minos.ZINB4)
lrtest(minos.ZINB3,minos.ZINB4)
lrtest(minos.ZINB1,minos.ZINB4)
vuong(minos.ZINB3,minos.ZINB4)
vuong(minos.ZINB1,minos.ZINB4)
AICtab(minos.ZINB3,minos.ZINB4)
AICtab(minos.ZINB1,minos.ZINB4)

sum(residuals(minos.ZINB3,type="pearson")^2)/(nrow(dat)-(length(coef(minos.ZINB3))+1))
#slightly overdispersed, so, this can lead to inflated significance

#8.7.2# comparing models
aic2<-AIC(minos.ZINB1,minos.ZINB2,minos.ZINB3,minos.ZINB4,k=2) # AIC
bic2<-AIC(minos.ZINB1,minos.ZINB2,minos.ZINB3,minos.ZINB4,k=log(dim(dat)[1])) # BIC

fm2<-list("M1"=minos.ZINB1,"M2"=minos.ZINB2,"M3"=minos.ZINB3,"M4"=minos.ZINB4)
Logliks2<-rbind(logLik=sapply(fm2,function (x) round(logLik(x),digits=0)),
                Df=sapply(fm2,function (x) attr(logLik(x),"df")))

ZERO2<-round(c("Obs"=sum(dat$Ntot<1),
               "M1"=sum(predict(minos.ZINB1,type="prob")[,1]),
               "M2"=sum(predict(minos.ZINB2,type="prob")[,1]),
               "M3"=sum(predict(minos.ZINB3,type="prob")[,1]),
               "M4"=sum(predict(minos.ZINB4,type="prob")[,1])))

MComp2<-data.frame(cbind(Logliks2[c(2,4,6,8)],bic2[[2]],aic2[[2]],Logliks2[c(1,3,5,7)],ZERO2[2:5]))
colnames(MComp2)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",ZERO2[[1]],")"))
MComp2$deltaAIC<-round(MComp2$AIC-min(MComp2$AIC),2)
MComp2$deltaBIC<-round(MComp2$BIC-min(MComp2$BIC),2)
MComp2$ModLikelihood<-round(exp(-MComp2$deltaAIC/2),2)
MComp2$AICweight<-round(MComp2$ModLikelihood/sum(MComp2$ModLikelihood),2)
MComp2<-MComp2[order(MComp2$deltaAIC),]
MComp2

summary(minos.ZINB3)
sum(residuals(minos.ZINB3,type="pearson")^2)/(nrow(dat)-(length(coef(minos.ZINB3))+1))
#still slightly overdispersed

# ----------------------------- #
#9# Bootstrapping the optimal model coefficients (zero & count parts)

#9.1# Function (add starting values to the model if needed!!)

dput(round(coef(minos.ZINB3,"count"),4))
dput(round(coef(minos.ZINB3,"zero"),4))

length(coef(minos.ZINB3))

boot.zinb<-function (data,i) {
  
  try(mod<-zeroinfl(Ntot~offset(log(OffSet))+
                      lGRT+
                      fyear*fZoneO+
                      poly(Julian,2)+
                      lDepth+
                      caladoNight+
                      Seafloor|lDepth+Seafloor,
                    dist="negbin",link="logit",data=data[i,]))
  
  if (exists("mod")) { 
    
    as.vector(t(do.call(rbind,coef(summary(mod)))[,1:2]))
    
  } else {rep(NA,times=114)}
  
}

#9.2# Coefficients (obtain CI of estimates excluding SE and theta)

RR<-10 # Number of resamples (reduce this number for testing the code)
zinb.boot.out<-boot(data=dat,statistic=boot.zinb,R=RR)
zinb.boot.out  # Basic output
plot(zinb.boot.out,index=1)

zinb.boot.out2<-as.data.frame(zinb.boot.out$t[,c(seq(1,103,2),107,109,111,113)]) # Coefficients of interest from the boot object matrix
colnames(zinb.boot.out2)<-names(coef(minos.ZINB3))
head(zinb.boot.out2)
write.table(x=zinb.boot.out2,file="wrassebootCoefs.txt",row.names=F)

parmsPERC<-t(sapply(c(seq(1,103,2),107,109,111,113),function (i) {
  out<-boot.ci(zinb.boot.out,index=i,type=c("perc"))
  with(out,c(Estimate=t0,pLow=percent[4],pUpp=percent[5]))
})) # Obtain intervals calculated using the bootstrap percentile method
row.names(parmsPERC)<-names(coef(minos.ZINB3))
head(parmsPERC)
write.table(x=as.data.frame(parmsPERC),file="wrassebootPERC.txt",row.names=T)

#BCA 

ci<-data.frame(cbind(parmsPERC,confint(minos.ZINB3)))
colnames(ci)<-c("Estimate","pLow","pUpp","2.5%","97.%")
str(ci)
ci
write.table(ci,file="ci_wrasse.txt",row.names=T)

#9.3# Example of histogram and qqplot for a given component

plot(zinb.boot.out,index=1)

#9.4# Histograms of all components

pdf(file="hist_bootwrasse.pdf",width=20,height=23)

par(mfrow=c(8,7))
for (i in 1:56) {
  hist(zinb.boot.out2[,i],breaks=50,col="light blue",main="",
       xlab=names(coef(minos.ZINB3))[i])
  abline(v=coef(minos.ZINB3)[i],col="red",lwd=2)
}

dev.off()


# ----------------------------- #
#10# Sketching results for the optimal model

#10.1# plotting nasa model
names(coef(minos.ZINB3))

rugsGRT<-data.frame(rugsGRT=unique(dat$GRT))
rugscalN<-data.frame(rugscalN=unique(dat$caladoNight))
rugsDoY<-data.frame(rugsDoY=unique(dat$Julian))
rugsDepth<-data.frame(rugsDepth=unique(dat$Depth))

#10.1.1#ZI part
prezero1<-predict(minos.ZINB3,type="zero",
                 newdata=data.frame(lGRT=mean(dat$lGRT,na.rm=T), # lGRT
                                    fyear="2006",
                                    Julian=183,
                                    lDepth=seq(min(dat$lDepth,na.rm=T),
                                               max(dat$lDepth,na.rm=T),length=100),
                                    caladoNight=mean(dat$caladoNight,na.rm=T),
                                    Seafloor="hard",
                                    fZoneO="1",
                                    OffSet=2500)) #2500 area media aprox. de miños

prezero1<-data.frame(prezero1,
                    seq(min(dat$lDepth),max(dat$lDepth),length.out = 100))
colnames(prezero1)<-c("ZeroPre","DepthSeq")
prezero1$Depth<-exp(prezero1$DepthSeq)
str(prezero1)
summary(prezero1)
prezero1

zeroDepth<-ggplot(data=prezero1,aes(x=Depth,y=ZeroPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Probability of false zeros",limits=c(-0.02,1))+
  scale_x_continuous("Depth (m)",limits=c(0,605))+
  geom_segment(data=rugsDepth,aes(x=rugsDepth,xend=rugsDepth,y=-0.02,yend=-0.01),stat="identity",lwd=0.1,col="gray50")
zeroDepth

prezero2<-predict(minos.ZINB3,type="zero",
                  newdata=data.frame(lGRT=rep(mean(dat$lGRT,na.rm=T),times=3), # lGRT
                                     fyear=rep("2006",times=3),
                                     Julian=rep(183,times=3),
                                     lDepth=rep(mean(dat$lDepth,na.rm=T),times=3),
                                     caladoNight=rep(mean(dat$caladoNight,na.rm=T),times=3),
                                     Seafloor=c("hard","mixed","soft"),
                                     fZoneO=rep("1",times=3),
                                     OffSet=2500)) #2500 area media aprox. de miños

prezero2<-data.frame(prezero2,c("hard","mixed","soft"))
colnames(prezero2)<-c("SeafloorPre","SeafloorSeq")
str(prezero2)
summary(prezero2)
prezero2

zeroSeafloor<-ggplot(prezero2, aes(x = factor(SeafloorSeq), y = SeafloorPre))+
  geom_bar(stat = "identity",fill="blue",col="grey")+
  scale_y_continuous("Standardized Index",limits=c(0,1))+
  scale_x_discrete("Sea-floor type",labels=c("Hard","Mixed","Soft"))
zeroSeafloor

pdf(file="wrasse_preZIzero.pdf",width=10,height=8)

multiplot(zeroDepth,zeroSeafloor,cols=2)

dev.off()

#10.1.2# Continuous variables

minos1<-predict(minos.ZINB3,type="count",
                 newdata=data.frame(lGRT=seq(min(dat$lGRT,na.rm=T),
                                             max(dat$lGRT,na.rm=T),length=100), # lGRT
                                    fyear="2006",
                                    Julian=183,
                                    lDepth=mean(dat$lDepth,na.rm=T),
                                    caladoNight=mean(dat$caladoNight,na.rm=T),
                                    Seafloor="hard",
                                    fZoneO="1",
                                    OffSet=2500)) #2500 area media aprox. de miños

datGRT<-data.frame(minos1,seq(min(dat$lGRT,na.rm=T),max(dat$lGRT,na.rm=T),length=100))
colnames(datGRT)<-c("GRTPre","GRTSeq")
datGRT$GRT<-exp(datGRT$GRTSeq)
head(datGRT)
summary(datGRT)

gGRT<-ggplot(data=datGRT,aes(x=GRT,y=GRTPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(0.25,1.5))+
  scale_x_continuous("GRT (Tons)")+
  #((1.5-0.25)*0.01)+0.25
  geom_segment(data=rugsGRT,aes(x=rugsGRT,xend=rugsGRT,y=0.25,yend=0.2625),stat="identity",lwd=0.5,col="gray50")
gGRT

minos2<-predict(minos.ZINB3,type="count",
                newdata=data.frame(lGRT=mean(dat$lGRT,na.rm=T), # lGRT
                                   fyear="2006",
                                   Julian=seq(1,365,1),
                                   lDepth=mean(dat$lDepth,na.rm=T),
                                   caladoNight=mean(dat$caladoNight,na.rm=T),
                                   Seafloor="hard",
                                   fZoneO="1",
                                   OffSet=2500)) #2500 area media aprox. de miños

datJulian<-data.frame(minos2,seq(1,365,1))
colnames(datJulian)<-c("JulianPre","JulianSeq")
head(datJulian)
summary(datJulian)

gJulian<-ggplot(data=datJulian,aes(x=JulianSeq,y=JulianPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(0.15,0.9))+
  scale_x_continuous("Day of the year")+
  #((0.9-0.15)*0.01)+0.15
  geom_segment(data=rugsDoY,aes(x=rugsDoY,xend=rugsDoY,y=0.15,yend=0.1575),stat="identity",lwd=0.5,col="gray50")
gJulian

minos3<-predict(minos.ZINB3,type="count",
                newdata=data.frame(lGRT=mean(dat$lGRT,na.rm=T), # lGRT
                                   fyear="2006",
                                   Julian=183,
                                   lDepth=seq(min(dat$lDepth,na.rm=T),
                                              max(dat$lDepth,na.rm=T),length=100),
                                   caladoNight=mean(dat$caladoNight,na.rm=T),
                                   Seafloor="hard",
                                   fZoneO="1",
                                   OffSet=2500)) #2500 area media aprox. de miños

datDepth<-data.frame(minos3,seq(min(dat$lDepth,na.rm=T),
                                max(dat$lDepth,na.rm=T),length=100))
colnames(datDepth)<-c("DepthPre","DepthSeq")
datDepth$Depth<-exp(datDepth$DepthSeq)
head(datDepth)
summary(datDepth)

gDepth<-ggplot(data=datDepth,aes(x=Depth,y=DepthPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(0,5))+
  scale_x_continuous("Depth (m)",limits=c(0,200))+
  #(5*0.01)
  geom_segment(data=rugsDepth,aes(x=rugsDepth,xend=rugsDepth,y=0,yend=0.05),stat="identity",lwd=0.5,col="gray50")
gDepth

minos4<-predict(minos.ZINB3,type="count",
                newdata=data.frame(lGRT=mean(dat$lGRT,na.rm=T), # lGRT
                                   fyear="2006",
                                   Julian=183,
                                   lDepth=mean(dat$lDepth,na.rm=T),
                                   caladoNight=seq(min(dat$caladoNight,na.rm=T),
                                                   max(dat$caladoNight,na.rm=T),length=100),
                                   Seafloor="hard",
                                   fZoneO="1",
                                   OffSet=2500)) #2500 area media aprox. de miños

datCal<-data.frame(minos4,seq(min(dat$caladoNight,na.rm=T),
                              max(dat$caladoNight,na.rm=T),length=100))
colnames(datCal)<-c("CalNPre","CalNSeq")
head(datCal)
summary(datCal)

gCalN<-ggplot(data=datCal,aes(x=CalNSeq,y=CalNPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(0,2))+
  scale_x_continuous("% Night time")+
  geom_segment(data=rugscalN,aes(x=rugscalN,xend=rugscalN,y=0,yend=0.02),stat="identity",lwd=0.5,col="gray50")
gCalN

pdf(file="ZINBcont-wrasse.pdf",width=10,height=10)
multiplot(gJulian,gGRT,gDepth,gCalN,cols=2)
dev.off()

#10.1.3# Categorical variables

datSeafloor<-data.frame(lGRT=rep(mean(dat$lGRT,na.rm=T),times=3), # lGRT
                                           fyear=rep("2006",times=3),
                                           Julian=rep(183,times=3),
                                           lDepth=rep(mean(dat$lDepth,na.rm=T),times=3),
                                           caladoNight=rep(mean(dat$caladoNight,na.rm=T),times=3),
                                           Seafloor=c("hard","mixed","soft"),
                                           fZoneO=rep("1",times=3),
                                           OffSet=2500) #2500 area media aprox. de miños

datSeafloor<-predict(minos.ZINB3,type="count",newdata=datSeafloor)
datSeafloor<-data.frame(cbind(datSeafloor,c("hard","mixed","soft")))
colnames(datSeafloor)<-c("SeafloorPre","SeafloorSeq")
datSeafloor$SeafloorPre<-round(as.numeric(levels(datSeafloor$SeafloorPre))[datSeafloor$SeafloorPre],digits = 2)
str(datSeafloor)
datSeafloor

gSeafloor<-ggplot(datSeafloor, aes(x = factor(SeafloorSeq), y = SeafloorPre))+
  geom_bar(stat = "identity",fill="blue",col="grey")+
  scale_y_continuous("Standardized Index",limits=c(0,1))+
  scale_x_discrete("Sea-floor type",labels=c("Hard","Mixed","Soft"))
gSeafloor

pdf(file="ZINBfitvar-wrasse.pdf",width=8,height=12)
multiplot(gJulian,gGRT,gSeafloor,gDepth,gCalN,cols=2)
dev.off()

newTrend<-data.frame(lGRT=rep(mean(dat$lGRT,na.rm=T),times=45), # lGRT
                     fyear=rep(c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013"),3),
                     Julian=rep(183,times=45),
                     lDepth=rep(mean(dat$lDepth,na.rm=T),times=45),
                     caladoNight=rep(mean(dat$caladoNight,na.rm=T),times=45),
                     Seafloor=rep("hard",times=45),
                     fZoneO=c(rep("1",times=15),rep("2",times=15),rep("3",times=15)),
                     OffSet=2500) #2500 area media aprox. de miños

abundInd<-predict(minos.ZINB3,newdata=newTrend,type="count")

abundInd<-data.frame(abundInd)
colnames(abundInd)<-c("Index")
abundInd$Year<-as.numeric(rep(seq(1999,2013,1),3))
abundInd$Zones<-c(rep("Rias Baixas",times=15),rep("Golfo Artabro",times=15),rep("Cantabrico",times=15))
abundInd$Zones <- ordered(abundInd$Zones,levels = c("Rias Baixas", "Golfo Artabro", "Cantabrico"))
str(abundInd)

pdf(file="ZINBTrend-wrasse.pdf",width=10,height=7)

ggplot(data=abundInd,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index",limits=c(0,2))+
  scale_x_continuous("Year",breaks=c(2000,2005,2010),labels=c("2000","2005","2010"))+
  facet_wrap(~Zones)

dev.off()

wrasse.abund<-lm(Index~Year*Zones,data=abundInd)
summary(wrasse.abund)
anova(wrasse.abund) #no decreasing trend, no zones diferencies

#11# Plot of abundance, nominal cpue, and landings

#11.1# Calculate nominal cpue and average for the Rias Baixas
head(dat)
dat

dat$cpue<-(dat$Ntot*2500)/dat$OffSet # Standardize at 2500 m2 lance

datRB<-dat[dat$ZoneO==1,]
datAR<-dat[dat$ZoneO==2,]
datCN<-dat[dat$ZoneO==3,]

cpues.dat.RB<-tapply(datRB$cpue,datRB$Year,mean,na.rm=T)
cpues.L.dat.RB<-tapply(datRB$cpue,datRB$Year,ci95Low)
cpues.U.dat.RB<-tapply(datRB$cpue,datRB$Year,ci95Upp)
cpues.dat.RB<-data.frame(cbind(cpues.dat.RB,cpues.L.dat.RB,cpues.U.dat.RB))
cpues.dat.RB$Zones<-rep("Rias Baixas",15)
cpues.dat.RB$Year<-as.numeric(seq(1999,2013,1))
colnames(cpues.dat.RB)<-c("Index","ciLow","ciUpp","Zones","Year")
cpues.dat.RB$ciLow<-ifelse(cpues.dat.RB$ciLow<0,0,cpues.dat.RB$ciLow)
str(cpues.dat.RB)
summary(cpues.dat.RB)

cpues.dat.AR<-tapply(datAR$cpue,datAR$Year,mean,na.rm=T)
cpues.L.dat.AR<-tapply(datAR$cpue,datAR$Year,ci95Low)
cpues.U.dat.AR<-tapply(datAR$cpue,datAR$Year,ci95Upp)
cpues.dat.AR<-data.frame(cbind(cpues.dat.AR,cpues.L.dat.AR,cpues.U.dat.AR))
cpues.dat.AR$Zones<-rep("Golfo Artabro",15)
cpues.dat.AR$Year<-as.numeric(seq(1999,2013,1))
colnames(cpues.dat.AR)<-c("Index","ciLow","ciUpp","Zones","Year")
cpues.dat.AR$ciLow<-ifelse(cpues.dat.AR$ciLow<0,0,cpues.dat.AR$ciLow)
str(cpues.dat.AR)
summary(cpues.dat.AR)

cpues.dat.CN<-tapply(datCN$cpue,datCN$Year,mean,na.rm=T)
cpues.L.dat.CN<-tapply(datCN$cpue,datCN$Year,ci95Low)
cpues.U.dat.CN<-tapply(datCN$cpue,datCN$Year,ci95Upp)
cpues.dat.CN<-data.frame(cbind(cpues.dat.CN,cpues.L.dat.CN,cpues.U.dat.CN))
cpues.dat.CN$Zones<-rep("Cantabrico",15)
cpues.dat.CN$Year<-as.numeric(seq(1999,2013,1))
colnames(cpues.dat.CN)<-c("Index","ciLow","ciUpp","Zones","Year")
str(cpues.dat.CN)
summary(cpues.dat.CN)

cpues.dat.wrasse<-data.frame(rbind(cpues.dat.RB,cpues.dat.AR,cpues.dat.CN))
str(cpues.dat.wrasse)
cpues.dat.wrasse
cpues.dat.wrasse$Zones <- ordered(cpues.dat.wrasse$Zones,levels = c("Rias Baixas", "Golfo Artabro", "Cantabrico"))

limits <- aes(ymax = ciUpp, ymin= ciLow)
ggplot(data=cpues.dat.wrasse,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_errorbar(limits, colour="blue", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year")+
  facet_wrap(~Zones)

#11.2# Official landings

landings<-read.table(file="wrasse_landings.txt",header=T,dec=".")
landings
str(landings)
zone1<-landings[,c(1,12)]
colnames(zone1)<-c("Year","Index")
zone2<-landings[,c(1,13)]
colnames(zone2)<-c("Year","Index")
zone3<-landings[,c(1,14)]
colnames(zone3)<-c("Year","Index")
zones<-rbind(zone1,zone2,zone3)
landings<-data.frame(cbind(zones,c(rep("Rias Baixas",18),rep("Golfo Artabro",18),rep("Cantabrico",18))))
colnames(landings)<-c("Year","Index","Zones")
str(landings)
head(landings)
landings$Zones <- ordered(landings$Zones,levels = c("Rias Baixas", "Golfo Artabro", "Cantabrico"))

ggplot(data=landings,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)

#11.3# Comparing abundance index, nominal cpue, landings

head(abundInd)

t1<-ggplot(data=abundInd,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)+
  ggtitle("Standardize abundance index (gillnet)")+
  theme(plot.title = element_text(size=26, face="bold"))
t1

t2<-ggplot(data=cpues.dat.wrasse,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  #geom_errorbar(limits, colour="blue", width=0)+
  scale_y_continuous("Nominal CPUE (nº/2500m2*h)")+
  scale_x_continuous("Year")+
  facet_wrap(~Zones)+  ggtitle("Nominal CPUE (gillnet)")+
  theme(plot.title = element_text(size=26, face="bold"))
t2

t3<-ggplot(data=landings,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Total landings (kg)")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)+facet_wrap(~Zones)+
  ggtitle("Oficial landings")+
  theme(plot.title = element_text(size=26, face="bold"))
t3

multiplot(t1,t2,t3,cols=1)

pdf(file="wrasse_trends.pdf",width=10,height=10)

multiplot(t1,t2,t3,cols=1)

dev.off()


#correlations trends

names(abundInd)
dim(abundInd)
names(cpues.dat.wrasse)
dim(cpues.dat.wrasse)
names(landings)
dim(landings)

landings2<-landings[landings$Year>1998&landings$Year<2014,]
dim(landings2)

#Rias Baixas
cor(abundInd[abundInd$Zones=="Rias Baixas",]$Index,
    cpues.dat.wrasse[cpues.dat.wrasse$Zones=="Rias Baixas",]$Index,method="spearman")
cor(abundInd[abundInd$Zones=="Rias Baixas",]$Index,
    landings2[landings2$Zones=="Rias Baixas",]$Index,method="spearman")
cor(cpues.dat.wrasse[cpues.dat.wrasse$Zones=="Rias Baixas",]$Index,
    landings2[landings2$Zones=="Rias Baixas",]$Index,method="spearman")

#Artabro
cor(abundInd[abundInd$Zones=="Golfo Artabro",]$Index,
    cpues.dat.wrasse[cpues.dat.wrasse$Zones=="Golfo Artabro",]$Index,method="spearman")
cor(abundInd[abundInd$Zones=="Golfo Artabro",]$Index,
    landings2[landings2$Zones=="Golfo Artabro",]$Index,method="spearman")
cor(cpues.dat.wrasse[cpues.dat.wrasse$Zones=="Golfo Artabro",]$Index,
    landings2[landings2$Zones=="Golfo Artabro",]$Index,method="spearman")

#Cantabrico
cor(abundInd[abundInd$Zones=="Cantabrico",]$Index,
    cpues.dat.wrasse[cpues.dat.wrasse$Zones=="Cantabrico",]$Index,method="spearman")
cor(abundInd[abundInd$Zones=="Cantabrico",]$Index,
    landings2[landings2$Zones=="Cantabrico",]$Index,method="spearman")
cor(cpues.dat.wrasse[cpues.dat.wrasse$Zones=="Cantabrico",]$Index,
    landings2[landings2$Zones=="Cantabrico",]$Index,method="spearman")

#ridiculous correlation with oficial landings :(

# ----------------------------- #
#12# Plotting ZINB model predicted means for YEAR with Bootstrapping
# of predictions (95% CI) by means of shuffling residuals (Thierry Onkelinx code)

newTrend2<-data.frame(lGRT=rep(mean(dat$lGRT,na.rm=T),times=45), # lGRT
                     fyear=rep(c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013"),3),
                     Julian=rep(183,times=45),
                     lDepth=rep(mean(dat$lDepth,na.rm=T),times=45),
                     caladoNight=rep(mean(dat$caladoNight,na.rm=T),times=45),
                     Seafloor=rep("hard",times=45),
                     fZoneO=c(rep("1",times=15),rep("2",times=15),rep("3",times=15)),
                     OffSet=2500) #2500 area media aprox. de miños

head(newTrend2)

Fit<-predict(minos.ZINB3,type="response")

Pearson<-residuals(minos.ZINB3,type="pearson") # Pearson residuals
VarComp<-residuals(minos.ZINB3,type="response")/Pearson # Raw residuals/Pearson residuals

lGRT<-dat$lGRT
fyear<-dat$fyear
Julian<-dat$Julian
lDepth<-dat$lDepth
caladoNight<-dat$caladoNight
fZoneO<-dat$fZoneO
Seafloor<-dat$Seafloor
OffSet<-dat$OffSet

#bootstraping residuals
RR<-100 # Number of resamples (reduce this number for testing the code)

bootstrap<-replicate(n=RR,{ 
  
  yStar<-pmax(round(Fit+sample(Pearson)*VarComp,0),0)
  try(mod<-zeroinfl(yStar~offset(log(OffSet))+
                      lGRT+
                      fyear*fZoneO+
                      poly(Julian,2)+
                      lDepth+
                      caladoNight+
                      Seafloor|lDepth+Seafloor,
                    dist="negbin",link="logit"))
  
  if (exists("mod")) {
    
    predict(mod,newdata=newTrend2,type="response")
    
  } else {rep(NA,times=45)} # If the above model crashes this fills in the gaps
  # with NA and the algorithm continues
  
})

CIs<-t(apply(X=bootstrap,MARGIN=1,FUN=quantile,c(0.025,0.975),na.rm=T))
colnames(CIs)<-c("ciLow","ciUpp")
newTrend2$fit<-predict(minos.ZINB3,newdata=newTrend2,type="response")
newTrend2<-cbind(newTrend2,CIs)
newTrend2$Year<-seq(1999,2013,1)
names(newTrend2)
str(newTrend2)
newTrend2

#12.2# Plot of abundance

levels(newTrend2$fZoneO)<-c("1"="Rias Baixas", "2"="Golfo Artabro","3"="Cantabrico")
levels(newTrend2$fZoneO)

pdf(file="wrasse_YearPredCI.pdf",width=10,height=8)

limits <- aes(ymax = ciUpp, ymin= ciLow)
ggplot(data=newTrend2,aes(x=Year,y=fit))+
  #geom_segment(aes(x=Year,y=fit,xend=Year,yend=ciUpp),col="black",lwd=0.5)+
  #geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=fit),col="black",lwd=0.5)+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_errorbar(limits, colour="black", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Standardized Index of Abundance")+
  scale_x_continuous("Year")+
  facet_wrap(~fZoneO,scales="free")

dev.off()

boot1<-ggplot(data=newTrend2,aes(x=Year,y=fit))+
  #geom_segment(aes(x=Year,y=fit,xend=Year,yend=ciUpp),col="black",lwd=0.5)+
  #geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=fit),col="black",lwd=0.5)+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_errorbar(limits, colour="black", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Standardized Index of Abundance")+
  scale_x_continuous("Year")+
  facet_wrap(~fZoneO,scales="free")+
  ggtitle("Standardized abundance index")+
  theme(plot.title = element_text(size=26, face="bold"))
boot1

cpues.dat.wrasse
limits <- aes(ymax = ciUpp, ymin= ciLow)
cpue1<-ggplot(data=cpues.dat.wrasse,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_errorbar(limits, colour="black", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Nominal CPUE (nº/2500m2*h)")+
  scale_x_continuous("Year")+
  facet_wrap(~Zones,scales="free_y")+
  ggtitle("Nominal CPUE")+
  theme(plot.title = element_text(size=26, face="bold"))
cpue1

landing1<-ggplot(data=landings,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Total landings (kg)")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones,scales="free_y")+
  ggtitle("Oficial landings")+
  theme(plot.title = element_text(size=26, face="bold"))
landing1

multiplot(boot1,cpue1,landing1,cols=1)

pdf(file="wrasse-ci_trends.pdf",width=15,height=20)

multiplot(boot1,cpue1,landing1,cols=1)

dev.off()
