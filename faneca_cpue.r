##########################################################
# Analysis of Trisopterus luscus catch rates in Galicia#
# using data sampled onboard fishing vessels by the UTPB #
##########################################################

#Date start: April 2014 by Jaime Otero (bulk of code from pollachius analysis)
#Last update: 04-02-2015 by Alex Alonso
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

source(file="C:\\Users\\alex\\Dropbox\\casgass\\HighStatLib.r")
source(file="C:\\Users\\alex\\Dropbox\\casgass\\multiple_ggplot.r")
source(file="C:\\Users\\alex\\Dropbox\\casgass\\CI_mean.r")

#set working directory for input and output
setwd("D:\\iim-csic\\proyectos\\ICES_CASGASS\\analysis\\pouting")

#to start since the last point
load("pouting.RData")

#1# Load species data
#el código está hecho en base al análisis de pollachius
#adaptar el código a cada especies
pouting<-read.table(file="faneca.txt",header=T,dec=".",sep=",")

head(pouting)
dim(pouting)
length(unique(pouting$Idflota)) # Number of vessels surveyed


# ----------------------------- #
#2# Evaluate fishing per gear, distribution of response and other elementary information

faneca.tot<-sum(pouting$Ntot,na.rm=T) # Total number of caught pouting
fanecas<-as.numeric(tapply(pouting$Ntot,pouting$Gear,sum,na.rm=T)) # Total number per gear
gear.hauls<-tapply(pouting$Ntot,pouting$Gear,length) # Hauls per gear
fish<-ifelse(pouting$Ntot==0,0,1) # Variable to classify 0/1 hauls
zero.fish<-gear.hauls-tapply(fish,pouting$Gear,sum,na.rm=T) # Zeros per gear
gear.zeros<-round((zero.fish*100)/gear.hauls,2) # Percentage of zeros per gear

basic.info<-data.frame(cbind(gear.hauls,gear.zeros,(fanecas*100)/faneca.tot))
colnames(basic.info)<-c("hauls","zeros","catch")
basic.info$gear<-c("Boliche","Bou de Man","Bou de Vara","Cerco bolo","Liña cordel",
	"Miños","Nasa centola","Nasa choco","Nasa nécora",
	"Nasa patexo","Nasa peixe","Nasa pulpo","Voitron","Palangrillo","Racu",
  "Rastro camaron","Rastro vieira","Trasmallo","Veta","Xeito")

pdf(file="C:\\Users\\alex\\Dropbox\\casgass\\gearsPouting.pdf",width=15,height=8)

ggplot(data=basic.info,aes(x=gear,y=catch))+
	geom_bar(stat="identity")+
	scale_y_continuous(limits=c(0,60),"Frequency")+
	scale_x_discrete("Gear")+
	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=14),
		  axis.text.y=element_text(size=14))+
	geom_text(aes(label=paste("[",hauls,",",zeros,"]"),vjust=-0.5),size=3)

dev.off()

par(mar=c(5,5,3,3))
plot(table(pouting$Ntot),ylab="Frequency",xlab="Number of pouting caught per haul")
par(mar=c(5,5,3,3))
plot(table(pouting$Ntot),ylab="Frequency",xlab="Number of pouting caught per haul",xlim=c(0,10))

basic.info

basic.info[order(-basic.info$catch), ]
basic.info[order(basic.info$zeros), ]
#seleccionamos Nasa peixes ("fanequeira") y vetas y miños en base
#al porcentage de capturas muestreadas
#y al porcentaje de zero hauls
#entre las dos representan casi el 70% y tienen porcentages de zeros <50%

# ----------------------------- #
#3# Select data set that only include the 2 main gears
pouting1<-pouting[pouting$Gear=="NASA-PEIXES" | pouting$Gear=="VETAS"| pouting$Gear=="MINOS",]

head(pouting1)
dim(pouting1)
length(unique(pouting1$Idflota))
length(unique(pouting1$Idlance))
length(unique(pouting1$Idjornada))


# ----------------------------- #
#4# Add variables and recode as factors when needed

#4.1# Variables already in the data

summary(pouting1)

pouting1<-pouting1[-which(is.na(pouting1$Lat)),]
pouting1<-pouting1[-which(is.na(pouting1$Seafloor)),]
pouting1<-pouting1[-which(is.na(pouting1$Depth)),]

str(pouting1)
summary(pouting1)

pouting1$fyear<-factor(pouting1$Year)
pouting1$fgear<-factor(pouting1$Gear)
pouting1$fcrew<-factor(ifelse(pouting1$Crew <= 3,1,2)) # Crew as factor
pouting1$fZoneA<-factor(pouting1$ZoneA)
pouting1$fZoneO<-factor(pouting1$ZoneO)

#4.2# Gear size

#4.2.1# Recorded gear size in each trip

gearSize<-read.csv2(file="C:\\Users\\alex\\Dropbox\\casgass\\dimensiones_artes_pesca.csv",header=T,dec=".",sep=";")

head(gearSize)

gearSize1<-gearSize[gearSize$ARTE=="NASA-PEIXES",] # Only nasa-peixe
gearSize2<-gearSize[gearSize$ARTE=="VETAS",] # Only Vetas
gearSize3<-gearSize[gearSize$ARTE=="MINOS",] # Only Vetas
gearSize1<-gearSize1[,c(1:5)]
gearSize2<-gearSize2[,c(1:5)]
gearSize3<-gearSize3[,c(1:5)]
head(gearSize1)
head(gearSize2)
head(gearSize3)

head(gearSize1)
dim(gearSize1)
head(gearSize2)
dim(gearSize2)
head(gearSize3)
dim(gearSize3)

#nasa
par(mfrow=c(1,3))
summary(gearSize1$longitud);hist(gearSize1$longitud) # Looks OK, units cm (I guess)
summary(gearSize1$ancho);hist(gearSize1$ancho) 
summary(gearSize1$altura);hist(gearSize1$altura) 

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

#join two 
gearSize4<-data.frame(rbind(gearSize1,gearSize2,gearSize3))
head(gearSize4)
tail(gearSize4)
dim(gearSize4)
  
#adding gear size information to main data set
pouting1<-merge(pouting1,gearSize4,by.x=c("Idlance","Gear"),by.y=c("Idlance","ARTE")) # Merge datasets
head(pouting1)

#treatment of gearsize by Gear independently
pouting2<-pouting1[pouting1$Gear=="NASA-PEIXES",]
pouting3<-pouting1[pouting1$Gear=="VETAS",]
pouting4<-pouting1[pouting1$Gear=="MINOS",]

dim(pouting2)
head(pouting2)
dim(pouting3)
head(pouting3)
dim(pouting4)
head(pouting4)

# Fill in NAs with vessel-average data
#nasas
#it does not require any treatment, no NA's
pouting2$longitud2<-pouting2$longitud
pouting2$ancho2<-pouting2$ancho
pouting2$altura2<-pouting2$altura

pouting2$longitud3<-pouting2$longitud
pouting2$ancho3<-pouting2$ancho
pouting2$altura3<-pouting2$altura

pouting2$longitud4<-pouting2$longitud
pouting2$ancho4<-pouting2$ancho
pouting2$altura4<-pouting2$altura

#vetas
ids<-unique(pouting3$Idflota) 
long<-length(ids)
lista<-as.list(1:long)
names(lista)<-ids

for (i in 1:long) {lista[[i]]<-pouting3[pouting3$Idflota==ids[i],]}
for (i in 1:long) {lista[[i]]$longitud2<-ifelse(sum(lista[[i]]$longitud,na.rm=T) > 0,mean(lista[[i]]$longitud,na.rm=T),NA)}
for (i in 1:long) {lista[[i]]$ancho2<-ifelse(sum(lista[[i]]$ancho,na.rm=T) > 0,mean(lista[[i]]$ancho,na.rm=T),NA)}
for (i in 1:long) {lista[[i]]$altura2<-ifelse(sum(lista[[i]]$altura,na.rm=T) > 0,mean(lista[[i]]$altura,na.rm=T),NA)}

pouting3<-do.call(rbind,lista)
head(pouting3)
summary(pouting3)

pouting3$longitud3<-ifelse(is.na(pouting3$longitud),pouting3$longitud2,pouting3$longitud)
pouting3$ancho3<-ifelse(is.na(pouting3$ancho),pouting3$ancho2,pouting3$ancho)
pouting3$altura3<-ifelse(is.na(pouting3$altura),pouting3$altura2,pouting3$altura)
head(pouting3)
tail(pouting3)
summary(pouting3)

pouting3$longitud4<-ifelse(is.na(pouting3$longitud3),
                           mean(pouting3$longitud3,na.rm=T),pouting3$longitud3) # Fill in with fleet average data for those vessels without any information
pouting3$ancho4<-ifelse(is.na(pouting3$ancho3),
                         mean(pouting3$ancho3,na.rm=T),pouting3$ancho3) 
pouting3$altura4<-ifelse(is.na(pouting3$altura3),
                         mean(pouting3$altura3,na.rm=T),pouting3$altura3)
summary(pouting3)

#miños
ids<-unique(pouting4$Idflota) 
long<-length(ids)
lista<-as.list(1:long)
names(lista)<-ids

for (i in 1:long) {lista[[i]]<-pouting4[pouting4$Idflota==ids[i],]}
for (i in 1:long) {lista[[i]]$longitud2<-ifelse(sum(lista[[i]]$longitud,na.rm=T) > 0,mean(lista[[i]]$longitud,na.rm=T),NA)}
for (i in 1:long) {lista[[i]]$ancho2<-ifelse(sum(lista[[i]]$ancho,na.rm=T) > 0,mean(lista[[i]]$ancho,na.rm=T),NA)}
for (i in 1:long) {lista[[i]]$altura2<-ifelse(sum(lista[[i]]$altura,na.rm=T) > 0,mean(lista[[i]]$altura,na.rm=T),NA)}

pouting4<-do.call(rbind,lista)
head(pouting4)
summary(pouting4)

pouting4$longitud3<-ifelse(is.na(pouting4$longitud),pouting4$longitud2,pouting4$longitud)
pouting4$ancho3<-ifelse(is.na(pouting4$ancho),pouting4$ancho2,pouting4$ancho)
pouting4$altura3<-ifelse(is.na(pouting4$altura),pouting4$altura2,pouting4$altura)
head(pouting4)
tail(pouting4)
summary(pouting4)

pouting4$longitud4<-ifelse(is.na(pouting4$longitud3),
                           mean(pouting4$longitud3,na.rm=T),pouting4$longitud3) # Fill in with fleet average data for those vessels without any information
pouting4$ancho4<-ifelse(is.na(pouting4$ancho3),
                        mean(pouting4$ancho3,na.rm=T),pouting4$ancho3) 
pouting4$altura4<-ifelse(is.na(pouting4$altura3),
                         mean(pouting4$altura3,na.rm=T),pouting4$altura3)
summary(pouting4)

pouting5<-data.frame(rbind(pouting2,pouting3,pouting4))
dim(pouting5)
summary(pouting5)

pouting5$Area<-ifelse(is.na(pouting5$ancho4),pouting5$Pieces*pouting5$longitud4*pouting5$altura4,NA) # MIÑ‘OS & VETAS Area (m2) 
summary(pouting5$Area)
hist(pouting5$Area)

pouting5$Vol<-ifelse(is.na(pouting5$ancho4),NA,
                     pouting5$Pieces*(pouting5$longitud4/100)*(pouting5$ancho4/100)*(pouting5$altura4/100)) # volumen NASAS (m3) 
summary(pouting5$Vol)
hist(pouting5$Vol)

#4.2.3# Differences between Nasas, Vetas and Miños

head(pouting5)

gear1<-lm(log(Area)~Gear,data=pouting5)
summary(gear1)
plot(allEffects(gear1))

gear2<-lm(log(Depth)~Gear,data=pouting5)
summary(gear2)
plot(allEffects(gear2))

gear3<-lm(log(Soak)~Gear,data=pouting5)
summary(gear3)
plot(allEffects(gear3))

gear4<-lm(Month~Gear,data=pouting5)
summary(gear4)
plot(allEffects(gear4))

gear5<-lm(GRT~Gear,data=pouting5)
summary(gear5)
plot(allEffects(gear5))

gear6<-lm(Crew~Gear,data=pouting5)
summary(gear6)
plot(allEffects(gear6))

gear7<-lm(ZoneA~Gear,data=pouting5)
summary(gear7)
plot(allEffects(gear7))

# ----------------------------- #
#5# Add environmental data (At this time I run models without AR structure)

#5.1# Load upwelling and remove seasonal cycle

qx<-read.csv2(file="Upwelling.csv",header=T,dec=",",sep=",")

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

sst<-read.table(file="sstGal.selection.txt",header=T,dec=".",sep=" ")

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

pouting5$QxM<-as.vector(rep(NA,dim(pouting5)[1]))

for (i in 1:dim(pouting5)[1]) {pouting5$QxM[i]<-mean(oceano$qxAno[(which(oceano$Year==pouting5$Year[i] &
	oceano$DoY==pouting5$Julian[i])-win.Up):(which(oceano$Year==pouting5$Year[i] &
	oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)}

pouting5$QyM<-as.vector(rep(NA,dim(pouting5)[1]))

for (i in 1:dim(pouting5)[1]) {pouting5$QyM[i]<-mean(oceano$qyAno[(which(oceano$Year==pouting5$Year[i] &
	oceano$DoY==pouting5$Julian[i])-win.Up):(which(oceano$Year==pouting5$Year[i] &
	oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)}

pouting5$sstM<-as.vector(rep(NA,dim(pouting5)[1]))

for (i in 1:dim(pouting5)[1]) {pouting5$sstM[i]<-
	
	if (pouting5$ZoneA[i]=="1") {mean(oceano$sstAnoZ1[(which(oceano$Year==pouting5$Year[i] &
		oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
		oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)} 
	
	else {if (pouting5$ZoneA[i]=="2") {mean(oceano$sstAnoZ2[(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)} 
		
	else {if (pouting5$ZoneA[i]=="3") {mean(oceano$sstAnoZ3[(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)} 
			
	else {if (pouting5$ZoneA[i]=="4") {mean(oceano$sstAnoZ4[(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)} 
				
	else {if (pouting5$ZoneA[i]=="5") {mean(oceano$sstAnoZ5[(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)} 
					
	else {if (pouting5$ZoneA[i]=="6") {mean(oceano$sstAnoZ6[(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)} 
						
	else {if (pouting5$ZoneA[i]=="7") {mean(oceano$sstAnoZ7[(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)} 
							
	else {if (pouting5$ZoneA[i]=="8") {mean(oceano$sstAnoZ8[(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
			oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)} 
								
	else {mean(oceano$sstAnoZ9[(which(oceano$Year==pouting5$Year[i] &
		oceano$DoY==pouting5$Julian[i])-win.Te):(which(oceano$Year==pouting5$Year[i] &
		oceano$DoY==pouting5$Julian[i])-1)],na.rm=T)}
								
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

head(pouting5)
summary(pouting5)

pouting5<-pouting5[-which(is.na(pouting5$QxM)),]

names(pouting5)
pouting6<-pouting5[,-c(30:38)]
dim(pouting6)
names(pouting6)

#6.1# Response variable

#6.1.1# Percentage of zeroes and other stuff

(length(which(pouting6$Ntot==0))*100)/nrow(pouting6) # Percentage of zeros
dim(pouting6)[1] # Number of hauls
length(unique(pouting6$Idflota)) # Number of vessels

pouting6$Gear<-factor(pouting6$Gear)
poll.dist<-as.data.frame(table(pouting6$Ntot,pouting6$Gear))
colnames(poll.dist)<-c("Count","Gear","Freq")
head(poll.dist)

pdf(file="respPouting.pdf",width=15,height=6)

ggplot(data=poll.dist,aes(x=Count,y=Freq))+facet_wrap(~Gear,nrow=3,scales="free_y")+
	geom_bar(stat="identity")+
	scale_y_continuous("Frequency")+
	scale_x_discrete("Number of pouting caught per haul")+
	theme(axis.text.x=element_text(size=7),
		  axis.text.y=element_text(size=8))

dev.off()

par(mar=c(5,5,3,3))
hist(pouting6$Ntot,prob=T,ylab="Probability",xlab="Nº of fished pouting",main="")

hist(log(pouting6$Wtot))

plot(log(Ntot)~log(Wtot),data=pouting6,ylab="Total number",xlab="Total biomass")

#6.1.2# Plots of spatial distribution of effort (m2 of net and fishing hour) 

p1<-ggplot(data=pouting6,aes(x=Lon,y=Lat))
p1+geom_point(aes(color=Area*Soak))+
	scale_color_gradient(low="blue",high="red")
p1+geom_point(aes(color=Vol*Soak))+
  scale_color_gradient(low="blue",high="red")

#6.1.3# Distribution of cases across explanatory variables

table(pouting6$Year,pouting6$Gear) # Few data in 1999
table(pouting6$Year,pouting6$Month,pouting6$Gear)
table(pouting6$fcrew,pouting6$Gear)
table(pouting6$fgear)
table(pouting6$ZoneA,pouting6$Gear)
table(pouting6$ZoneO,pouting6$Gear)

#problema de muestreo con NASAS, años submuestreados
#imposible obtener una serie temporal completa
#de todos modos haremos los modelos para nasas y enmalle

#6.2# Set up the offset

#6.2.1# Exploration of effort variables

head(pouting6)

par(mfrow=c(3,3))
plot(GRT~Crew,data=pouting6)
abline(lm(GRT~Crew,data=pouting6),col="red",lwd=2)
plot(GRT~Area,data=pouting6)
abline(lm(GRT~Area,data=pouting6),col="red",lwd=2)
plot(GRT~Soak,data=pouting6)
abline(lm(GRT~Vol,data=pouting6),col="red",lwd=2)
plot(GRT~Vol,data=pouting6)
abline(lm(GRT~Soak,data=pouting6),col="red",lwd=2)
plot(Crew~Area,data=pouting6)
abline(lm(Crew~Area,data=pouting6),col="red",lwd=2)
plot(Crew~Vol,data=pouting6)
abline(lm(Crew~Vol,data=pouting6),col="red",lwd=2)
plot(Crew~Soak,data=pouting6)
abline(lm(Crew~Soak,data=pouting6),col="red",lwd=2)
plot(Area~Soak,data=pouting6)
abline(lm(Area~Soak,data=pouting6),col="red",lwd=2)
plot(Vol~Soak,data=pouting6)
abline(lm(Vol~Soak,data=pouting6),col="red",lwd=2)

par(mfrow=c(1,3))
boxplot(Area~fcrew,data=pouting6,notch=T,xlab="Crew",ylab="Area")
abline(lm(Area~fcrew,data=pouting6),col="red",lwd=2)
boxplot(Vol~fcrew,data=pouting6,notch=T,xlab="Crew",ylab="Vol")
abline(lm(Vol~fcrew,data=pouting6),col="red",lwd=2)
boxplot(Soak~fcrew,data=pouting6,notch=T,xlab="Crew",ylab="Soak time")
abline(lm(Soak~fcrew,data=pouting6),col="red",lwd=2)

#6.2.2# possible offsets

par(mfrow=c(2,2))
#two for enmalle
pouting6$offs1<-pouting6$Area*(pouting6$Soak/60)
hist(log(pouting6$offs1),xlab="Offset.1",main="",col="grey")

# Offset (num pollack per m2 per h fishing). Considero que al calcular el area
# de la red en funcion del tamaño y numero de paños que tiene cada arte
# las artes se igualan y el esfuerzo en cuanto a diferencia entre artes
# queda ya contemplado, pero no la eficiencia de cada arte. Para ello
# podemos incluir arte como factor. En cuanto a GRT y crew, al tratarse de
# artes pasivas no deberian influir mucho, asi que se pueden incluir como
# explicativas en el modelo (factores o continuas) 

pouting6$offs2<-pouting6$Area
hist(log(pouting6$offs2),xlab="Offset.2",main="")

# Offset secundario. Consideramos solo el tamaÃ±o de la red y usarÃ­amos el tiempo
# de calado como variable explicativa

#three for nasa
pouting6$offs3<-pouting6$Vol*(pouting6$Soak/60)
hist(log(pouting6$offs3),xlab="Offset.3",main="",col="grey") #m3/h

pouting6$offs4<-pouting6$Vol
hist(log(pouting6$offs4),xlab="Offset.4",main="")

pouting6$offs5<-pouting6$Pieces*(pouting6$Soak/60)
hist(log(pouting6$offs5),xlab="Offset.5",main="")

#6.3# Calculation of nighttime fishing

library(lubridate) # Necesaria para obtener horas y minutos de objetos POSIXct

clock.M <- function (t) {hour(t)*60 + minute(t)} # Calcula el minuto entre 1 y 1440 de un dia

pho<-read.csv2(file="Sunrise_set.csv",header=T,dec=".",sep=";")

head(pho)
pho<-pho[,c(1:3,7:12)] # Nos quedamos con las columnas que interesan

#6.3.1# Convertir a POSIXct y obtener informacion sobre las fechas y horas

head(pouting6)

pouting6$Deployment<-as.POSIXct(pouting6$Deployment,format="%Y-%m-%d %H:%M:%S") 
pouting6$Retrieval<-as.POSIXct(pouting6$Retrieval,format="%Y-%m-%d %H:%M:%S")
pho$Rise<-as.POSIXct(pho$Rise,format="%d/%m/%Y %H:%M:%S")
pho$Set<-as.POSIXct(pho$Set,format="%d/%m/%Y %H:%M:%S")

pouting6$yLarg<-year(pouting6$Deployment) # Obtener año de largada
pouting6$mLarg<-month(pouting6$Deployment) # Obtener mes de largada
pouting6$dLarg<-day(pouting6$Deployment) # Obtener dia de largada
pouting6$yVir<-year(pouting6$Retrieval) # Obtener año de virada
pouting6$mVir<-month(pouting6$Retrieval) # Obtener mes de virada
pouting6$dVir<-day(pouting6$Retrieval) # Obtener dia de virada

pouting6$hLarg<-clock.M(pouting6$Deployment) # Calcula minuto largada
pouting6$hVir<-clock.M(pouting6$Retrieval) # Calcula minuto virada
pho$hRise<-clock.M(pho$Rise) # Calcula minuto Sunrise
pho$hSet<-clock.M(pho$Set) # Calcula minuto Sunset

#6.3.2# Obtener informacion cruzando los datos de pesca con fotoperiodo

# Obtener dia unico para la largada de entre todos los dias posibles desde
# el 1/1/1999 (dia 1) hasta el 31/12/2013 (dia 5479) 

pouting6$tLarg<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$tLarg[i]<-pho$Timer[which(pouting6$yLarg[i]==pho$Year & pouting6$mLarg[i]==pho$Month & pouting6$dLarg[i]==pho$Day)]} 

# Obtener dia unico para la virada de entre todos los dias posibles desde
# el 1/1/1999 (dia 1) hasta el 31/12/2013 (dia 5479) 

pouting6$tVir<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$tVir[i]<-pho$Timer[which(pouting6$yVir[i]==pho$Year & pouting6$mVir[i]==pho$Month & pouting6$dVir[i]==pho$Day)]}  

# Obtener minuto Sunrise el dia de largada

pouting6$hRiseL<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$hRiseL[i]<-pho$hRise[which(pouting6$tLarg[i]==pho$Timer)]} 

# Obtener minuto Sunset el dia de largada

pouting6$hSetL<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$hSetL[i]<-pho$hSet[which(pouting6$tLarg[i]==pho$Timer)]} 

# Obtener minuto Sunrise el dia de virada

pouting6$hRiseV<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$hRiseV[i]<-pho$hRise[which(pouting6$tVir[i]==pho$Timer)]} 

# Obtener minuto Sunset el dia de virada

pouting6$hSetV<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$hSetV[i]<-pho$hSet[which(pouting6$tVir[i]==pho$Timer)]} 

# Obtener minutos nocturnos transcurridos entre dias en los casos en que
# pasÃ³ mas de un dia entre largada y virada

pouting6$minNigTot<-as.vector(rep(NA,dim(pouting6)[1]))

for (i in 1:dim(pouting6)[1]) {pouting6$minNigTot[i]<-if (pouting6$tVir[i]-pouting6$tLarg[i]<=1) {0}
	else {sum(pho$Night[which(pho$Timer==pouting6$tLarg[i]+1) : which(pho$Timer== pouting6$tVir[i])])}
	} 

#6.3.3# Calcular minutos nocturnos

#6.3.3.1# Minutos nocturnos si la largada y la virada tienen lugar el mismo dia (6 combinaciones posibles) 

pouting6$minNig1<-as.vector(rep(NA,dim(pouting6)[1]))

for (i in 1:dim(pouting6)[1]) {pouting6$minNig1[i]<-if (pouting6$tVir[i]-pouting6$tLarg[i]==0) {
		
		# Largada y Virada ocurren antes del Sunrise
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] <= pouting6$hRiseL[i],pouting6$hVir[i]-pouting6$hLarg[i],
		
		# Largada ocurre antes del Sunrise y Virada entre Sunrise y Sunset 
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] > pouting6$hRiseL[i] & pouting6$hVir[i] <= pouting6$hSetL[i],pouting6$hRiseL[i]-pouting6$hLarg[i], 
		
		# Largada ocurre antes del Sunrise y Virada despues del Sunset
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] > pouting6$hSetL[i],(pouting6$hRiseL[i]-pouting6$hLarg[i])+(pouting6$hVir[i]-pouting6$hSetL[i]),
		
		# Largada ocurre entre Sunrise y Sunset y Virada ocurre entre Sunrise y Sunset 
		ifelse(pouting6$hLarg[i] >= pouting6$hRiseL[i] & pouting6$hLarg[i] <= pouting6$hSetL & pouting6$hVir[i] >= pouting6$hRiseL & pouting6$hVir[i] <= pouting6$hSetL[i],0,
		
		# Largada ocurre entre Sunrise y Sunset y Virada ocurre despues del Sunset
		ifelse(pouting6$hLarg[i] >= pouting6$hRiseL[i] & pouting6$hLarg[i] <= pouting6$hSetL & pouting6$hVir[i] > pouting6$hSetL[i],pouting6$hVir[i]-pouting6$hSetL[i],
		
		# Largada y Virada ocurren despues del Sunset
		pouting6$hVir[i]-pouting6$hLarg[i])))))
		
	} else {0}

} 

#6.3.3.2# Minutos nocturnos si la virada tiene lugar al dia siguiente de la largada (9 combinaciones posibles)

pouting6$minNig2<-as.vector(rep(NA,dim(pouting6)[1]))

for (i in 1:dim(pouting6)[1]) {pouting6$minNig2[i]<-if (pouting6$tVir[i]-pouting6$tLarg[i]==1) {
		
		# Largada ocurre antes del Sunrise y Virada ocurre antes del Sunrise del dia siguiente 
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] <= pouting6$hRiseV[i],(pouting6$hRiseL[i]-pouting6$hLarg[i])+((1440-pouting6$hSetL[i])+pouting6$hVir[i]),
		
		# Largada ocurre antes del Sunrise y Virada ocurre entre el Sunrise y el Sunset del dia siguiente
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] > pouting6$hRiseV[i] & pouting6$hVir[i] <= pouting6$hSetV[i],(pouting6$hRiseL[i]-pouting6$hLarg[i])+((1440-pouting6$hSetL[i])+pouting6$hRiseV[i]),
		
		# Largada ocurre antes del Sunrise y Virada ocurre despues del Sunset del dia siguiente
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] > pouting6$hSetV[i],(pouting6$hRiseL[i]-pouting6$hLarg[i])+(1440-pouting6$hSetL[i])+pouting6$hRiseV[i]+(pouting6$hVir[i]-pouting6$hSetV[i]),
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre antes del Sunrise del dia siguiente 
		ifelse(pouting6$hLarg[i] >= pouting6$hRiseL[i] & pouting6$hLarg[i] <= pouting6$hSetL[i] & pouting6$hVir[i] <= pouting6$hRiseV[i],(1440-pouting6$hSetL[i])+pouting6$hVir[i],
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre entre el Sunrise y el Sunset del dia siguiente 
		ifelse(pouting6$hLarg[i] >= pouting6$hRiseL[i] & pouting6$hLarg[i] <= pouting6$hSetL[i] & pouting6$hVir[i] >= pouting6$hRiseV[i] & pouting6$hVir[i] <= pouting6$hSetV[i],(1440-pouting6$hSetL[i])+pouting6$hRiseV[i],
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre despues del Sunset del dia siguiente
		ifelse(pouting6$hLarg[i] >= pouting6$hRiseL[i] & pouting6$hLarg[i] <= pouting6$hSetL[i] & pouting6$hVir[i] > pouting6$hSetV[i],(1440-pouting6$hSetL[i])+pouting6$hRiseV[i]+(pouting6$hVir[i]-pouting6$hSetV[i]),
		
		# Largada ocurre despues del Sunset y virada ocurre antes del Sunrise del dia siguiente
		ifelse(pouting6$hLarg[i] > pouting6$hSetL[i] & pouting6$hVir[i] <= pouting6$hRiseV[i],(1440-pouting6$hLarg[i])+pouting6$hVir[i],
		
		# Largada ocurre despues del Sunset y virada ocurre entre el Sunrise y el Sunset del dia siguiente
		ifelse(pouting6$hLarg[i] > pouting6$hSetL[i] & pouting6$hVir[i] >= pouting6$hRiseV[i] & pouting6$hVir[i] <= pouting6$hSetV[i],(1440-pouting6$hLarg[i])+pouting6$hRiseV[i],
		
		# Largada ocurre despues del Sunset y virada ocurre despues del Sunset del dia siguiente
		(1440-pouting6$hLarg[i])+pouting6$hRiseV[i]+(pouting6$hVir[i]-pouting6$hSetV[i])))))))))
		
	} else {0}

}

#6.3.3.3# Minutos nocturnos si entre la largada y la virada pasa mas de un dia (las mismas 9 combinaciones posibles)

pouting6$minNig3<-as.vector(rep(NA,dim(pouting6)[1]))

for (i in 1:dim(pouting6)[1]) {pouting6$minNig3[i]<-if (pouting6$tVir[i]-pouting6$tLarg[i] > 1) {
		
		# Largada ocurre antes del Sunrise y Virada ocurre antes del Sunrise varios dias despues
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] <= pouting6$hRiseV[i],(pouting6$hRiseL[i]-pouting6$hLarg[i])+(pouting6$minNigTot[i]-(pouting6$hRiseV[i]-pouting6$hVir[i])),
		
		# Largada ocurre antes del Sunrise y Virada ocurre entre el Sunrise y el Sunset varios dias despues
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] > pouting6$hRiseV[i] & pouting6$hVir[i] <= pouting6$hSetV[i],(pouting6$hRiseL[i]-pouting6$hLarg[i])+pouting6$minNigTot[i],
		
		# Largada ocurre antes del Sunrise y Virada ocurre despues del Sunset varios dias despues
		ifelse(pouting6$hLarg[i] < pouting6$hRiseL[i] & pouting6$hVir[i] > pouting6$hSetV[i],(pouting6$hRiseL[i]-pouting6$hLarg[i])+pouting6$minNigTot[i]+(pouting6$hVir[i]-pouting6$hSetV[i]),
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre antes del Sunrise varios dias despues
		ifelse(pouting6$hLarg[i] >= pouting6$hRiseL[i] & pouting6$hLarg[i] <= pouting6$hSetL[i] & pouting6$hVir[i] <= pouting6$hRiseV[i],pouting6$minNigTot[i]-(pouting6$hRiseV[i]-pouting6$hVir[i]),
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre entre el Sunrise y el Sunset varios dias despues
		ifelse(pouting6$hLarg[i] >= pouting6$hRiseL[i] & pouting6$hLarg[i] <= pouting6$hSetL[i] & pouting6$hVir[i] >= pouting6$hRiseV[i] & pouting6$hVir[i] <= pouting6$hSetV[i],pouting6$minNigTot[i],
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre despues del Sunset varios dias despues
		ifelse(pouting6$hLarg[i] >= pouting6$hRiseL[i] & pouting6$hLarg[i] <= pouting6$hSetL[i] & pouting6$hVir[i] > pouting6$hSetV[i],pouting6$minNigTot[i]+(pouting6$hVir[i]-pouting6$hSetV[i]),
		
		# Largada ocurre despues del Sunset y virada ocurre antes del Sunrise varios dias despues
		ifelse(pouting6$hLarg[i] > pouting6$hSetL[i] & pouting6$hVir[i] <= pouting6$hRiseV[i],pouting6$minNigTot[i]-((pouting6$hLarg[i]-pouting6$hSetL[i])-(pouting6$hRiseV[i]-pouting6$hVir[i])),
		
		# Largada ocurre despues del Sunset y virada ocurre entre el Sunrise y el Sunset varios dias despues
		ifelse(pouting6$hLarg[i] > pouting6$hSetL[i] & pouting6$hVir[i] >= pouting6$hRiseV[i] & pouting6$hVir[i] <= pouting6$hSetV[i],pouting6$minNigTot[i]-(pouting6$hLarg[i]-pouting6$hSetL[i]),
		
		# Largada ocurre despues del Sunset y virada ocurre despues del Sunset varios dias despues
		(pouting6$minNigTot[i]-(pouting6$hLarg[i]-pouting6$hSetL[i]))+(pouting6$hVir[i]-pouting6$hSetV[i])))))))))
	
	} else {0}

}

#6.3.4# Nueva variable '% minutos nocturnos'

#6.3.4.1# Check that soak from UTPB is correct

# pouting6$Soaktime<-as.numeric(difftime(pouting6$Retrieval,pouting6$Deployment,units="mins"))

# plot(pouting6$Soak~Soaktime)
# soakTest<-ifelse(pouting6$Soak==Soaktime,0,1)
# sum(soakTest) # Los dos valores son iguales asi que usamos el suyo
# nn<-which(soakTest==1,)
# pouting6[nn,]

#6.3.4.2# Calculate proportion of night time

pouting6$caladoNight<-as.vector(rep(NA,dim(pouting6)[1])) 

for (i in 1:dim(pouting6)[1]) {pouting6$caladoNight[i]<-if (pouting6$tVir[i]==pouting6$tLarg[i]) {
	
	# % minutos nocturnos si largada y virada ocurren el mismo dia
	round(pouting6$minNig1[i]/pouting6$Soak[i],2)
	
	} else {if (pouting6$tVir[i]-pouting6$tLarg[i]==1) {
		
		# % minutos nocturnos si virada ocurre al dia siguiente de largada
		round(pouting6$minNig2[i]/pouting6$Soak[i],2)
		
		} else {
			
			# % minutos nocturnos si virada ocurre varios dias despues de largada
			round(pouting6$minNig3[i]/pouting6$Soak[i],2)
		}
	}
}

summary(pouting6$caladoNight)
hist(pouting6$caladoNight)

#6.3.5# Standarizar hora de largada segun fotoperiodo

#6.3.5.1# Obtener minuto Sunrise el dia siguiente de largada

names(pouting6)

pouting6$hRiseL1<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$hRiseL1[i]<-pho$hRise[which(pouting6$tLarg[i]+1==pho$Timer)]} 

#6.3.5.2# Corregir hora de largada segun fotoperiodo: si la hora de largada ocurre antes del
# sunset restamos hora de largada-hora sunrise (con lo que largadas antes de sunrise
# tendran valores negativos, y largadas despues de sunrise hasta sunset tedran valores
# positivos. Si la largada ocurre despues del sunset los valores seran muy negativos
# por que los contaremos respecto del sunrise del dia siguiente)

pouting6$hLargC<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$hLargC[i]<-if (pouting6$hLarg[i] <= pouting6$hSetL[i]) {
		
		pouting6$hLarg[i] - pouting6$hRiseL[i]
		
		} else {
			
			((1440 - pouting6$hLarg[i]) + pouting6$hRiseL1[i])*-1
			
		}
	}

#6.3.6# Standarizar hora de virada segun fotoperiodo

#6.3.6.1# Obtener minuto Sunrise el dia siguiente de largada

names(pouting6)

pouting6$hSetV1<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$hSetV1[i]<-pho$hSet[which(pouting6$tVir[i]-1==pho$Timer)]} 

#6.3.6.2# Corregir hora de virada segun fotoperiodo: si la hora de virada ocurre antes del
# sunrise restamos hora de virada-hora sunset del dia anterior (con lo que viradas mas
# lejos del sunset tienen valores muy negativos). Si la virada ocurre despues del sunrise
# los valores seran positivos

pouting6$hVirC<-as.vector(rep(NA,dim(pouting6)[1]))
for (i in 1:dim(pouting6)[1]) {pouting6$hVirC[i]<-if (pouting6$hVir[i] <= pouting6$hRiseV[i]) {
		
		((1440 - pouting6$hSetV1[i]) + pouting6$hVir[i])*-1
		
		} else {
			
			pouting6$hVir[i] - pouting6$hRiseV[i]
			
		}
	}

#6.3.7# Uso de factores (vale 0 si la proporcion de calado nocturno es
# mayor que 75%, si esta entre el 25 y 75% vale 1, y si es menor o igual que
# el 25%, esto es, el lance es practicamente todo diurno, vale 2)

pouting6$Period<-factor(ifelse(pouting6$caladoNight <= 0.25,2,
	ifelse(pouting6$caladoNight > 0.25 & pouting6$caladoNight <= 0.75,1,0)))

#6.4 added by Alex Alonso
#añadimos una nueva variable indicativa del tamaño promedio de los individuos capturados
#Wtot/Ntot puede ser utilizada como indicativo de cambios de distribución ontogénicos???
summary(pouting6$Ntot)
summary(pouting6$Wtot)
plot(pouting6$Ntot,pouting6$Wtot)
#necesitamos eliminar NA de los datos Wtot para poder calcular esta nueva variable
pouting6$Wtot<-ifelse(pouting6$Ntot==0,0,pouting6$Wtot)
par(mfrow=c(1,2))
dotchart(pouting6$Wtot)
hist(pouting6$Wtot,breaks=100)
summary(pouting6$Wtot)
plot(pouting6$Ntot,pouting6$Wtot)
#calculamos la relación
pouting6$AvgSize<-ifelse(pouting6$Ntot==0,0,pouting6$Wtot/pouting6$Ntot)
hist(pouting6$AvgSize,breaks=100)
summary(pouting6$AvgSize)
#does it make sense 0 values for this variable?
plot(pouting6$Ntot,pouting6$AvgSize)


# ----------------------------- #
#7# Exploratory steps before modelling

#7.1# Distribution of all potential explanatory variables

#likely outlier in enmalle >500
which(pouting6$Depth>500) # Solo un lance por encima de 200 m (utpb dice que es real; pero no representativo de la pesquería)
pouting7<-pouting6[pouting6$Depth<500,]

d1<-ggplot(data=pouting7,aes(x=Crew))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Crew")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))
d2<-ggplot(data=pouting7,aes(x=GRT))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("GRT")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d3<-ggplot(data=pouting7,aes(x=Area))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous(expression(paste("Gear area"," ","(",m^2,")")))+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d4<-ggplot(data=pouting7,aes(x=Soak/1440))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Soak time (days)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d5<-ggplot(data=pouting7,aes(x=hLarg/60))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Deployment (Local time)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d6<-ggplot(data=pouting7,aes(x=hVir/60))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Retrieval (Local time)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d7<-ggplot(data=pouting7,aes(x=Julian))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Day of the Year")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d8<-ggplot(data=pouting7,aes(x=Lat))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Latitude (ºN)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d9<-ggplot(data=pouting7,aes(x=Lon))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Longitude (ºW)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d10<-ggplot(data=pouting7,aes(x=Depth))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Depth (m)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d11<-ggplot(data=pouting7,aes(x=QxM))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous(expression(paste("Qx"," ","(",m^3,s^-1,km^-1,") Ã—",10^3)))+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d12<-ggplot(data=pouting7,aes(x=sstM))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("SST (ºC)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")
d13<-ggplot(data=pouting7,aes(x=AvgSize))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Individual size (kg)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  theme(legend.position="none")

pdf(file="explPouting.pdf",width=20,height=8)

multiplot(d1,d7,d2,d8,d3,d9,d4,d10,d5,d11,d6,d12,cols=6)

dev.off()

ggplot(data=pouting7,aes(x=caladoNight))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Night-time soak (%)")+
  scale_y_continuous("Density")+scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))

ggplot(data=pouting7,aes(x=hLargC))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Corrected deployment time")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))

ggplot(data=pouting7,aes(x=hVirC))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Corrected retrieval time")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("#3182bd","orange","light blue"),labels=c("Miño","Nasa","Veta"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))

summary(pouting7[,c("GRT","Area","hLarg","hVir","Depth","Soak")])

#7.2 Simple relationships

names(pouting7)
write.table(pouting7,"pouting_dat.txt",col.names = TRUE,sep=",")

pouting.dat<-read.table("pouting_dat.txt",header=T,dec=".",sep=",")

names(pouting.dat)
dim(pouting.dat)
summary(pouting.dat)

names(pouting.dat)

#data base for nasa
nasa.pouting<-pouting.dat[pouting.dat$Gear=="NASA-PEIXES",]
dim(nasa.pouting)
names(nasa.pouting)
nasa.pouting<-nasa.pouting[,-c(33,38,39,43:48,51:60,62,64)]

#data base for enmalle
enmalle.pouting<-pouting.dat[pouting.dat$Gear=="MINOS"|pouting.dat$Gear=="VETAS",]
dim(enmalle.pouting)
names(enmalle.pouting)
enmalle.pouting<-enmalle.pouting[,-c(31,34,40:48,51:60,62,64)]

par(mfrow=c(4,4))
plot(log(Ntot+1)~log(GRT),data=nasa.pouting)
boxplot(log(Ntot+1)~fcrew,data=nasa.pouting)
boxplot(log(Ntot+1)~fyear,data=nasa.pouting)
plot(log(Ntot+1)~Julian,data=nasa.pouting)
plot(log(Ntot+1)~Depth,data=nasa.pouting)
plot(log(Ntot+1)~QxM,data=nasa.pouting)
plot(log(Ntot+1)~QyM,data=nasa.pouting)
plot(log(Ntot+1)~sstM,data=nasa.pouting)
boxplot(log(Ntot+1)~ZoneA,data=nasa.pouting)
boxplot(log(Ntot+1)~ZoneO,data=nasa.pouting)
plot(log(Ntot+1)~Soak,data=nasa.pouting)
hist(nasa.pouting$GRT)
hist(nasa.pouting$Depth)
hist(nasa.pouting$QxM)
hist(nasa.pouting$sstM)

par(mfrow=c(4,4))
plot(log(Ntot+1)~log(GRT),data=enmalle.pouting)
boxplot(log(Ntot+1)~fcrew,data=enmalle.pouting)
boxplot(log(Ntot+1)~fyear,data=enmalle.pouting)
plot(log(Ntot+1)~Julian,data=enmalle.pouting)
plot(log(Ntot+1)~Depth,data=enmalle.pouting)
plot(log(Ntot+1)~QxM,data=enmalle.pouting)
plot(log(Ntot+1)~QyM,data=enmalle.pouting)
plot(log(Ntot+1)~sstM,data=enmalle.pouting)
boxplot(log(Ntot+1)~ZoneA,data=enmalle.pouting)
boxplot(log(Ntot+1)~ZoneO,data=enmalle.pouting)
plot(log(Ntot+1)~Soak,data=enmalle.pouting)
hist(enmalle.pouting$GRT,breaks=50)
hist(enmalle.pouting$Depth,breaks=50)
hist(enmalle.pouting$QxM,breaks=50)
hist(enmalle.pouting$sstM,breaks=50)

#7.3# Collinearity

head(nasa.pouting)
dim(nasa.pouting)
head(enmalle.pouting)
dim(enmalle.pouting)

#transform depth and GRT for enmalles
nasa.pouting$lGRT<-log(nasa.pouting$GRT)
nasa.pouting$lDepth<-log(nasa.pouting$Depth)

enmalle.pouting$lGRT<-log(enmalle.pouting$GRT)
enmalle.pouting$lDepth<-log(enmalle.pouting$Depth)

#VIF
vifs1<-c("lGRT","Crew","Year","Julian","Depth","QxM","sstM","caladoNight")
vifs2<-c("lGRT","Crew","Year","Julian","lDepth","QxM","sstM","caladoNight")
corvif(nasa.pouting[,vifs1])
corvif(enmalle.pouting[,vifs2])


# ----------------------------- #
#8# Modelling of catches using non-mixed models

#factores
names(nasa.pouting)
nasa.pouting$fyear<-as.factor(nasa.pouting$fyear)
nasa.pouting$fcrew<-as.factor(nasa.pouting$fcrew)
nasa.pouting$fZoneA<-as.factor(nasa.pouting$fZoneA)
nasa.pouting$fZoneO<-as.factor(nasa.pouting$fZoneO)

names(enmalle.pouting)
enmalle.pouting$fyear<-as.factor(enmalle.pouting$fyear)
enmalle.pouting$fcrew<-as.factor(enmalle.pouting$fcrew)
enmalle.pouting$fZoneA<-as.factor(enmalle.pouting$fZoneA)
enmalle.pouting$fZoneO<-as.factor(enmalle.pouting$fZoneO)

#8.1#NASA modelling
#para simplificar el cáculo de offset asumimos uniformidad en las nasas y utilizamo el offset5<-nasas*hora

#8.1.1# GAM Poisson modelling

#many problems in sample distribution in nasa fanequeira
table(nasa.pouting$Month,nasa.pouting$Year) #this is specially critical in 2012
table(nasa.pouting$Year,nasa.pouting$fZoneO)
table(nasa.pouting$Year,nasa.pouting$Seafloor)
table(nasa.pouting$fZoneO,nasa.pouting$Seafloor)
#this will limit a lot conclusions, so, simplify as much as posible the model to extract robust and general conclusions

#lGrt pocos datos con extremos muy influyentes => coufounding factor con boat
table(nasa.pouting$Idflota,nasa.pouting$lGRT) #we will force this variable as linear!!!
#caladoNight also as linear!!!

nasa.poiss0<-gam(Ntot~offset(log(offs5))+
                   s(lGRT,k=3)+
                   fyear+s(Julian,k=6,bs="cc")+
                   s(Depth,k=3)+
                   s(sstM,k=3)+
                   s(caladoNight,k=3),
                 family=poisson,data=nasa.pouting)

summary(nasa.poiss0)
anova(nasa.poiss0)
plot(nasa.poiss0,pages=1,all.terms=T,scale=0,shade="T")
sum(residuals(nasa.poiss0,type="pearson")^2)/nasa.poiss0$df.resid # Overdispersion

#8.1.2# GLM Poisson modelling

nasa.poiss1<-glm(Ntot~offset(log(offs5))+
				 lGRT+
			 	 fyear+poly(Julian,2)+
			 	 Depth+
			 	 sstM+
			 	 caladoNight,
			 	 family=poisson,data=nasa.pouting)

summary(nasa.poiss1)
Anova(nasa.poiss1)
plot(allEffects(nasa.poiss1))
sum(residuals(nasa.poiss1,type="pearson")^2)/nasa.poiss1$df.resid # Overdispersion

#8.1.3# Negative Binomial modelling

nasa.negbin1<-glm.nb(Ntot~offset(log(offs5))+
                       lGRT+
                       fyear+poly(Julian,2)+
                       Depth+
                       sstM+
                       caladoNight,
                     data=nasa.pouting)

summary(nasa.negbin1)
Anova(nasa.negbin1)
plot(allEffects(nasa.negbin1))
par(mfrow=c(1,3))
hist(residuals(nasa.negbin1,type="deviance"))
plot(residuals(nasa.negbin1,type="deviance")~log(fitted(nasa.negbin1)))
boxplot(residuals(nasa.negbin1,type="deviance")~nasa.pouting$Idflota)
abline(0,0,col="red")
sum(residuals(nasa.negbin1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin1))+1))

#8.1.3.1. Mixed model NB

#http://glmmadmb.r-forge.r-project.org/glmmADMB.html
#https://rpubs.com/bbolker/glmmchapter

nasa.pouting$Idflota<-as.factor(nasa.pouting$Idflota)
nasa.pouting$Nasaoffset<-log(nasa.pouting$offs5)

nasa.negbinran1 <- glmmadmb(Ntot~lGRT+
                              fyear+poly(Julian,2)+
                              Depth+
                              sstM+
                              caladoNight+
                              offset(Nasaoffset)+
                              (1|Idflota), 
                            data=nasa.pouting,
                         zeroInflation=FALSE, 
                         family="nbinom")
summary(nasa.negbinran1)
Anova(nasa.negbinran1)
sum(residuals(nasa.negbinran1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbinran1))+1)) 
par(mfrow=c(1,3))
hist(residuals(nasa.negbinran1,type="pearson"))
plot(residuals(nasa.negbinran1,type="pearson")~log(fitted(nasa.negbinran1)))
boxplot(residuals(nasa.negbinran1,type="pearson")~nasa.pouting$Idflota)
abline(0,0,col="red")

#8.1.4# comparing model approaches

AICtab(nasa.poiss0,nasa.poiss1,nasa.negbinran1,nasa.negbin1)
BICtab(nasa.poiss0,nasa.poiss1,nasa.negbinran1,nasa.negbin1)

#mejoramos el modelo con los random effects, what shall we do?
#el efecto barco está recogiendo algo de variabilidad

#8.2#ENMALLE modelling

#8.2.1# GAM Poisson modelling

enmalle.poiss0<-gam(Ntot~offset(log(offs1))+
                      s(lGRT,by=Gear,k=3)+
                      fyear*fZoneO+
                      s(Julian,k=6,bs="cc")+
                      s(lDepth,by=Gear,k=3)+
                      s(sstM,k=3)+
                      s(caladoNight,by=Gear,k=3)+
                      Seafloor,
                    family=poisson,data=enmalle.pouting)

summary(enmalle.poiss0)
anova(enmalle.poiss0)
plot(enmalle.poiss0,pages=1,all.terms=T,scale=0,shade="T")
sum(residuals(enmalle.poiss0,type="pearson")^2)/enmalle.poiss0$df.resid # Overdispersion

#8.2.2# GLM Poisson modelling

enmalle.poiss1<-glm(Ntot~offset(log(offs1))+
                      Gear*lGRT+
                      fyear*fZoneO+
                      poly(Julian,2)+
                      Gear*lDepth+
                      sstM+
                      Gear*caladoNight+
                      Seafloor,
                 family=poisson,data=enmalle.pouting)

summary(enmalle.poiss1)
Anova(enmalle.poiss1)
plot(allEffects(enmalle.poiss1))
sum(residuals(enmalle.poiss1,type="pearson")^2)/enmalle.poiss1$df.resid # Overdispersion

#8.2.3# Negative Binomial modelling

enmalle.negbin1<-glm.nb(Ntot~offset(log(offs1))+
                          Gear*lGRT+
                          fyear*fZoneO+
                          poly(Julian,2)+
                          Gear*lDepth+
                          sstM+
                          Gear*caladoNight+
                          Seafloor,
                        data=enmalle.pouting)

summary(enmalle.negbin1)
Anova(enmalle.negbin1)
plot(allEffects(enmalle.negbin1))
par(mfrow=c(1,3))
hist(residuals(enmalle.negbin1,type="deviance"))
plot(residuals(enmalle.negbin1,type="deviance")~log(fitted(enmalle.negbin1)))
boxplot(residuals(enmalle.negbin1,type="deviance")~enmalle.pouting$Idflota)
abline(0,0,col="red")
sum(residuals(enmalle.negbin1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbin1))+1)) 

#8.2.3.1. Mixed model NB

#http://glmmadmb.r-forge.r-project.org/glmmADMB.html

enmalle.pouting$Idflota<-as.factor(enmalle.pouting$Idflota)
enmalle.pouting$Enmalleoffset<-log(enmalle.pouting$offs1)

enmalle.negbinran1 <- glmmadmb(Ntot~Gear*lGRT+
                                 fyear*fZoneO+
                                 poly(Julian,2)+
                                 Gear*lDepth+
                                 sstM+
                                 Gear*caladoNight+
                                 Seafloor+
                                 offset(Enmalleoffset)+
                                 (1|Idflota),
                               data=enmalle.pouting,
                               zeroInflation=FALSE,
                               family="nbinom")
summary(enmalle.negbinran1)
Anova(enmalle.negbinran1)
par(mfrow=c(1,3))
hist(residuals(enmalle.negbinran1,type="pearson"))
plot(residuals(enmalle.negbinran1,type="pearson")~log(fitted(enmalle.negbinran1)))
boxplot(residuals(enmalle.negbinran1,type="pearson")~enmalle.pouting$Idflota)
abline(0,0,col="red")
sum(residuals(enmalle.negbinran1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbinran1))+1)) 

#8.2.4# comparing model approaches

AICtab(enmalle.poiss0,enmalle.poiss1,enmalle.negbinran1,enmalle.negbin1)
BICtab(enmalle.poiss0,enmalle.poiss1,enmalle.negbinran1,enmalle.negbin1)

#mejoramos el modelo con los random effects, what shall we do?

#según los parámetros de sobredispersión el mixed negative binomial es el más adecuado
# parece que no será necesario recurrir a los zero inflated en ninguno de los casos

#8.3# en cualquier caso probamos con un modelo zero inflated simple para ambos tipos de arte
#it does not make sense to try ZI with nasa

#8.3.1# exploring binomial 
enmalle.pouting$binNtot<-ifelse(enmalle.pouting$Ntot==0,0,1)

enmalle.bin<-gam(binNtot~Gear+s(lGRT,k=3)+
                fyear+fZoneO+
                s(Julian,k=6,bs="cc")+
                s(lDepth,k=3)+
                s(sstM,k=3)+
                s(caladoNight,k=3)+
                Seafloor,
              family=binomial,data=enmalle.pouting)
summary(enmalle.bin)
anova(enmalle.bin)
plot(enmalle.bin,pages=1,all.terms=T,scale=0,shade="T")

#depth and calado night result in positive effects on pouting catch

#8.3.2# enmalle
enmalle.zipoiss1<-zeroinfl(Ntot~offset(log(offs1))+
                             Gear*lGRT+
                             fyear*fZoneO+
                             poly(Julian,2)+
                             Gear*lDepth+
                             sstM+
                             Gear*caladoNight+
                             Seafloor|
                             lDepth,
                        dist="poisson",link="logit",data=enmalle.pouting)

summary(enmalle.zipoiss1)
sum(residuals(enmalle.zipoiss1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zipoiss1))))

enmalle.zinb1<-zeroinfl(Ntot~offset(log(offs1))+
                          Gear*lGRT+
                          fyear*fZoneO+
                          poly(Julian,2)+
                          Gear*lDepth+
                          sstM+
                          Gear*caladoNight+
                          Seafloor|
                          lDepth,
                     dist="negbin",link="logit",data=enmalle.pouting)

summary(enmalle.zinb1)
sum(residuals(enmalle.zinb1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zinb1))+1))

#glmmadmb does not allow to include explanmatory variables in ZI part
enmalle.zinbran1 <- glmmadmb(Ntot~Gear*lGRT+
                               fyear*fZoneO+
                               poly(Julian,2)+
                               Gear*lDepth+
                               sstM+
                               Gear*caladoNight+
                               Seafloor+
                               offset(Enmalleoffset)+
                            (1|Idflota), 
                          data=enmalle.pouting,
                          zeroInflation=TRUE, 
                          family="nbinom")

#it does not work properly
# :(

###############################################
# Hurdle modelling (assumption: two step model)
#we do not considered to include this in the model comparison
#but we keep it in the code just to see other posibilities

# nasa
#nasa.hurdle1<-hurdle(Ntot~offset(log(offs5))+fcrew+GRT+fyear+ns(Julian,2)+Depth+ns(Julian,2):Depth+
#                       QxM+sstM+caladoNight+fZoneO+Seafloor|Depth,
#                     dist="poisson",link="logit",data=nasa.pouting)
#summary(nasa.hurdle1)

#nasa.hurdle2<-hurdle(Ntot~offset(log(offs5))+fcrew+GRT+fyear+ns(Julian,2)+Depth+ns(Julian,2):Depth+
#                       QxM+sstM+caladoNight+fZoneO+Seafloor|Depth,
#                     dist="negbin",link="logit",data=nasa.pouting)
#summary(nasa.hurdle2)

# enmalle
#enmalle.hurdle1<-hurdle(Ntot~offset(log(offs1))+fcrew+GRT+fyear+ns(Julian,2)+Depth+ns(Julian,2):Depth+
#                       QxM+sstM+caladoNight+fZoneO+Seafloor|Depth,
#                     dist="poisson",link="logit",data=enmalle.pouting)
#summary(enmalle.hurdle1)

#enmalle.hurdle2<-hurdle(Ntot~offset(log(offs1))+fcrew+GRT+fyear+ns(Julian,2)+Depth+ns(Julian,2):Depth+
#                       QxM+sstM+caladoNight+fZoneO+Seafloor|Depth,
#                     dist="negbin",link="logit",data=enmalle.pouting)
#summary(enmalle.hurdle2)
###############################################

#8.4# full model comparison
# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

#8.4.1# AIC/BIC

#AIC
AICtab(nasa.poiss0,nasa.poiss1,nasa.negbinran1,nasa.negbin1)
AICtab(enmalle.poiss0,enmalle.poiss1,enmalle.negbinran1,enmalle.negbin1,enmalle.zipoiss1,enmalle.zinb1)

#BIC
nasa.bic=AIC(nasa.poiss0,nasa.poiss1,nasa.negbinran1,nasa.negbin1,k=log(dim(nasa.pouting)[1]))
nasa.bic[order(nasa.bic$AIC), ]
enmalle.bic=AIC(enmalle.poiss0,enmalle.poiss1,enmalle.negbinran1,enmalle.negbin1,enmalle.zipoiss1,enmalle.zinb1,k=log(dim(enmalle.pouting)[1]))
enmalle.bic[order(enmalle.bic$AIC), ]

#8.4.2# overdispersion parameter

#nasa
nasa.GAMpoisson=sum(residuals(nasa.poiss0,type="pearson")^2)/nasa.poiss0$df.resid
nasa.GLMpoisson=sum(residuals(nasa.poiss1,type="pearson")^2)/nasa.poiss1$df.resid
nasa.NB=sum(residuals(nasa.negbin1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin1))+1)) 
nasa.MixedNB=sum(residuals(nasa.negbinran1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbinran1))+1)) 

nasa.model=c("nasa.GAMpoisson","nasa.GLMpoisson","nasa.NB","nasa.MixedNB")
nasa.theta=c(nasa.GAMpoisson,nasa.GLMpoisson,nasa.NB,nasa.MixedNB)
nasa.M=data.frame(cbind(nasa.model,nasa.theta))
nasa.M=nasa.M[order(nasa.M$nasa.theta),]
nasa.M

#enmalles
enmalle.GAMpoisson=sum(residuals(enmalle.poiss0,type="pearson")^2)/enmalle.poiss0$df.resid
enmalle.GLMpoisson=sum(residuals(enmalle.poiss1,type="pearson")^2)/enmalle.poiss1$df.resid
enmalle.NB=sum(residuals(enmalle.negbin1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbin1))+1)) 
enmalle.MixedNB=sum(residuals(enmalle.negbinran1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbinran1))+1)) 
enmalle.ZIpoisson=sum(residuals(enmalle.zipoiss1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zipoiss1))))
enmalle.ZINB=sum(residuals(enmalle.zinb1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zinb1))+1))

enmalle.model=c("enmalle.GAMpoisson","enmalle.GLMpoisson","enmalle.NB","enmalle.MixedNB","enmalle.ZIpoisson","enmalle.ZINB")
enmalle.theta=c(enmalle.GAMpoisson,enmalle.GLMpoisson,enmalle.NB,enmalle.MixedNB,enmalle.ZIpoisson,enmalle.ZINB)
enmalle.M=data.frame(cbind(enmalle.model,enmalle.theta))
enmalle.M=enmalle.M[order(enmalle.M$enmalle.theta),]
enmalle.M

#If the dispersion parameter is significantly greater than one, indicating overdispersion (variance greater than the mean),
#then the scale parameter should be used to adjust the variance.  
#Failing to account for the overdispersion can result in inflated test statistics. 
#However, when the dispersion parameter is less than one, then the test statistics become more conservative,
#which is not considered as much of a problem.

#8.4.3# loglikelihood

fm.nasa1<-list("ML-Pois"=nasa.poiss1,"NB"=nasa.negbin1,"NB-mix"=nasa.negbinran1)
fm.enmalle1<-list("ML-Pois"=enmalle.poiss1,"NB"=enmalle.negbin1,"NB-mix"=enmalle.negbinran1,"ZIP"=enmalle.zipoiss1,"ZINB"=enmalle.zinb1)

logliks.nasa1<-rbind(logLik=sapply(fm.nasa1,function (x) round(logLik(x),digits=0)),
               Df=sapply(fm.nasa1,function (x) attr(logLik(x),"df")))
logliks.nasa1
logliks.enmalle1<-rbind(logLik=sapply(fm.enmalle1,function (x) round(logLik(x),digits=0)),
                     Df=sapply(fm.enmalle1,function (x) attr(logLik(x),"df")))
logliks.enmalle1

#8.4.4# Vuong test to compare NB vs ZINB
# do not run with random models

vuong(enmalle.zinb1,enmalle.negbin1)

#8.4.5# Predicting probabilities

#nasa
phat.pois.nasa<-predprob(nasa.poiss1)
phat.pois.mn.nasa<-apply(phat.pois.nasa,2,mean)
phat.nb.nasa<-predprob(nasa.negbin1)
phat.nb.mn.nasa<-apply(phat.nb.nasa,2,mean) 

with(nasa.pouting,{
  hist(Ntot,prob=TRUE,col="grey",breaks=seq(-0.5,246,1),xlab="",main="")
  lines(x=seq(0,245,1),y=phat.pois.mn.nasa,type="l",lwd=2,col="red")
  lines(x=seq(0,245,1),y=phat.nb.mn.nasa,type="l",lwd=2,col="blue")
})

#nasa
phat.pois.enmalle<-predprob(enmalle.poiss1)
phat.pois.mn.enmalle<-apply(phat.pois.enmalle,2,mean)
phat.nb.enmalle<-predprob(enmalle.negbin1)
phat.nb.mn.enmalle<-apply(phat.nb.enmalle,2,mean)
phat.zipoiss.enmalle<-predprob(enmalle.zipoiss1)
phat.zipoiss.mn.enmalle<-apply(phat.zipoiss.enmalle,2,mean)
phat.zinb.enmalle<-predprob(enmalle.zinb1)
phat.zinb.mn.enmalle<-apply(phat.zinb.enmalle,2,mean)

with(enmalle.pouting,{
  hist(Ntot,prob=TRUE,col="grey",breaks=seq(-0.5,1450,1),xlab="",
       main="",ylim=c(0,0.5),xlim=c(-0.5,30))
  lines(x=seq(0,1448,1),y=phat.pois.mn.enmalle,type="l",lwd=2,col="red")
  lines(x=seq(0,1448,1),y=phat.nb.mn.enmalle,type="l",lwd=2,col="blue")
  lines(x=seq(0,1448,1),y=phat.zipoiss.mn.enmalle,type="l",lwd=2,col="yellow")
  lines(x=seq(0,1448,1),y=phat.zinb.mn.enmalle,type="l",lwd=2,col="green")
})

#8.4.6# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

aics.nasa<-AIC(nasa.poiss1,nasa.negbin1,k=2) # AIC
bics.nasa<-AIC(nasa.poiss1,nasa.negbin1,k=log(dim(nasa.pouting)[1])) # BIC
aics.enmalle<-AIC(enmalle.poiss1,enmalle.negbin1,enmalle.zipoiss1,enmalle.zinb1,k=2) # AIC
bics.enmalle<-AIC(enmalle.poiss1,enmalle.negbin1,enmalle.zipoiss1,enmalle.zinb1,k=log(dim(enmalle.pouting)[1])) # BIC

fm.nasa2<-list("ML-Pois"=nasa.poiss1,"NB"=nasa.negbin1)
logliks.nasa2<-rbind(logLik=sapply(fm.nasa2,function (x) round(logLik(x),digits=0)),
               Df=sapply(fm.nasa2,function (x) attr(logLik(x),"df")))
fm.enmalle2<-list("ML-Pois"=enmalle.poiss1,"NB"=enmalle.negbin1,"ZIP"=enmalle.zipoiss1,"ZINB"=enmalle.zinb1)
logliks.enmalle2<-rbind(logLik=sapply(fm.enmalle2,function (x) round(logLik(x),digits=0)),
                     Df=sapply(fm.enmalle2,function (x) attr(logLik(x),"df")))

zeroes.nasa<-round(c("Obs"=sum(nasa.pouting$Ntot<1),
                "ML-Pois"=sum(dpois(0,fitted(nasa.poiss1))),
                "NB"=sum(dnbinom(0,mu=fitted(nasa.negbin1),size=nasa.negbin1$theta))))
zeroes.enmalle<-round(c("Obs"=sum(enmalle.pouting$Ntot<1),
                     "ML-Pois"=sum(dpois(0,fitted(enmalle.poiss1))),
                     "NB"=sum(dnbinom(0,mu=fitted(enmalle.negbin1),size=enmalle.negbin1$theta)),
                     "ZIP"=sum(predict(enmalle.zipoiss1,type="prob")[,1]),
                     "ZINB"=sum(predict(enmalle.zinb1,type="prob")[,1])))

modelsComp.nasa<-data.frame(cbind(logliks.nasa2[c(2,4)],bics.nasa[[2]],aics.nasa[[2]],logliks.nasa2[c(1,3)],zeroes.nasa[2:3]))
colnames(modelsComp.nasa)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes.nasa[[1]],")"))
modelsComp.nasa$deltaAIC<-round(modelsComp.nasa$AIC-min(modelsComp.nasa$AIC),2)
modelsComp.nasa$deltaBIC<-round(modelsComp.nasa$BIC-min(modelsComp.nasa$BIC),2)
modelsComp.nasa$ModLikelihood<-round(exp(-modelsComp.nasa$deltaAIC/2),2)
modelsComp.nasa$AICweight<-round(modelsComp.nasa$ModLikelihood/sum(modelsComp.nasa$ModLikelihood),2)
modelsComp.nasa

modelsComp.enmalle<-data.frame(cbind(logliks.enmalle2[c(2,4,6,8)],bics.enmalle[[2]],aics.enmalle[[2]],logliks.enmalle2[c(1,3,5,7)],zeroes.enmalle[2:5]))
colnames(modelsComp.enmalle)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes.enmalle[[1]],")"))
modelsComp.enmalle$deltaAIC<-round(modelsComp.enmalle$AIC-min(modelsComp.enmalle$AIC),2)
modelsComp.enmalle$deltaBIC<-round(modelsComp.enmalle$BIC-min(modelsComp.enmalle$BIC),2)
modelsComp.enmalle$ModLikelihood<-round(exp(-modelsComp.enmalle$deltaAIC/2),2)
modelsComp.enmalle$AICweight<-round(modelsComp.enmalle$ModLikelihood/sum(modelsComp.enmalle$ModLikelihood),2)
modelsComp.enmalle

#final candidate models

#NO random effects
nasa.negbin1
enmalle.zinb1

#8.5# Selection of predictors

#8.5.1# nasa model selection
nasa.NB1<-nasa.negbin1
summary(nasa.NB1)
Anova(nasa.NB1)

# rm Depth
nasa.NB2<-glm.nb(Ntot~offset(log(offs5))+
                      lGRT+
                      fyear+poly(Julian,2)+
                      sstM+
                      caladoNight,
                    data=nasa.pouting)
summary(nasa.NB2)
Anova(nasa.NB2)
lrtest(nasa.NB1,nasa.NB2)
AICtab(nasa.NB1,nasa.NB2)

# rm poly(Julian, 2)
nasa.NB3<-glm.nb(Ntot~offset(log(offs5))+
                      lGRT+
                      fyear+
                      sstM+
                      caladoNight,
                    data=nasa.pouting)
summary(nasa.NB3)
Anova(nasa.NB3)
lrtest(nasa.NB2,nasa.NB3)
AICtab(nasa.NB2,nasa.NB3)

# rm lGRT
nasa.NB4<-glm.nb(Ntot~offset(log(offs5))+
                      fyear+
                      sstM+
                      caladoNight,
                    data=nasa.pouting)
summary(nasa.NB4)
Anova(nasa.NB4)
lrtest(nasa.NB3,nasa.NB4)
AICtab(nasa.NB3,nasa.NB4)

# rm sstM
nasa.NB5<-glm.nb(Ntot~offset(log(offs5))+
                      fyear+
                      caladoNight,
                    data=nasa.pouting)
summary(nasa.NB5)
Anova(nasa.NB5)
lrtest(nasa.NB4,nasa.NB5)
AICtab(nasa.NB4,nasa.NB5)

AICtab(nasa.NB1,nasa.NB2,nasa.NB3,nasa.NB4,nasa.NB5)
BICtab(nasa.NB1,nasa.NB2,nasa.NB3,nasa.NB4,nasa.NB5)

coefplot2(nasa.NB2)
sum(residuals(nasa.NB2,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.NB2))+1))

#comparing models
AICnasa<-AIC(nasa.NB1,nasa.NB2,nasa.NB3,nasa.NB4,nasa.NB5,k=2) # AIC
BICnasa<-AIC(nasa.NB1,nasa.NB2,nasa.NB3,nasa.NB4,nasa.NB5,k=log(dim(nasa.pouting)[1])) # BIC

FMnasa<-list("M1"=nasa.NB1,"M2"=nasa.NB2,"M3"=nasa.NB3,"M4"=nasa.NB4,"M4"=nasa.NB4)
Logliksnasa<-rbind(logLik=sapply(FMnasa,function (x) round(logLik(x),digits=0)),
                     Df=sapply(FMnasa,function (x) attr(logLik(x),"df")))

ZEROnasa<-round(c("Obs"=sum(nasa.pouting$Ntot<1),
                  "M1"=sum(dnbinom(0,mu=fitted(nasa.NB1),size=nasa.NB1$theta)),
                  "M2"=sum(dnbinom(0,mu=fitted(nasa.NB2),size=nasa.NB2$theta)),
                  "M3"=sum(dnbinom(0,mu=fitted(nasa.NB3),size=nasa.NB3$theta)),
                  "M4"=sum(dnbinom(0,mu=fitted(nasa.NB4),size=nasa.NB4$theta)),
                  "M5"=sum(dnbinom(0,mu=fitted(nasa.NB5),size=nasa.NB5$theta))
                  ))

MCompnasa<-data.frame(cbind(Logliksnasa[c(2,4,6,8,10)],BICnasa[[2]],AICnasa[[2]],Logliksnasa[c(1,3,5,7,9)],ZEROnasa[2:6]))
colnames(MCompnasa)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes.nasa[[1]],")"))
MCompnasa$deltaAIC<-round(MCompnasa$AIC-min(MCompnasa$AIC),2)
MCompnasa$deltaBIC<-round(MCompnasa$BIC-min(MCompnasa$BIC),2)
MCompnasa$ModLikelihood<-round(exp(-MCompnasa$deltaAIC/2),2)
MCompnasa$AICweight<-round(MCompnasa$ModLikelihood/sum(MCompnasa$ModLikelihood),2)
MCompnasa

#8.5.2# enmalle model selection
enmalle.ZINB1 <- enmalle.zinb1
summary(enmalle.ZINB1)

# rm sstM
enmalle.ZINB2 <- zeroinfl(Ntot~offset(log(offs1))+
                            Gear*lGRT+
                            fyear*fZoneO+
                            poly(Julian,2)+
                            Gear*lDepth+
                            Gear*caladoNight+
                            Seafloor|
                            lDepth,
                          dist="negbin",link="logit",data=enmalle.pouting)
summary(enmalle.ZINB2)
lrtest(enmalle.ZINB1,enmalle.ZINB2)
AICtab(enmalle.ZINB1,enmalle.ZINB2)

# rm ploy(2)
enmalle.ZINB3 <- zeroinfl(Ntot~offset(log(offs1))+
                            Gear*lGRT+
                            fyear*fZoneO+
                            Julian+
                            Gear*lDepth+
                            Gear*caladoNight+
                            Seafloor|
                            lDepth,
                          dist="negbin",link="logit",data=enmalle.pouting)
summary(enmalle.ZINB3)
lrtest(enmalle.ZINB2,enmalle.ZINB3)
AICtab(enmalle.ZINB2,enmalle.ZINB3)

# rm GearVETAS:caladoNight
enmalle.ZINB4 <- zeroinfl(Ntot~offset(log(offs1))+
                            Gear*lGRT+
                            fyear*fZoneO+
                            Julian+
                            Gear*lDepth+
                            caladoNight+
                            Seafloor|
                            lDepth,
                          dist="negbin",link="logit",data=enmalle.pouting)
summary(enmalle.ZINB4)
lrtest(enmalle.ZINB3,enmalle.ZINB4)
AICtab(enmalle.ZINB3,enmalle.ZINB4)

# rm GearVETAS:caladoNight
enmalle.ZINB5 <- zeroinfl(Ntot~offset(log(offs1))+
                            Gear*lGRT+
                            fyear*fZoneO+
                            Julian+
                            lDepth+
                            caladoNight+
                            Seafloor|
                            lDepth,
                          dist="negbin",link="logit",data=enmalle.pouting)
summary(enmalle.ZINB5)
lrtest(enmalle.ZINB4,enmalle.ZINB5)
AICtab(enmalle.ZINB4,enmalle.ZINB5)

#comparing models
AICenmalle<-AIC(enmalle.ZINB1,enmalle.ZINB2,enmalle.ZINB3,enmalle.ZINB4,enmalle.ZINB5,k=2) # AIC
BICenmalle<-AIC(enmalle.ZINB1,enmalle.ZINB2,enmalle.ZINB3,enmalle.ZINB4,enmalle.ZINB5,k=log(dim(enmalle.pouting)[1])) # BIC

FMenmalle<-list("M1"=enmalle.ZINB1,"M2"=enmalle.ZINB2,"M3"=enmalle.ZINB3,"M4"=enmalle.ZINB4,"M5"=enmalle.ZINB5)
Logliksenmalle<-rbind(logLik=sapply(FMenmalle,function (x) round(logLik(x),digits=0)),
                        Df=sapply(FMenmalle,function (x) attr(logLik(x),"df")))

ZEROenmalle<-round(c("Obs"=sum(enmalle.pouting$Ntot<1),
                     "M1"=sum(predict(enmalle.ZINB1,type="prob")[,1]),
                     "M2"=sum(predict(enmalle.ZINB2,type="prob")[,1]),
                     "M3"=sum(predict(enmalle.ZINB3,type="prob")[,1]),
                     "M4"=sum(predict(enmalle.ZINB4,type="prob")[,1]),
                     "M5"=sum(predict(enmalle.ZINB5,type="prob")[,1])))

MCompenmalle<-data.frame(cbind(Logliksenmalle[c(2,4,6,8,10)],BICenmalle[[2]],AICenmalle[[2]],Logliksenmalle[c(1,3,5,7,9)],ZEROenmalle[2:6]))
colnames(MCompenmalle)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes.enmalle[[1]],")"))
MCompenmalle$deltaAIC<-round(MCompenmalle$AIC-min(MCompenmalle$AIC),2)
MCompenmalle$deltaBIC<-round(MCompenmalle$BIC-min(MCompenmalle$BIC),2)
MCompenmalle$ModLikelihood<-round(exp(-MCompenmalle$deltaAIC/2),2)
MCompenmalle$AICweight<-round(MCompenmalle$ModLikelihood/sum(MCompenmalle$ModLikelihood),2)
MCompenmalle

MCompnasa
MCompenmalle

nasa.NB2
enmalle.ZINB3

#8.6# Basic model checking

nasa.pouting$residNB<-residuals(nasa.NB2,type="pearson")
nasa.pouting$fitNB<-fitted(nasa.NB2,type="response")

enmalle.pouting$residZINB<-residuals(enmalle.ZINB3,type="pearson")
enmalle.pouting$fitZINB<-fitted(enmalle.ZINB3,type="response")

#8.6.1# Normality (not really relevant) and heterogeneity

#nasa
par(mfrow=c(1,2))
hist(nasa.pouting$residNB,col="grey",breaks=50) # More negative residuals
plot(residNB~log(fitNB),data=nasa.pouting,xlab="Fitted values",ylab="Pearson residuals",pch=16)
abline(h=0,col="grey")

resPlot1.nasa<-ggplot(data=nasa.pouting,aes(x=log(fitNB),y=residNB))+
	geom_point(col="#3182bd",size=2)+
	scale_y_continuous("Pearson residuals")+
	scale_x_continuous("Fitted values")+
	ggtitle("Residual vs. Fitted") # No appaerent lack of fit
resPlot1.nasa

#enmalle
par(mfrow=c(1,2))
hist(enmalle.pouting$residZINB,col="grey",breaks=50) # More negative residuals
plot(residZINB~log(fitZINB),data=enmalle.pouting,xlab="Fitted values",ylab="Pearson residuals",pch=16)
abline(h=0,col="grey")

resPlot1.enmalle<-ggplot(data=enmalle.pouting,aes(x=log(fitZINB),y=residZINB))+
  geom_point(col="#3182bd",size=2)+
  scale_y_continuous("Pearson residuals")+
  scale_x_continuous("Fitted values")+
  ggtitle("Residual vs. Fitted") # No appaerent lack of fit
resPlot1.enmalle

#8.6.2# Residuals vs predictors

#nasa
par(mfcol=c(2,6))
boxplot(residNB~fgear,data=nasa.pouting)
abline(h=0,col="red")
boxplot(residNB~fcrew,data=nasa.pouting)
abline(h=0,col="red")
boxplot(residNB~fyear,data=nasa.pouting)
abline(h=0,col="red")
boxplot(residNB~fZoneO,data=nasa.pouting)
abline(h=0,col="red")
boxplot(residNB~Seafloor,data=nasa.pouting)
abline(h=0,col="red")
plot(residNB~lGRT,data=nasa.pouting)
abline(h=0,col="red")
plot(residNB~Julian,data=nasa.pouting)
abline(h=0,col="red")
plot(residNB~lDepth,data=nasa.pouting)
abline(h=0,col="red")
plot(residNB~QxM,data=nasa.pouting)
abline(h=0,col="red")
plot(residNB~sstM,data=nasa.pouting)
abline(h=0,col="red")
plot(residNB~caladoNight,data=nasa.pouting)
abline(h=0,col="red")
boxplot(residNB~Idflota,data=nasa.pouting)
abline(h=0,col="red")

#enmalle
par(mfcol=c(2,6))
boxplot(residZINB~fgear,data=enmalle.pouting)
abline(h=0,col="red")
boxplot(residZINB~fcrew,data=enmalle.pouting)
abline(h=0,col="red")
boxplot(residZINB~fyear,data=enmalle.pouting)
abline(h=0,col="red")
boxplot(residZINB~fZoneO,data=enmalle.pouting)
abline(h=0,col="red")
boxplot(residZINB~Seafloor,data=enmalle.pouting)
abline(h=0,col="red")
plot(residZINB~lGRT,data=enmalle.pouting)
abline(h=0,col="red")
plot(residZINB~Julian,data=enmalle.pouting)
abline(h=0,col="red")
plot(residZINB~lDepth,data=enmalle.pouting)
abline(h=0,col="red")
plot(residZINB~QxM,data=enmalle.pouting)
abline(h=0,col="red")
plot(residZINB~sstM,data=enmalle.pouting)
abline(h=0,col="red")
plot(residZINB~caladoNight,data=enmalle.pouting)
abline(h=0,col="red")
boxplot(residZINB~Idflota,data=enmalle.pouting)
abline(h=0,col="red")

#8.6.3# Spatial patterns

#nasa
ggplot(data=nasa.pouting,aes(x=Lon,y=Lat))+
	geom_point(aes(colour=residNB))+
	scale_colour_gradient(low="blue",high="red")

nasa.pouting.Spat<-data.frame(nasa.pouting$residNB,nasa.pouting$Lon,nasa.pouting$Lat)
coordinates(nasa.pouting.Spat)<- ~nasa.pouting.Lon+nasa.pouting.Lat
vario1.nasa<-variogram(nasa.pouting.residNB~1,data=nasa.pouting.Spat)
plot(vario1.nasa,pch=16,col=1,cex=1.5) # No spatial autocorrelation using the default
# classical method of moments estimate of the variogram (this assumes normality)
vario2.nasa<-variogram(nasa.pouting.residNB~1,data=nasa.pouting.Spat,cressie=TRUE)
plot(vario2.nasa,pch=16,col=1,cex=1.5) # No clear spatial autocorrelation using the
# robust variogram estimate

var1.nasa<-ggplot(data=vario1.nasa,aes(x=dist,y=gamma))+
	geom_point(aes(size=np),col="#3182bd")+
	scale_y_continuous("Sample variogram",limits=c(0.5,1.75))+
	scale_x_continuous("Distance")+
	scale_size_continuous("Number of \npoint pairs")+
	ggtitle("Classical variogram")+
	theme(legend.position=c(0,0),legend.justification=c(-0.3,-0.15))
var1.nasa
	
var2.nasa<-ggplot(data=vario2.nasa,aes(x=dist,y=gamma))+
	geom_point(aes(size=np),col="#3182bd")+
	scale_y_continuous("Sample variogram",limits=c(0.5,1))+
	scale_x_continuous("Distance")+
	scale_size_continuous("Number of \npoint pairs")+
	ggtitle("Robust variogram")+
	theme(legend.position="none")
var2.nasa

pdf(file="NBresid-nasapouting.pdf",width=13,height=4)

multiplot(resPlot1.nasa,var1.nasa,var2.nasa,cols=3)

dev.off()

#enmalle
ggplot(data=enmalle.pouting,aes(x=Lon,y=Lat))+
  geom_point(aes(colour=residZINB))+
  scale_colour_gradient(low="blue",high="red")

enmalle.pouting.Spat<-data.frame(enmalle.pouting$residZINB,enmalle.pouting$Lon,enmalle.pouting$Lat)
coordinates(enmalle.pouting.Spat)<- ~enmalle.pouting.Lon+enmalle.pouting.Lat
vario1.enmalle<-variogram(enmalle.pouting.residZINB~1,data=enmalle.pouting.Spat)
plot(vario1.enmalle,pch=16,col=1,cex=1.5) # No spatial autocorrelation using the default
# classical method of moments estimate of the variogram (this assumes normality)
vario2.enmalle<-variogram(enmalle.pouting.residZINB~1,data=enmalle.pouting.Spat,cressie=TRUE)
plot(vario2.enmalle,pch=16,col=1,cex=1.5) # No clear spatial autocorrelation using the
# robust variogram estimate

var1.enmalle<-ggplot(data=vario1.enmalle,aes(x=dist,y=gamma))+
  geom_point(aes(size=np),col="#3182bd")+
  scale_y_continuous("Sample variogram",limits=c(0.5,1.75))+
  scale_x_continuous("Distance")+
  scale_size_continuous("Number of \npoint pairs")+
  ggtitle("Classical variogram")+
  theme(legend.position=c(0,0),legend.justification=c(-0.3,-0.15))
var1.enmalle

var2.enmalle<-ggplot(data=vario2.enmalle,aes(x=dist,y=gamma))+
  geom_point(aes(size=np),col="#3182bd")+
  scale_y_continuous("Sample variogram",limits=c(0.3,0.4))+
  scale_x_continuous("Distance")+
  scale_size_continuous("Number of \npoint pairs")+
  ggtitle("Robust variogram")+
  theme(legend.position="none")
var2.enmalle

pdf(file="ZINBresid-enmallepouting.pdf",width=13,height=4)

multiplot(resPlot1.enmalle,var1.enmalle,var2.enmalle,cols=3)

dev.off()

#8.7# Model coefficients interpretation (care with ln-transformation of Depth)??

nasa.NB2
enmalle.ZINB3

exp(coef(nasa.NB2))

exp(coef(enmalle.ZINB3)[c(1:55)]) #count part
exp(coef(enmalle.ZINB3)[c(56,57)]) #zero part

# ----------------------------- #
#9# Bootstrapping the optimal model coefficients (zero & count parts)

#9.1# Function (add starting values to the model if needed!!)

nasa.NB2<-glm.nb(Ntot~offset(log(offs5))+
                   lGRT+
                   fyear+poly(Julian,2)+
                   sstM+
                   caladoNight,
                 data=nasa.pouting)
#dput(round(coef(nasa.NB2),4))

enmalle.ZINB3 <- zeroinfl(Ntot~offset(log(offs1))+
                            Gear*lGRT+
                            fyear*fZoneO+
                            Julian+
                            Gear*lDepth+
                            Gear*caladoNight+
                            Seafloor|
                            lDepth,
                          dist="negbin",link="logit",data=enmalle.pouting)
#dput(round(coef(enmalle.ZINB3,"count"),4))
#dput(round(coef(enmalle.ZINB3,"zero"),4))

boot.nb.nasa<-function (data,i) {
	
	try(mod<-glm.nb(Ntot~offset(log(offs5))+
	                  lGRT+
	                  fyear+poly(Julian,2)+
	                  sstM+
	                  caladoNight,
	                data=data[i,]))
						   			  
	if (exists("mod")) { 
						   			 
	as.vector(t(coef(summary(mod))[,1:2]))
	
	} else {rep(NA,times=2)} # If the above model crashes this fills in the gaps
    						  # with NA and the algorithm continues
                  # times=length(as.vector(t(coef(summary(mod))[,1:2])))/2
		
}

boot.zinb.enmalle<-function (data,i) {
  
  try(mod<-zeroinfl(Ntot~offset(log(offs1))+
                      Gear*lGRT+
                      fyear*fZoneO+
                      Julian+
                      Gear*lDepth+
                      Gear*caladoNight+
                      Seafloor|
                      lDepth,
                    dist="negbin",link="logit",data=data[i,]))
  
  if (exists("mod")) { 
    
    as.vector(t(do.call(rbind,coef(summary(mod)))[,1:2]))
    
  } else {rep(NA,times=58)} #times=length(as.vector(t(do.call(rbind,coef(summary(enmalle.ZINB3)))[,1:2])))/2 
  
}

#9.2# Coefficients (obtain CI of estimates excluding SE and theta)

#nasa
RR<-10000 # Number of resamples (reduce this number for testing the code)
nasa.nb.boot.out<-boot(data=nasa.pouting,statistic=boot.nb.nasa,R=RR)
nasa.nb.boot.out  # Basic output

parms.nasa<-t(sapply(c(1,3,5,7,9,11,13,15,17,19,21,23),
                     function (i) {
                       out<-boot.ci(nasa.nb.boot.out,index=c(i,i+1),type=c("perc","bca"))
                       with(out,c(Estimate=t0,pLow=percent[4],pUpp=percent[5],bcaLow=bca[4],bcaUpp=bca[5]))
                       })) # Care!!! BCA calculation sometimes crashes (if RR is small)... 

row.names(parms.nasa)<-names(coef(nasa.NB2))
parms.nasa

summary(nasa.NB2)
confint(nasa.NB2)

#enmalle
RR<-10000 # Number of resamples (reduce this number for testing the code)
enmalle.zinb.boot.out<-boot(data=enmalle.pouting,statistic=boot.zinb.enmalle,R=RR)
enmalle.zinb.boot.out  # Basic output

parms.enmalle<-t(sapply(c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57),
                function (i) {
                  out<-boot.ci(zinb.boot.out,index=c(i,i+1),type=c("perc","bca"))
                  with(out,c(Estimate=t0,pLow=percent[4],pUpp=percent[5],bcaLow=bca[4],bcaUpp=bca[5]))
                  }))
row.names(parms.enmalle)<-names(coef(enmalle.ZINB3))
parms.enmalle

summary(enmalle.ZINB3)
confint(enmalle.ZINB3)

#9.3# Example of histogram and qqplot for a given component

plot(nasa.nb.boot.out,index=1)
plot(enmalle.zinb.boot.out,index=1)

#9.4# Histograms of all components

nasa.nb.boot.out2<-nasa.nb.boot.out$t[,c(1,3,5,7,9,11,13,15,17,19,21,23)] # Coefficients of interest from the boot object matrix
par(mfrow=c(4,3))
for (i in 1:12) {
  hist(nasa.nb.boot.out2[,i],breaks=50,col="light blue",main="",
      xlab=names(coef(nasa.NB2))[i])
  abline(v=coef(nasa.NB2)[i],col="red",lwd=2)
  }

enmalle.nb.boot.out2<-nasa.nb.boot.out$t[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57)] # Coefficients of interest from the boot object matrix
par(mfrow=c(6,10))
for (i in 1:58) {
  hist(enmalle.zinb.boot.out[,i],breaks=50,col="light blue",main="",
       xlab=names(coef(enmalle.ZINB3))[i])
  abline(v=coef(enmalle.ZINB3)[i],col="red",lwd=2)
}


# ----------------------------- #
#10# Sketching results for the optimal model

#10.1# plotting nasa model
names(coef(nasa.NB2))
#[1] "(Intercept)"      "lGRT"             "fyear2005"        "fyear2006"        "fyear2007"        "fyear2008"       
#[7] "fyear2009"        "fyear2012"        "poly(Julian, 2)1" "poly(Julian, 2)2" "sstM"             "caladoNight"  

#10.1.1# Continuous variables

nbnasa1<-predict(nasa.NB2,
            newdata=data.frame(lGRT=seq(min(nasa.pouting$lGRT,na.rm=T),
                                        max(nasa.pouting$lGRT,na.rm=T),length=100), # lGRT
                               fyear="2006",
                               Julian=183,
                               sstM=0,
                               caladoNight=mean(nasa.pouting$caladoNight,na.rm=T),
                               offs5=50))

nasaGRT<-data.frame(nbnasa1,seq(min(nasa.pouting$lGRT),max(nasa.pouting$lGRT),length=100))
colnames(nasaGRT)<-c("nbnasa1","GRTSeq")
nasaGRT$GRT<-exp(nasaSST$SSTSeq)
head(nasaGRT)

nasaf1<-ggplot(data=nasaGRT,aes(x=GRT,y=nbnasa1))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(4.8,5.25))+
  scale_x_continuous("GRT",limits=c(0,5))
nasaf1

nbnasa2<-predict(nasa.NB2,
                 newdata=data.frame(lGRT=mean(nasa.pouting$lGRT,na.rm=T),
                                    fyear="2006",
                                    Julian=183,
                                    sstM=seq(min(nasa.pouting$sstM),max(nasa.pouting$sstM),length=100), #sstM
                                    caladoNight=mean(nasa.pouting$caladoNight,na.rm=T),
                                    offs5=50))

nasaSST<-data.frame(nbnasa2,seq(min(nasa.pouting$sstM),max(nasa.pouting$sstM),length=100))
colnames(nasaSST)<-c("nbnasa2","SSTSeq")
head(nasaSST)

nasaf2<-ggplot(data=nasaSST,aes(x=SSTSeq,y=nbnasa2))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(4.85,5.25))+
  scale_x_continuous("res-SST",limits=c(-2.5,1.55))
nasaf2

nbnasa3<-predict(nasa.NB2,
                 newdata=data.frame(lGRT=mean(nasa.pouting$lGRT,na.rm=T),
                                    fyear="2006",
                                    Julian=183,
                                    sstM=0, 
                                    caladoNight=seq(min(nasa.pouting$caladoNight),max(nasa.pouting$caladoNight),length=100), #caladonight
                                    offs5=50))

nasaCalN<-data.frame(nbnasa3,seq(min(nasa.pouting$caladoNight),max(nasa.pouting$caladoNight),length=100))
colnames(nasaCalN)<-c("nbnasa3","CalNSeq")
head(nasaCalN)

nasaf3<-ggplot(data=nasaCalN,aes(x=CalNSeq,y=nbnasa3))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(2.45,5.35))+
  scale_x_continuous("% Night operation",limits=c(0,1))
nasaf3

nbnasa4<-predict(nasa.NB2,
                 newdata=data.frame(lGRT=mean(nasa.pouting$lGRT,na.rm=T),
                                    fyear="2006",
                                    Julian=seq(1,365,1), #Julian
                                    sstM=0, 
                                    caladoNight=mean(nasa.pouting$caladoNight,na.rm=T),
                                    offs5=50))

nasaJulian<-data.frame(nbnasa4,seq(1,365,1))
colnames(nasaJulian)<-c("nbnasa4","JulianSeq")
head(nasaJulian)

nasaf4<-ggplot(data=nasaJulian,aes(x=JulianSeq,y=nbnasa4))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(5,5.85))+
  scale_x_continuous("Day of the Year",limits=c(1,365))
nasaf4

pdf(file="NBcont-nasapouting.pdf",width=10,height=10)
multiplot(nasaf1,nasaf2,nasaf3,nasaf4,cols=2)
dev.off()

#10.1.2# Categorical variables

nasa.newTrend<-data.frame(lGRT=rep(mean(nasa.pouting$lGRT,na.rm=T),times=7),
                     fyear=c("2002","2005","2006","2007","2008","2009","2012"),
                     Julian=rep(183,times=7),
                     sstM=rep(0,times=7),
                     caladoNight=rep(mean(nasa.pouting$caladoNight,na.rm=T),times=7),
                     offs5=rep(50,times=7))

nasa.abundInd<-predict(nasa.NB2,newdata=nasa.newTrend)
years<-c(2002,2005,2006,2007,2008,2009,2012)

nasa.poutingAbund<-data.frame(cbind(nasa.abundInd,years))
colnames(nasa.poutingAbund)<-c("Index","Year")
str(nasa.poutingAbund)

pdf(file="NBTrend-nasapouting.pdf",width=10,height=7)

ggplot(data=nasa.poutingAbund,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  scale_y_continuous("Standardized Index",limits=c(3.5,5.5))+
  scale_x_continuous("Year",limits=c(1999,2013))

dev.off()

nasapouting.abund<-lm(Index~Year,data=nasa.poutingAbund)
summary(nasapouting.abund) #decresing trend!

#10.2# plotting enmalle model
names(coef(enmalle.ZINB3))
#[1] "count_(Intercept)"           "count_GearVETAS"             "count_lGRT"                  "count_fyear2000"            
#[5] "count_fyear2001"             "count_fyear2002"             "count_fyear2003"             "count_fyear2004"            
#[9] "count_fyear2005"             "count_fyear2006"             "count_fyear2007"             "count_fyear2008"            
#[13] "count_fyear2009"             "count_fyear2010"             "count_fyear2011"             "count_fyear2012"            
#[17] "count_fyear2013"             "count_fZoneO2"               "count_fZoneO3"               "count_Julian"               
#[21] "count_lDepth"                "count_caladoNight"           "count_Seafloormixed"         "count_Seafloorsoft"         
#[25] "count_GearVETAS:lGRT"        "count_fyear2000:fZoneO2"     "count_fyear2001:fZoneO2"     "count_fyear2002:fZoneO2"    
#[29] "count_fyear2003:fZoneO2"     "count_fyear2004:fZoneO2"     "count_fyear2005:fZoneO2"     "count_fyear2006:fZoneO2"    
#[33] "count_fyear2007:fZoneO2"     "count_fyear2008:fZoneO2"     "count_fyear2009:fZoneO2"     "count_fyear2010:fZoneO2"    
#[37] "count_fyear2011:fZoneO2"     "count_fyear2012:fZoneO2"     "count_fyear2013:fZoneO2"     "count_fyear2000:fZoneO3"    
#[41] "count_fyear2001:fZoneO3"     "count_fyear2002:fZoneO3"     "count_fyear2003:fZoneO3"     "count_fyear2004:fZoneO3"    
#[45] "count_fyear2005:fZoneO3"     "count_fyear2006:fZoneO3"     "count_fyear2007:fZoneO3"     "count_fyear2008:fZoneO3"    
#[49] "count_fyear2009:fZoneO3"     "count_fyear2010:fZoneO3"     "count_fyear2011:fZoneO3"     "count_fyear2012:fZoneO3"    
#[53] "count_fyear2013:fZoneO3"     "count_GearVETAS:lDepth"      "count_GearVETAS:caladoNight" "zero_(Intercept)"           
#[57] "zero_lDepth 

#10.2.1# ZI part

newDepth<-seq(min(enmalle.pouting$lDepth),max(enmalle.pouting$lDepth),length=100)
logistDepth<-coef(enmalle.ZINB3,model="zero")[1]+(coef(enmalle.ZINB3,model="zero")[2]*newDepth)
probDepth<-exp(logistDepth)/(1+exp(logistDepth))
newDepth2<-exp(newDepth)

par(mar=c(5,5,4,4))
plot(newDepth2,probDepth,xlab="ln-Depth (m)",ylab="Probability of false zeros",type="n",cex.lab=1.4,cex.axis=1.3,ylim=c(0,1))
lines(newDepth2,probDepth,col="blue",lwd=3)
grid()

#predict(zi) 
zienmalle1<-predict(enmalle.ZINB3,type="zero",
	newdata=data.frame(Gear="VETAS",
					   lGRT=mean(enmalle.pouting$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(enmalle.pouting$lDepth),
					   			  max(enmalle.pouting$lDepth),length=100),
					   caladoNight=mean(enmalle.pouting$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200)) # Estandarizacion a num por 200 m2 por h (200 m2 es el tamaÃ±o legal aproximado de un paÃ±o: 50 m de largo por 4 de alto)

prob0<-data.frame(zienmalle1,seq(min(enmalle.pouting$lDepth),max(enmalle.pouting$lDepth),length=100))
colnames(prob0)<-c("l1","DepthSeq")

pdf(file="ZIprob-enmallepouting.pdf")

f0<-ggplot(data=prob0,aes(x=DepthSeq,y=l1))
f0+geom_line(colour="blue",lwd=1)+
	scale_y_continuous("Probability of false zeros",limits=c(0,0.45))+
	scale_x_continuous("ln-Depth (m)",limits=c(0,5.5))

dev.off()

#10.2.2# Continuous variables

zinbenmalle1<-predict(enmalle.ZINB3,type="count",
                      newdata=data.frame(Gear="VETAS",
                                         lGRT=mean(enmalle.pouting$lGRT,na.rm=T),
                                         fyear="2006",
                                         Julian=seq(1,365,1), #Julian
                                         lDepth=mean(enmalle.pouting$lDepth,na.rm=T),
                                         caladoNight=mean(enmalle.pouting$caladoNight,na.rm=T),
                                         fZoneO="1",
                                         Seafloor="hard",
                                         offs1=200))

enmalleJulian<-data.frame(zinbenmalle1,seq(1,365,1))
colnames(enmalleJulian)<-c("zinbenmalle1","JulianSeq")

enmallef1<-ggplot(data=enmalleJulian,aes(x=JulianSeq,y=zinbenmalle1))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance")+
  scale_x_continuous("Day of the Year")
enmallef1 #I do not like this seasonal trend ¿lineal? Fuck off

levels(enmalle.pouting$Gear)
zinbenmalle2<-predict(enmalle.ZINB3,type="count",
                      newdata=data.frame(Gear=c(rep("VETAS",100),rep("MINOS",100)),
                                         lGRT=mean(enmalle.pouting$lGRT,na.rm=T),
                                         fyear="2006",
                                         Julian=183,
                                         lDepth=c(rep(seq(min(enmalle.pouting$lDepth),max(enmalle.pouting$lDepth),length.out=100),2)),
                                         caladoNight=mean(enmalle.pouting$caladoNight,na.rm=T),
                                         fZoneO="1",
                                         Seafloor="hard",
                                         offs1=200))

enmalleDepth<-data.frame(zinbenmalle2,
                         c(rep(seq(min(enmalle.pouting$lDepth),max(enmalle.pouting$lDepth),length.out=100),2)),
                         c(rep("VETAS",100),rep("MINOS",100)))
colnames(enmalleDepth)<-c("zinbenmalle2","DepthSeq","GearSeq")
enmalleDepth$Depth<-exp(enmalleDepth$DepthSeq)
head(enmalleDepth)

enmallef2<-ggplot(data=enmalleDepth,aes(x=Depth,y=zinbenmalle2))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance")+
  scale_x_continuous("Depth (m)")+
  facet_wrap(~GearSeq)
enmallef2

zinbenmalle3<-predict(enmalle.ZINB3,type="count",
                      newdata=data.frame(Gear=c(rep("VETAS",100),rep("MINOS",100)),
                                         lGRT=mean(enmalle.pouting$lGRT,na.rm=T),
                                         fyear="2006",
                                         Julian=183,
                                         lDepth=mean(enmalle.pouting$lDepth,na.rm=T),
                                         caladoNight=c(rep(seq(min(enmalle.pouting$caladoNight),max(enmalle.pouting$caladoNight),length.out=100),2)),
                                         fZoneO="1",
                                         Seafloor="hard",
                                         offs1=200))

enmalleCalN<-data.frame(zinbenmalle3,
                         c(rep(seq(min(enmalle.pouting$caladoNight),max(enmalle.pouting$caladoNight),length.out=100),2)),
                         c(rep("VETAS",100),rep("MINOS",100)))
colnames(enmalleCalN)<-c("zinbenmalle3","CalNSeq","GearSeq")
head(enmalleCalN)

enmallef3<-ggplot(data=enmalleCalN,aes(x=CalNSeq,y=zinbenmalle3))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance")+
  scale_x_continuous("% Night soak time")+
  facet_wrap(~GearSeq)
enmallef3

zinbenmalle4<-predict(enmalle.ZINB3,type="count",
                      newdata=data.frame(Gear=c(rep("VETAS",100),rep("MINOS",100)),
                                         lGRT=c(rep(seq(min(enmalle.pouting$caladoNight),max(enmalle.pouting$caladoNight),length.out=100),2)),
                                         fyear="2006",
                                         Julian=183,
                                         lDepth=mean(enmalle.pouting$lDepth,na.rm=T),
                                         caladoNight=mean(enmalle.pouting$lDepth,na.rm=T),
                                         fZoneO="1",
                                         Seafloor="hard",
                                         offs1=200))

enmalleGRT<-data.frame(zinbenmalle4,
                        c(rep(seq(min(enmalle.pouting$caladoNight),max(enmalle.pouting$caladoNight),length.out=100),2)),
                        c(rep("VETAS",100),rep("MINOS",100)))
colnames(enmalleGRT)<-c("zinbenmalle4","GRTSeq","GearSeq")
enmalleGRT$GRT<-exp(enmalleGRT$GRTSeq)
head(enmalleGRT)

enmallef4<-ggplot(data=enmalleGRT,aes(x=GRT,y=zinbenmalle4))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance")+
  scale_x_continuous("GRT")+
  facet_wrap(~GearSeq)
enmallef4

#10.2.3# Categorical variables

enmalle.newTrend<-data.frame(Gear=rep("VETAS",times=45),
                             lGRT=rep(mean(enmalle.pouting$lGRT,na.rm=T),times=45),
                             fyear=rep(c("1999","2000","2001","2002","2003","2004","2005",
                                     "2006","2007","2008","2009","2010","2011","2012","2013"),3),
                             Julian=rep(183,times=45),
                             lDepth=rep(mean(enmalle.pouting$lDepth,na.rm=T),times=45),
                             caladoNight=rep(mean(enmalle.pouting$caladoNight,na.rm=T),times=45),
                             fZoneO=c(rep("1",times=15),rep("2",times=15),rep("3",times=15)),
                             Seafloor=rep("hard",times=45),
                             offs1=rep(200,times=45))

enmalle.abundInd<-predict(enmalle.ZINB3,newdata=enmalle.newTrend,type="count")

enmalle.poutingAbund<-data.frame(cbind(enmalle.abundInd))
colnames(enmalle.poutingAbund)<-c("Index")
enmalle.poutingAbund$Year<-as.numeric(rep(seq(1999,2013,1),3))
enmalle.poutingAbund$Zones<-c(rep("Rias Baixas",times=15),rep("Artabro",times=15),rep("Cantabrico",times=15))
enmalle.poutingAbund$Zones <- ordered(enmalle.poutingAbund$Zones,
                                      levels = c("Rias Baixas", "Artabro", "Cantabrico"))
str(enmalle.poutingAbund)

pdf(file="ZINBTrend-enmallepouting.pdf",width=10,height=7)

ggplot(data=enmalle.poutingAbund,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)

dev.off()

enmallepouting.abund<-lm(Index~Year+Zones,data=enmalle.poutingAbund)
summary(enmallepouting.abund)
anova(enmallepouting.abund) #no decreasing trend, only zones diferencies

newSeafloor<-data.frame(Gear=rep("VETAS",times=3),
                        lGRT=rep(mean(enmalle.pouting$lGRT,na.rm=T),times=3),
                        fyear=rep("2006",times=3),
                        Julian=rep(183,times=3),
                        lDepth=rep(mean(enmalle.pouting$lDepth,na.rm=T),times=3),
                        caladoNight=rep(mean(enmalle.pouting$caladoNight,na.rm=T),times=3),
                        fZoneO=rep("1",times=3),
                        Seafloor=c("hard","mixed","soft"),
                        offs1=rep(200,times=3))

seafloorP<-predict(enmalle.ZINB3,newdata=newSeafloor,type="count")
vals<-c(1,2,3)
probSeafloor<-data.frame(cbind(seafloorP,vals))

enmalleSFT<-ggplot(data=probSeafloor,aes(x=vals,y=seafloorP))+
  geom_bar(stat="identity",fill="blue",col="blue")+
  scale_y_continuous("")+scale_x_continuous("Seafloor type",breaks=c(1,2,3),
                     labels=c("Hard","Mixed","Soft"))
enmalleSFT

pdf(file="ZINB-enmallepouting.pdf",width=20,height=30)
multiplot(enmallef4,enmallef2,enmallef3,enmallef1,enmalleSFT,cols=2)
dev.off()

#11# Plot of abundance, nominal cpue, and landings

#11.1# Calculate nominal cpue and average for the Rias Baixas

#nasa
head(nasa.pouting)

nasa.pouting$cpue<-(nasa.pouting$Ntot*50)/nasa.pouting$offs5 # Standardize at 50 pieces
cpues.nasa<-tapply(nasa.pouting$cpue,nasa.pouting$Year,mean,na.rm=T)
cpues.L.nasa<-tapply(nasa.pouting$cpue,nasa.pouting$Year,ci95Low)
cpues.U.nasa<-tapply(nasa.pouting$cpue,nasa.pouting$Year,ci95Upp)

cpues.nasa.pouting<-data.frame(cbind(cpues.nasa,cpues.L.nasa,cpues.U.nasa,c(2002,2005,2006,2007,2008,2009,2012)))
colnames(cpues.nasa.pouting)<-c("Index","ciLow","ciUpp","Year")
str(cpues.nasa.pouting)

limits <- aes(ymax = ciUpp, ymin= ciLow)
ggplot(data=cpues.nasa.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  geom_errorbar(limits, colour="blue", width=0)+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))

#enmalle
head(enmalle.pouting)

enmalle.pouting$cpue<-(enmalle.pouting$Ntot*200)/enmalle.pouting$offs1 # Standardize at 200 m2 paño

enmalleRB<-enmalle.pouting[enmalle.pouting$ZoneO==1,]
enmalleAR<-enmalle.pouting[enmalle.pouting$ZoneO==2,]
enmalleCN<-enmalle.pouting[enmalle.pouting$ZoneO==3,]

cpues.enmalle.RB<-tapply(enmalleRB$cpue,enmalleRB$Year,mean,na.rm=T)
cpues.L.enmalle.RB<-tapply(enmalleRB$cpue,enmalleRB$Year,ci95Low)
cpues.U.enmalle.RB<-tapply(enmalleRB$cpue,enmalleRB$Year,ci95Upp)
cpues.enmalle.RB<-data.frame(cbind(cpues.enmalle.RB,cpues.L.enmalle.RB,cpues.U.enmalle.RB))
cpues.enmalle.RB$Zones<-rep("Rias Baixas",15)
cpues.enmalle.RB$Year<-as.numeric(seq(1999,2013,1))
colnames(cpues.enmalle.RB)<-c("Index","ciLow","ciUpp","Zones","Year")
str(cpues.enmalle.RB)

cpues.enmalle.AR<-tapply(enmalleAR$cpue,enmalleAR$Year,mean,na.rm=T)
cpues.L.enmalle.AR<-tapply(enmalleAR$cpue,enmalleAR$Year,ci95Low)
cpues.U.enmalle.AR<-tapply(enmalleAR$cpue,enmalleAR$Year,ci95Upp)
cpues.enmalle.AR<-data.frame(cbind(cpues.enmalle.AR,cpues.L.enmalle.AR,cpues.U.enmalle.AR))
cpues.enmalle.AR$Zones<-rep("Artabro",15)
cpues.enmalle.AR$Year<-as.numeric(seq(1999,2013,1))
colnames(cpues.enmalle.AR)<-c("Index","ciLow","ciUpp","Zones","Year")
str(cpues.enmalle.AR)

cpues.enmalle.CN<-tapply(enmalleCN$cpue,enmalleCN$Year,mean,na.rm=T)
cpues.L.enmalle.CN<-tapply(enmalleCN$cpue,enmalleCN$Year,ci95Low)
cpues.U.enmalle.CN<-tapply(enmalleCN$cpue,enmalleCN$Year,ci95Upp)
cpues.enmalle.CN<-data.frame(cbind(cpues.enmalle.CN,cpues.L.enmalle.CN,cpues.U.enmalle.CN))
cpues.enmalle.CN$Zones<-rep("Cantabrico",15)
cpues.enmalle.CN$Year<-as.numeric(seq(1999,2013,1))
colnames(cpues.enmalle.CN)<-c("Index","ciLow","ciUpp","Zones","Year")
str(cpues.enmalle.CN)

cpues.enmalle.pouting<-data.frame(rbind(cpues.enmalle.RB,cpues.enmalle.AR,cpues.enmalle.CN))
str(cpues.enmalle.pouting)
cpues.enmalle.pouting
cpues.enmalle.pouting$Zones <- ordered(cpues.enmalle.pouting$Zones,levels = c("Rias Baixas", "Artabro", "Cantabrico"))

limits <- aes(ymax = ciUpp, ymin= ciLow)
ggplot(data=cpues.enmalle.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  geom_errorbar(limits, colour="blue", width=0)+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year")+
  facet_wrap(~Zones)

#11.2# Official landings

landings.pouting<-read.table(file="faneca_landings.txt",header=T,dec=".")
landings.pouting
str(landings.pouting)
zone1<-landings.pouting[,c(1,12)]
colnames(zone1)<-c("Year","Index")
zone2<-landings.pouting[,c(1,13)]
colnames(zone2)<-c("Year","Index")
zone3<-landings.pouting[,c(1,14)]
colnames(zone3)<-c("Year","Index")
zones<-rbind(zone1,zone2,zone3)
landings.pouting<-data.frame(cbind(zones,c(rep("Rias Baixas",18),rep("Artabro",18),rep("Cantabrico",18))))
colnames(landings.pouting)<-c("Year","Index","Zones")
str(landings.pouting)
head(landings.pouting)
landings.pouting$Zones <- ordered(landings.pouting$Zones,levels = c("Rias Baixas", "Artabro", "Cantabrico"))

ggplot(data=landings.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)

#11.3# Comparing abundance index, nominal cpue, landings

t1a<-ggplot(data=nasa.poutingAbund,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  ggtitle("Standardize abundance index (trap)")+
  theme(plot.title = element_text(size=26, face="bold"))
t1a

t1b<-ggplot(data=enmalle.poutingAbund,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)+
  ggtitle("Standardize abundance index (gillnet)")+
  theme(plot.title = element_text(size=26, face="bold"))
t1b

t2a<-ggplot(data=cpues.nasa.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  #geom_errorbar(limits, colour="blue", width=0)+
  scale_y_continuous("Nominal CPUE (nº/50pieces*h)")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  ggtitle("Nominal CPUE (trap)")+
  theme(plot.title = element_text(size=26, face="bold"))
t2a

t2b<-ggplot(data=cpues.enmalle.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  #geom_errorbar(limits, colour="blue", width=0)+
  scale_y_continuous("Nominal CPUE (nº/200m2*h)")+
  scale_x_continuous("Year")+
  facet_wrap(~Zones)+  ggtitle("Nominal CPUE (gillnet)")+
  theme(plot.title = element_text(size=26, face="bold"))
t2b

t3<-ggplot(data=landings.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="blue")+
  scale_y_continuous("Total landings (kg)")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)+facet_wrap(~Zones)+
  ggtitle("Oficial landings")+
  theme(plot.title = element_text(size=26, face="bold"))
t3

multiplot(t1a,t2a,t3,t1b,t2b,t3,cols=2)

pdf(file="pouting_trends.pdf",width=20,height=15)

multiplot(t1a,t2a,t3,t1b,t2b,t3,cols=2)

dev.off()

#correlations trends

names(enmalle.poutingAbund)
names(cpues.enmalle.pouting)
names(landings.pouting)
landings.pouting2<-landings.pouting[landings.pouting$Year>1998&landings.pouting$Year<2014,]

#Rias Baixas
cor(enmalle.poutingAbund[enmalle.poutingAbund$Zones=="Rias Baixas",]$Index,
    cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Rias Baixas",]$Index,method="spearman")
cor(enmalle.poutingAbund[enmalle.poutingAbund$Zones=="Rias Baixas",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Rias Baixas",]$Index,method="spearman")
cor(cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Rias Baixas",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Rias Baixas",]$Index,method="spearman")

#Artabro
cor(enmalle.poutingAbund[enmalle.poutingAbund$Zones=="Artabro",]$Index,
    cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Artabro",]$Index,method="spearman")
cor(enmalle.poutingAbund[enmalle.poutingAbund$Zones=="Artabro",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Artabro",]$Index,method="spearman")
cor(cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Artabro",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Artabro",]$Index,method="spearman")

#Cantabrico
cor(enmalle.poutingAbund[enmalle.poutingAbund$Zones=="Cantabrico",]$Index,
    cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Cantabrico",]$Index,method="spearman")
cor(enmalle.poutingAbund[enmalle.poutingAbund$Zones=="Cantabrico",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Cantabrico",]$Index,method="spearman")
cor(cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Cantabrico",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Cantabrico",]$Index,method="spearman")

#ridiculous correlation with oficial landings :(

# ----------------------------- #
#12# Plotting ZINB model predicted means for YEAR with Bootstrapping
# of predictions (95% CI) by means of shuffling residuals (Thierry Onkelinx code)

#12.1# Bootstrap

newTrend2<-enmalle.newTrend
head(newTrend2)

Fit<-predict(enmalle.ZINB3,type="response")

Pearson<-residuals(enmalle.ZINB3,type="pearson") # Pearson residuals
VarComp<-residuals(enmalle.ZINB3,type="response")/Pearson # Raw residuals/Pearson residuals

Gear<-enmalle.pouting$Gear
lGRT<-enmalle.pouting$lGRT
fyear<-enmalle.pouting$fyear
Julian<-enmalle.pouting$Julian
lDepth<-enmalle.pouting$lDepth
caladoNight<-enmalle.pouting$caladoNight
fZoneO<-enmalle.pouting$fZoneO
Seafloor<-enmalle.pouting$Seafloor
offs1<-enmalle.pouting$offs1

RR<-10 # Number of resamples (reduce this number for testing the code)

bootstrap<-replicate(n=RR,{ 
	
    yStar<-pmax(round(Fit+sample(Pearson)*VarComp,0),0)
    try(mod<-zeroinfl(yStar~offset(log(offs1))+
                        Gear*lGRT+
                        fyear*fZoneO+
                        Julian+
                        Gear*lDepth+
                        Gear*caladoNight+
                        Seafloor|
                        lDepth,
                      dist="negbin",link="logit"))
    
    if (exists("mod")) {
    	
    	predict(mod,newdata=newTrend2,type="response")
    	
    } else {rep(NA,times=15)} # If the above model crashes this fills in the gaps
    						  # with NA and the algorithm continues

})

CIs<-t(apply(X=bootstrap,MARGIN=1,FUN=quantile,c(0.025,0.975),na.rm=T))
colnames(CIs)<-c("ciLow","ciUpp")
newTrend2$fit<-predict(enmalle.ZINB3,newdata=newTrend2,type="response")
newTrend2<-cbind(newTrend2,CIs)
newTrend2$Year<-seq(1999,2013,1)

#12.2# Plot of abundance
levels(newTrend2$fZoneO)<-c("1"="Rias Baixas", "2"="Artabro","3"="Cantabrico")
levels(newTrend2$fZoneO)

pdf(file="pouting_YearPredCI.pdf",width=10,height=8)

ggplot(data=newTrend2,aes(x=Year,y=fit))+
	geom_segment(aes(x=Year,y=fit,xend=Year,yend=ciUpp),col="blue",lwd=0.5)+
	geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=fit),col="blue",lwd=0.5)+
	geom_line(lwd=0.5,linetype="dotted",col="blue")+
	geom_point(size=5,col="gray90")+
	geom_point(size=3,col="blue")+
	scale_y_continuous("Standardized Index of Abundance")+
	scale_x_continuous("Year",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))+
  facet_wrap(~fZoneO)

dev.off()

boot1<-ggplot(data=newTrend2,aes(x=Year,y=fit))+
  geom_segment(aes(x=Year,y=fit,xend=Year,yend=ciUpp),col="blue",lwd=0.5)+
  geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=fit),col="blue",lwd=0.5)+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="blue")+
  scale_y_continuous("Standardized Index of Abundance")+
  scale_x_continuous("Year")+
  facet_wrap(~fZoneO,scales="free_y")+
  ggtitle("Standardized abundance index")+
  theme(plot.title = element_text(size=26, face="bold"))
boot1

cpues.enmalle.pouting2<-cpues.enmalle.pouting
min(cpues.enmalle.pouting2$ciLow)
cpues.enmalle.pouting2$ciLow<-ifelse(cpues.enmalle.pouting2$ciLow<0,0,cpues.enmalle.pouting2$ciLow)
limits <- aes(ymax = ciUpp, ymin= ciLow)
cpue1<-ggplot(data=cpues.enmalle.pouting2,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="blue")+
  geom_errorbar(limits, colour="blue", width=0)+
  scale_y_continuous("Nominal CPUE (nº/200m2*h)")+
  scale_x_continuous("Year")+
  facet_wrap(~Zones,scales="free_y")+
  ggtitle("Nominal CPUE")+
  theme(plot.title = element_text(size=26, face="bold"))
cpue1

landing1<-ggplot(data=landings.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="blue")+
  scale_y_continuous("Total landings (kg)")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones,scales="free_y")+
  ggtitle("Oficial landings")+
  theme(plot.title = element_text(size=26, face="bold"))
landing1

multiplot(boot1,cpue1,landing1,cols=1)

pdf(file="pouting-enmalle_trends.pdf",width=15,height=20)

multiplot(boot1,cpue1,landing1,cols=1)

dev.off()
