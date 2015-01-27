##########################################################
# Analysis of Trisopterus luscus catch rates in Galicia#
# using data sampled onboard fishing vessels by the UTPB #
##########################################################

#Date start: April 2014 by Jaime Otero (bulk of code from pollachius analysis)
#Last update: 27-1-2015 by Alex Alonso
#many code string were carried out by Jaime Otero during POollachius pollachius analysis
#input data comes from "preparation_data_base_casgass.R" code that generate faneca.txt file
#General objective: to analyse variation of catch rates from monitoring program data of UTPB
#
# Project: CASGASS
#
# Authors: Jaime Otero (jotero@iim.csic.es) and Alex Alonso (alex@iim.csic.es)
# Further modifications by: Alexandre Alonso (alex@iim.csic.es)
# Institute of Marine Research (Vigo)

#set working directory for input and output
setwd("C:\\Users\\alex\\Dropbox\\casgass") 

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

source(file="HighStatLib.r")
source(file="multiple_ggplot.r")
source(file="CI_mean.r")

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
	scale_y_continuous(limits=c(0,100),"Frequency")+
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
qx$QY<-qx$QY/1000 # Change scale

g1<-gam(QX~s(DoY,k=6,bs="cc")+s(Timer,k=4,bs="cr"),data=qx,na.action=na.exclude)

summary(g1)
#plot(g1,pages=1)

g2<-gam(QY~s(DoY,k=6,bs="cc")+s(Timer,k=4,bs="cr"),data=qx,na.action=na.exclude)

summary(g2)
#plot(g2,pages=1)

qx$qxAno<-residuals(g1)
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
enmalle.pouting$lGRT<-log(enmalle.pouting$GRT)
enmalle.pouting$lDepth<-log(enmalle.pouting$Depth)

#VIF
vifs1<-c("GRT","Crew","Year","Julian","Depth","QxM","sstM","caladoNight")
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

nasa.poiss0<-gam(Ntot~offset(log(offs5))+
                   fcrew+s(GRT,k=3)+
                   fyear+s(Julian,by=Depth,k=6,bs="cc")+
                   s(QxM,k=3)+s(sstM,k=3)+
                   caladoNight+
                   fZoneO+Seafloor,
                 family=poisson,data=nasa.pouting)

summary(nasa.poiss0)
sum(residuals(nasa.poiss0,type="pearson")^2)/nasa.poiss0$df.resid # Overdispersion
plot(nasa.poiss0,pages=1,all.terms=T,scale=0,shade="T")

#8.1.2# GLM Poisson modelling

nasa.poiss1<-glm(Ntot~offset(log(offs5))+
				 fcrew+GRT+
			 	 fyear+ns(Julian,3)+
			 	 Depth+ns(Julian,3):Depth+
			 	 QxM+sstM+
			 	 caladoNight+
			 	 fZoneO+Seafloor,
			 	 family=poisson,data=nasa.pouting)

summary(nasa.poiss1)
sum(residuals(nasa.poiss1,type="pearson")^2)/nasa.poiss1$df.resid # Overdispersion
plot(allEffects(nasa.poiss1))

#8.1.3# Negative Binomial modelling

nasa.negbin1<-glm.nb(Ntot~offset(log(offs5))+
          fcrew+GRT+
          fyear+ns(Julian,2)+
          Depth+ns(Julian,2):Depth+
          QxM+sstM+
          caladoNight+
          fZoneO+Seafloor,
			 		 data=nasa.pouting)

summary(nasa.negbin1)
sum(residuals(nasa.negbin1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin1))+1)) 
hist(residuals(nasa.negbin1,type="deviance"))
plot(residuals(nasa.negbin1,type="deviance")~log(fitted(nasa.negbin1)))
boxplot(residuals(nasa.negbin1,type="deviance")~nasa.pouting$Idflota)
abline(0,0,col="red")
  table(nasa.pouting$Idjornada,nasa.pouting$Idflota)
plot(allEffects(nasa.negbin1))
plot(effect("ns(Julian,2):Depth",nasa.negbin1))
plot(effect("Julian",nasa.negbin1))
plot(effect("Depth",nasa.negbin1))

#8.1.3.1. Mixed model NB

#http://glmmadmb.r-forge.r-project.org/glmmADMB.html

nasa.pouting$Idflota<-as.factor(nasa.pouting$Idflota)
nasa.pouting$Nasaoffset<-log(nasa.pouting$offs5)

nasa.negbinran1 <- glmmadmb(Ntot~fcrew+GRT+
                              fyear+ns(Julian,2)+
                              Depth+ns(Julian,2):Depth+
                              QxM+sstM+
                              caladoNight+
                              fZoneO+Seafloor+
                              offset(Nasaoffset)+
                              (1|Idflota), 
                            data=nasa.pouting,
                         zeroInflation=FALSE, 
                         family="nbinom")
summary(nasa.negbinran1)
sum(residuals(nasa.negbinran1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbinran1))+1)) 
Anova(nasa.negbinran1)
hist(residuals(nasa.negbinran1,type="pearson"))
plot(residuals(nasa.negbinran1,type="pearson")~log(fitted(nasa.negbinran1)))
boxplot(residuals(nasa.negbinran1,type="pearson")~nasa.pouting$Idflota)
abline(0,0,col="red")

#8.1.4# comparing model approaches

AICtab(nasa.poiss0,nasa.poiss1,nasa.negbinran1,nasa.negbin1)
BICtab(nasa.poiss0,nasa.poiss1,nasa.negbinran1,nasa.negbin1)

#mejoramos el modelo con los random effects, what shall we do?

#8.2#ENMALLE modelling

#8.2.1# GAM Poisson modelling

enmalle.poiss0<-gam(Ntot~offset(log(offs1))+
                   Gear*Depth+
                   fcrew+s(GRT,k=3)+
                   fyear*fZoneO+
                     s(Julian,by=Depth,k=6,bs="cc")+
                   s(QxM,k=3)+s(sstM,k=3)+
                   Gear*caladoNight+
                   Seafloor,
                 family=poisson,data=enmalle.pouting)

summary(enmalle.poiss0)
sum(residuals(enmalle.poiss0,type="pearson")^2)/enmalle.poiss0$df.resid # Overdispersion
plot(enmalle.poiss0,pages=1,all.terms=T,scale=0,shade="T")

#8.2.2# GLM Poisson modelling

enmalle.poiss1<-glm(Ntot~offset(log(offs1))+
                   Gear*Depth+
                   fcrew+GRT+
                   fyear*fZoneO+
                   ns(Julian,3)+
                   Depth+ns(Julian,3):Depth+
                   QxM+sstM+
                   Gear*caladoNight+
                   Seafloor,
                 family=poisson,data=enmalle.pouting)

summary(enmalle.poiss1)
sum(residuals(enmalle.poiss1,type="pearson")^2)/enmalle.poiss1$df.resid # Overdispersion
plot(allEffects(enmalle.poiss1))

#8.2.3# Negative Binomial modelling

enmalle.negbin1<-glm.nb(Ntot~offset(log(offs1))+
                          Gear*Depth+
                          fcrew+GRT+
                          fyear*fZoneO+
                          ns(Julian,3)+
                          Depth+ns(Julian,3):Depth+
                          QxM+sstM+
                          Gear*caladoNight+
                          Seafloor,
                     data=enmalle.pouting)

summary(enmalle.negbin1)
Anova(enmalle.negbin1)
sum(residuals(enmalle.negbin1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbin1))+1)) 
hist(residuals(enmalle.negbin1,type="deviance"))
plot(residuals(enmalle.negbin1,type="deviance")~log(fitted(enmalle.negbin1)))
boxplot(residuals(enmalle.negbin1,type="deviance")~enmalle.pouting$Idflota)
abline(0,0,col="red")
plot(allEffects(enmalle.negbin1))
plot(effect("Julian",enmalle.negbin1))
plot(effect("Depth",enmalle.negbin1))

#8.2.3.1. Mixed model NB

#http://glmmadmb.r-forge.r-project.org/glmmADMB.html

enmalle.pouting$Idflota<-as.factor(enmalle.pouting$Idflota)
enmalle.pouting$Enmalleoffset<-log(enmalle.pouting$offs1)

enmalle.negbinran1 <- glmmadmb(Ntot~Gear*Depth+
                              fcrew+GRT+
                              fyear*fZoneO+
                              ns(Julian,3)+
                              Depth+ns(Julian,3):Depth+
                              QxM+sstM+
                              Gear*caladoNight+
                              Seafloor+
                              offset(Enmalleoffset)+
                              (1|Idflota), 
                            data=enmalle.pouting,
                            zeroInflation=FALSE, 
                            family="nbinom")
summary(enmalle.negbinran1)
sum(residuals(enmalle.negbinran1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbinran1))+1)) 
Anova(enmalle.negbinran1)
hist(residuals(enmalle.negbinran1,type="pearson"))
plot(residuals(enmalle.negbinran1,type="pearson")~log(fitted(enmalle.negbinran1)))
boxplot(residuals(enmalle.negbinran1,type="pearson")~enmalle.pouting$Idflota)
abline(0,0,col="red")

#8.2.4# comparing model approaches

AICtab(enmalle.poiss0,enmalle.poiss1,enmalle.negbinran1,enmalle.negbin1)
BICtab(enmalle.poiss0,enmalle.poiss1,enmalle.negbinran1,enmalle.negbin1)

#mejoramos el modelo con los random effects, what shall we do?

#según los parámetros de sobredispersión el mixed negative binomial es el más adecuado
# parece que no será necesario recurrir a los zero inflated en ninguno de los casos

#8.3# en cualquier caso probamos con un modelo zero inflated simple para ambos tipos de arte

#8.3.1# nasa
nasa.zipoiss1<-zeroinfl(Ntot~offset(log(offs5))+
                     fcrew+GRT+
                     fyear+ns(Julian,2)+
                     Depth+ns(Julian,2):Depth+
                     QxM+sstM+
                     caladoNight+
                     fZoneO+Seafloor|
                       Depth,
                   dist="poisson",link="logit",data=nasa.pouting)

summary(nasa.zipoiss1)
sum(residuals(nasa.zipoiss1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zipoiss1))))

nasa.zinb1<-zeroinfl(Ntot~offset(log(offs5))+
                       fcrew+GRT+
                       fyear+ns(Julian,2)+
                       Depth+ns(Julian,2):Depth+
                       QxM+sstM+
                       caladoNight+
                       fZoneO+Seafloor|
                       Depth,
                dist="negbin",link="logit",data=nasa.pouting)

summary(nasa.zinb1)
sum(residuals(nasa.zinb1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zinb1))+1))

nasa.zinbran1 <- glmmadmb(Ntot~fcrew+GRT+
                              fyear+ns(Julian,2)+
                              Depth+ns(Julian,2):Depth+
                              QxM+sstM+
                              caladoNight+
                              fZoneO+Seafloor+
                              offset(Nasaoffset)+
                              (1|Idflota), 
                            data=nasa.pouting,
                            zeroInflation=TRUE, 
                            family="nbinom")
summary(nasa.zinbran1)
sum(residuals(nasa.zinbran1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zinbran1))+1)) 
Anova(nasa.zinbran1)
hist(residuals(nasa.zinbran1,type="pearson"))
plot(residuals(nasa.zinbran1,type="pearson")~log(fitted(nasa.zinbran1)))
boxplot(residuals(nasa.zinbran1,type="pearson")~nasa.pouting$Idflota)
abline(0,0,col="red")

#8.3.2# enmalle
enmalle.zipoiss1<-zeroinfl(Ntot~offset(log(offs1))+
                          fcrew+GRT+
                          fyear+ns(Julian,2)+
                          Depth+ns(Julian,2):Depth+
                          QxM+sstM+
                          caladoNight+
                          fZoneO+Seafloor|
                          Depth,
                        dist="poisson",link="logit",data=enmalle.pouting)

summary(enmalle.zipoiss1)
sum(residuals(enmalle.zipoiss1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zipoiss1)))) # Overdispersed

enmalle.zinb1<-zeroinfl(Ntot~offset(log(offs1))+
                       fcrew+GRT+
                       fyear+ns(Julian,2)+
                       Depth+ns(Julian,2):Depth+
                       QxM+sstM+
                       caladoNight+
                       fZoneO+Seafloor|
                       Depth,
                     dist="negbin",link="logit",data=enmalle.pouting)

summary(enmalle.zinb1)
sum(residuals(enmalle.zinb1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zinb1))+1)) # Slightly overdispersed

#8.4# full model comparison

#8.4.1# AIC/BIC

#AIC
AICtab(nasa.poiss0,nasa.poiss1,nasa.negbinran1,nasa.negbin1,nasa.zipoiss1,nasa.zinb1)
AICtab(enmalle.poiss0,enmalle.poiss1,enmalle.negbinran1,enmalle.negbin1,enmalle.zipoiss1,enmalle.zinb1)

#BIC
nasa.bic=BIC(nasa.poiss0,nasa.poiss1,nasa.negbinran1,nasa.negbin1,nasa.zipoiss1,nasa.zinb1)
nasa.bic[order(nasa.bic$BIC), ]
enmalle.bic=BIC(enmalle.poiss0,enmalle.poiss1,enmalle.negbinran1,enmalle.negbin1,enmalle.zipoiss1,enmalle.zinb1)
enmalle.bic[order(enmalle.bic$BIC), ]

#8.4.2# overdispersion parameter

#nasa
nasa.GAMpoisson=sum(residuals(nasa.poiss0,type="pearson")^2)/nasa.poiss0$df.resid
nasa.GLMpoisson=sum(residuals(nasa.poiss1,type="pearson")^2)/nasa.poiss1$df.resid
nasa.NB=sum(residuals(nasa.negbin1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin1))+1)) 
nasa.MixedNB=sum(residuals(nasa.negbinran1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbinran1))+1)) 
nasa.ZIpoisson=sum(residuals(nasa.zipoiss1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zipoiss1))))
nasa.ZINB=sum(residuals(nasa.zinb1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zinb1))+1))

nasa.model=c("nasa.GAMpoisson","nasa.GLMpoisson","nasa.NB","nasa.MixedNB","nasa.ZIpoisson","nasa.ZINB")
nasa.theta=c(nasa.GAMpoisson,nasa.GLMpoisson,nasa.NB,nasa.MixedNB,nasa.ZIpoisson,nasa.ZINB)
nasa.M=data.frame(cbind(nasa.model,nasa.theta))
nasa.M=nasa.M[order(nasa.M$nasa.theta),]

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

# Hurdle modelling (assumption: two step model)
#we do not considered to include this in the model comparison
#but we keep it in the code just to see other posibilities

#nasa
nasa.hurdle1<-hurdle(Ntot~offset(log(offs5))+fcrew+GRT+fyear+ns(Julian,2)+Depth+ns(Julian,2):Depth+
                       QxM+sstM+caladoNight+fZoneO+Seafloor|Depth,
                     dist="poisson",link="logit",data=nasa.pouting)
summary(nasa.hurdle1)

nasa.hurdle2<-hurdle(Ntot~offset(log(offs5))+fcrew+GRT+fyear+ns(Julian,2)+Depth+ns(Julian,2):Depth+
                       QxM+sstM+caladoNight+fZoneO+Seafloor|Depth,
                     dist="negbin",link="logit",data=nasa.pouting)
summary(nasa.hurdle2)
#nemalle
enmalle.hurdle1<-hurdle(Ntot~offset(log(offs1))+fcrew+GRT+fyear+ns(Julian,2)+Depth+ns(Julian,2):Depth+
                       QxM+sstM+caladoNight+fZoneO+Seafloor|Depth,
                     dist="poisson",link="logit",data=enmalle.pouting)
summary(enmalle.hurdle1)

enmalle.hurdle2<-hurdle(Ntot~offset(log(offs1))+fcrew+GRT+fyear+ns(Julian,2)+Depth+ns(Julian,2):Depth+
                       QxM+sstM+caladoNight+fZoneO+Seafloor|Depth,
                     dist="negbin",link="logit",data=enmalle.pouting)
summary(enmalle.hurdle2)



#8.7# Selection of predictors

zinb2<-zeroinfl(Ntot~offset(log(offs1))+
				   	 fgear+
				   	 fcrew+
				   	 fyear+poly(Julian,2)+
				   	 lDepth+
				   	 QxM+sstM+
				   	 poly(caladoNight,2)+
				   	 fZoneO+Seafloor |
				   	 lDepth,
				   	 dist="negbin",link="logit",data=pollack1)

lrtest(zinb1,zinb2) # Removing of lGRT

zinb3<-zeroinfl(Ntot~offset(log(offs1))+
				   	 fgear+
				   	 fcrew+
				   	 fyear+poly(Julian,2)+
				   	 lDepth+
				   	 sstM+
				   	 poly(caladoNight,2)+
				   	 fZoneO+Seafloor |
				   	 lDepth,
				   	 dist="negbin",link="logit",data=pollack1)

lrtest(zinb2,zinb3) # Removing of QxM
summary(zinb3)

zinb4<-zeroinfl(Ntot~offset(log(offs1))+
				   	 fgear+
				   	 fcrew+
				   	 fyear+poly(Julian,2)+
				   	 lDepth+
				   	 poly(caladoNight,2)+
				   	 fZoneO+Seafloor |
				   	 lDepth,
				   	 dist="negbin",link="logit",data=pollack1)

lrtest(zinb3,zinb4) # Removing of sstM
summary(zinb4)

#8.8# Basic model checking

pollack1$residZINB<-residuals(zinb1,type="pearson")
pollack1$fitZINB<-fitted(zinb1,type="response")

#8.8.1# Normality (not really relevant) and heterogeneity

hist(pollack1$residZINB) # More negative residuals
plot(residZINB~fitZINB,data=pollack1,xlab="Fitted values",ylab="Pearson residuals")

resPlot1<-ggplot(data=pollack1,aes(x=log(fitZINB),y=residZINB))+
	geom_point(col="#3182bd",size=2)+
	scale_y_continuous("Pearson residuals")+
	scale_x_continuous("Fitted values")+
	ggtitle("Residual vs. Fitted") # No appaerent lack of fit

#8.8.2# Residuals vs predictors

par(mfcol=c(2,6))
boxplot(residZINB~fgear,data=pollack1)
boxplot(residZINB~fcrew,data=pollack1)
boxplot(residZINB~fyear,data=pollack1)
boxplot(residZINB~fZoneO,data=pollack1)
boxplot(residZINB~Seafloor,data=pollack1)
plot(residZINB~lGRT,data=pollack1)
plot(residZINB~Julian,data=pollack1)
plot(residZINB~lDepth,data=pollack1)
plot(residZINB~QxM,data=pollack1)
plot(residZINB~sstM,data=pollack1)
plot(residZINB~caladoNight,data=pollack1)

#8.8.3# Spatial patterns

ggplot(data=pollack1,aes(x=Lon,y=Lat))+
	geom_point(aes(colour=residZINB))+
	scale_colour_gradient(low="blue",high="red")

pollack.Spat<-data.frame(pollack1$residZINB,pollack1$Lon,pollack1$Lat)
coordinates(pollack.Spat)<- ~pollack1.Lon+pollack1.Lat
vario1<-variogram(pollack1.residZINB~1,data=pollack.Spat)
# plot(vario1,pch=16,col=1,cex=1.5) # No spatial autocorrelation using the default
# classical method of moments estimate of the variogram (this assumes normality)
vario2<-variogram(pollack1.residZINB~1,data=pollack.Spat,cressie=TRUE)
# plot(vario2,pch=16,col=1,cex=1.5) # No clear spatial autocorrelation using the
# robust variogram estimate

var1<-ggplot(data=vario1,aes(x=dist,y=gamma))+
	geom_point(aes(size=np),col="#3182bd")+
	scale_y_continuous("Sample variogram",limits=c(0,2))+
	scale_x_continuous("Distance")+
	scale_size_continuous("Number of \npoint pairs")+
	ggtitle("Classical variogram")+
	theme(legend.position=c(0,0),legend.justification=c(-0.1,0))
	
var2<-ggplot(data=vario2,aes(x=dist,y=gamma))+
	geom_point(aes(size=np),col="#3182bd")+
	scale_y_continuous("Sample variogram",limits=c(0,0.3))+
	scale_x_continuous("Distance")+
	scale_size_continuous("Number of \npoint pairs")+
	ggtitle("Robust variogram")+
	theme(legend.position="none")

pdf(file="/Users/jaimeoterovillar/Desktop/ZIresid.pdf",width=13,height=4)

multiplot(resPlot1,var1,var2,cols=3)

dev.off()

#8.9# Model selection and predicted probabilities

#8.9.1# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

aics<-AIC(poiss1,negbin1,hurdle1,hurdle2,zipoiss1,zinb1,k=2) # AIC
bics<-AIC(poiss1,negbin1,hurdle1,hurdle2,zipoiss1,zinb1,k=log(dim(pollack1)[1])) # BIC

fm<-list("ML-Pois"=poiss1,"NB"=negbin1,"Hurdle-P"=hurdle1,"Hurdle-NB"=hurdle2,"ZIP"=zipoiss1,"ZINB"=zinb1)

logliks<-rbind(logLik=sapply(fm,function (x) round(logLik(x),digits=0)),
	Df=sapply(fm,function (x) attr(logLik(x),"df"))) 

zeroes<-round(c("Obs"=sum(pollack1$Ntot<1),
		"ML-Pois"=sum(dpois(0,fitted(poiss1))),
		"NB"=sum(dnbinom(0,mu=fitted(negbin1),size=negbin1$theta)),
		"Hurdle-P"=sum(predict(hurdle1,type="prob")[,1]),
		"Hurdle-NB"=sum(predict(hurdle2,type="prob")[,1]),
		"ZIP"=sum(predict(zipoiss1,type="prob")[,1]),
		"ZINB"=sum(predict(zinb1,type="prob")[,1])))

modelsComp<-data.frame(cbind(logliks[c(2,4,6,8,10,12)],bics[[2]],aics[[2]],logliks[c(1,3,5,7,9,11)],zeroes[2:7]))
colnames(modelsComp)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes[[1]],")"))

modelsComp$deltaAIC<-round(modelsComp$AIC-min(modelsComp$AIC),2)
modelsComp$deltaBIC<-round(modelsComp$BIC-min(modelsComp$BIC),2)
modelsComp$ModLikelihood<-round(exp(-modelsComp$deltaAIC/2),2)
modelsComp$AICweight<-round(modelsComp$ModLikelihood/sum(modelsComp$ModLikelihood),2)
modelsComp

#8.9.2# Vuong test to compare NB vs ZINB

vuong(negbin1,zinb1) # zinb model is sligthly better than negbin model

#8.9.3# Predicting probabilities

phat.pois<-predprob(poiss1)
phat.pois.mn<-apply(phat.pois,2,mean)
phat.nb<-predprob(negbin1)
phat.nb.mn<-apply(phat.nb,2,mean) 

with(pollack1,{
  hist(Ntot,prob=TRUE,col="yellow",breaks=seq(-0.5,247.5,1),xlab="",
     main="")
  lines(x=seq(0,247,1),y=phat.pois.mn,type="b",lwd=2,col="red")
  lines(x=seq(0,247,1),y=phat.nb.mn,type="b",lwd=2,col="blue")
})

#8.10# Model coefficients interpretation (care with ln-transformation of Depth)??

exp(coef(zinb1)[c(30,31)]) # Zero-inflated part
# Baseline odds of no catching a pollack is 4.077108e-06; the odds is
# increased by one unit increase in Depth by 1.278464e+01 (for a one-unit
# change in Depth there is a ~12% increase in the odds of non-catching a pollack)

exp(coef(zinb1)[c(1:29)]) # Count part
# The baseline number of pollack catch is 2.421547e+01. A unit increase in
# Depth decrease 0.48 times the expected catch of pollack (for a one-unit
# change in Depth there is a ~52% decrease in the expected catch of pollack)
# Fishing in soft seafloor decreases the expected catch by ~73% 100*(2.709858e-01-1)...


# ----------------------------- #
#9# Bootstrapping the optimal model coefficients (zero & count parts)

#9.1# Function (add starting values to the model if needed!!)

#dput(round(coef(zinb1,"count"),4))
#dput(round(coef(zinb1,"zero"),4))

boot.zinb<-function (data,i) {
	
	try(mod<-zeroinfl(Ntot~offset(log(offs1))+fgear+fcrew+lGRT+
						   fyear+poly(Julian,2)+lDepth+QxM+sstM+
						   poly(caladoNight,2)+fZoneO+Seafloor |
						   lDepth,data=data[i,],dist="negbin",link="logit"))
						   			  
	if (exists("mod")) { 
						   			 
	as.vector(t(do.call(rbind,coef(summary(mod)))[,1:2]))
	
	} else {rep(NA,times=64)} # If the above model crashes this fills in the gaps
    						  # with NA and the algorithm continues
		
}

#9.2# Coefficients (obtain CI of estimates excluding SE and theta)

RR<-10000 # Number of resamples (reduce this number for testing the code)
zinb.boot.out<-boot(data=pollack1,statistic=boot.zinb,R=RR,parallel="snow",ncpus=4)
zinb.boot.out  # Basic output

parms<-t(sapply(c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,
	37,39,41,43,45,47,49,51,53,55,57,61,63),function (i) {
	out<-boot.ci(zinb.boot.out,index=c(i,i+1),type=c("perc","bca"))
	with(out,c(Estimate=t0,pLow=percent[4],pUpp=percent[5],
	bcaLow=bca[4],bcaUpp=bca[5]))
})) # Care!!! BCA calculation sometimes crashes (if RR is small)... 

row.names(parms)<-names(coef(zinb1))
parms

summary(zinb1)
confint(zinb1)

#9.3# Example of histogram and qqplot for a given component

plot(zinb.boot.out,index=1)

#9.4# Histograms of all components

zinb.boot.out2<-zinb.boot.out$t[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,
	33,35,37,39,41,43,45,47,49,51,53,55,57,61,63)] # Coefficients of interest from the boot object matrix

par(mfrow=c(4,8))
for (i in 1:31) {
  hist(zinb.boot.out2[,i],col="light blue",main="",
      xlab=names(coef(zinb1))[i])
  abline(v=coef(zinb1)[i],lwd=2)
  }


# ----------------------------- #
#10# Sketching results for the optimal model

#10.1# Plotting false zeros (binomial) curve from ZINB

#10.1.1# Hand-made

# newDepth<-seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100)

# logistDepth<-coef(zinb1,model="zero")[1]+(coef(zinb1,model="zero")[2]*newDepth)

# probDepth<-exp(logistDepth)/(1+exp(logistDepth))

# par(mar=c(5,5,4,4))
# plot(newDepth,probDepth,xlab="ln-Depth (m)",ylab="Probability of false zeros",type="n",cex.lab=1.4,cex.axis=1.3,ylim=c(0,1))
# lines(newDepth,probDepth,col="blue",lwd=3)
# grid()

#10.1.2# Using zeroinfl.predict

l1<-predict(zinb1,type="zero",
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(pollack1$lDepth),
					   			  max(pollack1$lDepth),length=100),
					   QxM=0,
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200)) # Estandarizacion a num por 200 m2 por h (200 m2 es el tamaÃ±o legal aproximado de un paÃ±o: 50 m de largo por 4 de alto)

prob0<-data.frame(l1,seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100))
colnames(prob0)<-c("l1","DepthSeq")

pdf(file="/Users/jaimeoterovillar/Desktop/ZIprob.pdf")

f0<-ggplot(data=prob0,aes(x=DepthSeq,y=l1))
f0+geom_line(colour="red",lwd=1)+
	scale_y_continuous("Probability of false zeros",limits=c(0,1))+
	scale_x_continuous("ln-Depth (m)",limits=c(0,7))

dev.off()

#10.2# Plotting count (NB) curves from ZINB (this coding could be shortend.
# However, I wrote it this way to keep track on what I'm doing!!)

#10.2.1# Continuous variables

z1<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=seq(min(pollack1$lGRT,na.rm=T),
					   			max(pollack1$lGRT,na.rm=T),length=100), # lGRT
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=0,
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probGRT<-data.frame(z1,seq(min(pollack1$lGRT),max(pollack1$lGRT),length=100))
colnames(probGRT)<-c("z1","GRTSeq")

z2<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=seq(from=1,to=366,by=1), # Julian
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=0,
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDOY<-data.frame(z2,seq(from=1,to=366,by=1))
colnames(probDOY)<-c("z2","DOYSeq")

z3<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(pollack1$lDepth,na.rm=T),
					   			  max(pollack1$lDepth,na.rm=T),length=100), # lDepth
					   QxM=0,
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDepth<-data.frame(z3,seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100))
colnames(probDepth)<-c("z3","DepthSeq")

z4<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=seq(min(pollack1$QxM,na.rm=T),
					   			max(pollack1$QxM,na.rm=T),length=100), # Qx
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probQX<-data.frame(z4,seq(min(pollack1$QxM),max(pollack1$QxM),length=100))
colnames(probQX)<-c("z4","QXSeq")

z5<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=0,
					   sstM=seq(min(pollack1$sstM,na.rm=T),
					   			max(pollack1$sstM,na.rm=T),length=100), # sst
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probSST<-data.frame(z5,seq(min(pollack1$sstM),max(pollack1$sstM),length=100))
colnames(probSST)<-c("z5","SSTSeq")

z6<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=0,
					   sstM=0,
					   caladoNight=seq(min(pollack1$caladoNight,na.rm=T),
					   			   max(pollack1$caladoNight,na.rm=T),length=100), # caladoNight
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probNight<-data.frame(z6,seq(min(pollack1$caladoNight),max(pollack1$caladoNight),length=100))
colnames(probNight)<-c("z6","NightSeq")

f1<-ggplot(data=probGRT,aes(x=GRTSeq,y=z1))+
	geom_line(colour="blue",lwd=1)+
	scale_y_continuous("Standardized Abundance",limits=c(0,0.2))+
	scale_x_continuous("ln-GRT",limits=c(-1,4))

f2<-ggplot(data=probDOY,aes(x=DOYSeq,y=z2))+
	geom_line(colour="blue",lwd=1)+
	scale_y_continuous("",limits=c(0,0.2))+
	scale_x_continuous("Day of the Year",limits=c(1,366),breaks=c(1,100,200,300))

f3<-ggplot(data=probDepth,aes(x=DepthSeq,y=z3))+
	geom_line(colour="blue",lwd=1)+
	scale_y_continuous("Standardized Abundance",limits=c(0,1.7))+
	scale_x_continuous("ln-Depth (m)",limits=c(0,7))

f4<-ggplot(data=probQX,aes(x=QXSeq,y=z4))+
	geom_line(colour="blue",lwd=1)+
	scale_y_continuous("",limits=c(0,0.2))+
	scale_x_continuous(expression(paste(italic(-Q[X])," ","(",m^3," ",s^-1," ",km^-1,") Ã—",10^3)),limits=c(-5,2))

f5<-ggplot(data=probSST,aes(x=SSTSeq,y=z5))+
	geom_line(colour="blue",lwd=1)+
	scale_y_continuous("",limits=c(0,0.2))+
	scale_x_continuous("SST (ÂºC)",limits=c(-4,3))

f6<-ggplot(data=probNight,aes(x=NightSeq,y=z6))+
	geom_line(colour="blue",lwd=1)+
	scale_y_continuous("",limits=c(0,0.5))+
	scale_x_continuous("Night soak (%)",limits=c(0,1))

pdf(file="/Users/jaimeoterovillar/Desktop/ZIprobNB.pdf",width=10,height=6)

multiplot(f1,f3,f2,f4,f6,f5,cols=3)

dev.off()

#10.2.2# Categorical variables

#10.2.2.1# Year trend

newTrend<-data.frame(fgear=rep("VETAS",times=15),fcrew=rep("2",times=15),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=15),
					fyear=c("1999","2000","2001","2002","2003","2004","2005",
						"2006","2007","2008","2009","2010","2011","2012","2013"),
					Julian=rep(183,times=15),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=15),
					QxM=rep(0,times=15),sstM=rep(0,times=15),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=15),
					fZoneO=rep("1",times=15),Seafloor=rep("hard",times=15),
					offs1=rep(200,times=15))

abundInd<-predict(zinb1,newdata=newTrend,type="count")
years<-seq(1999,2013,1)

pollackAbund<-data.frame(cbind(abundInd,years))
colnames(pollackAbund)<-c("Index","Year")

pdf(file="/Users/jaimeoterovillar/Desktop/ZIprobTrend.pdf",width=10,height=8)

ggplot(data=pollackAbund,aes(x=Year,y=Index))+
	geom_line(lwd=0.5,linetype="dotted",col="blue")+
	geom_point(size=5,col="orange")+
	scale_y_continuous("Standardized Index")+
	scale_x_continuous("Year",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

dev.off()

t.abund<-lm(Index~Year,data=pollackAbund)
summary(t.abund)

#10.2.2.2# Gear

newGear<-data.frame(fgear=c("VETAS","MINOS"),fcrew=rep("2",times=2),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=2),
					fyear=rep("2006",times=2),
					Julian=rep(183,times=2),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=2),
					QxM=rep(0,times=2),sstM=rep(0,times=2),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=2),
					fZoneO=rep("1",times=2),Seafloor=rep("hard",times=2),
					offs1=rep(200,times=2))

gearP<-predict(zinb1,newdata=newGear,type="count")
vals<-c(1,2)
probGear<-data.frame(cbind(gearP,vals))

#10.2.2.3# Crew

newCrew<-data.frame(fgear=rep("VETAS",times=2),fcrew=c("1","2"),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=2),
					fyear=rep("2006",times=2),
					Julian=rep(183,times=2),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=2),
					QxM=rep(0,times=2),sstM=rep(0,times=2),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=2),
					fZoneO=rep("1",times=2),Seafloor=rep("hard",times=2),
					offs1=rep(200,times=2))

crewP<-predict(zinb1,newdata=newCrew,type="count")
probCrew<-data.frame(cbind(crewP,vals))

#10.2.2.4# ZoneO

newZoneO<-data.frame(fgear=rep("VETAS",times=3),fcrew=rep("2",times=3),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=3),
					fyear=rep("2006",times=3),
					Julian=rep(183,times=3),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=3),
					QxM=rep(0,times=3),sstM=rep(0,times=3),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=3),
					fZoneO=c("1","2","3"),Seafloor=rep("hard",times=3),
					offs1=rep(200,times=3))

zoneP<-predict(zinb1,newdata=newZoneO,type="count")
vals<-c(1,2,3)
probZone<-data.frame(cbind(zoneP,vals))

#10.2.2.5# Seafloor

newSeafloor<-data.frame(fgear=rep("VETAS",times=3),fcrew=rep("2",times=3),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=3),
					fyear=rep("2006",times=3),
					Julian=rep(183,times=3),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=3),
					QxM=rep(0,times=3),sstM=rep(0,times=3),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=3),
					fZoneO=rep("1",times=3),Seafloor=c("hard","mixed","soft"),
					offs1=rep(200,times=3))

seafloorP<-predict(zinb1,newdata=newSeafloor,type="count")
probSeafloor<-data.frame(cbind(seafloorP,vals))

#10.2.2.6# Plots

f7<-ggplot(data=probGear,aes(x=vals,y=gearP))+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("Standardized Abundance",limits=c(0,0.2))+
	scale_x_continuous("Gear",breaks=c(1,2),
		labels=c("Veta","MiÃ±o"),limits=c(0.75,2.25))
	
f8<-ggplot(data=probCrew,aes(x=vals,y=crewP))+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("",limits=c(0,0.2))+
	scale_x_continuous("Crew",breaks=c(1,2),
		labels=c("1-3","4-6"),limits=c(0.75,2.25))
		
f9<-ggplot(data=probZone,aes(x=vals,y=zoneP))+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("Standardized Abundance",limits=c(0,0.34))+
	scale_x_continuous("Oceanographic Zone",breaks=c(1,2,3),
		labels=c("RÃ­as Baixas","Arco Ãrtabro","CantÃ¡brico"),limits=c(0.75,3.25))

f10<-ggplot(data=probSeafloor,aes(x=vals,y=seafloorP))+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("",limits=c(0,0.2))+
	scale_x_continuous("Seafloor type",breaks=c(1,2,3),
		labels=c("Hard","Mixed","Soft"),limits=c(0.75,3.25))

pdf(file="/Users/jaimeoterovillar/Desktop/ZIprobFacts.pdf")

multiplot(f7,f9,f8,f10,cols=2)

dev.off()


# ----------------------------- #
#11# Plotting ZINB model predicted means for YEAR with Bootstrapping
# of predictions (95% CI) by means of shuffling residuals (Thierry Onkelinx code)

#11.1# Bootstrap

newTrend2<-newTrend

Fit<-predict(zinb1,type="response")

Pearson<-residuals(zinb1,type="pearson") # Pearson residuals
VarComp<-residuals(zinb1,type="response")/Pearson # Raw residuals/Pearson residuals

fgear<-pollack1$fgear
fcrew<-pollack1$fcrew
lGRT<-pollack1$lGRT
fyear<-pollack1$fyear
Julian<-pollack1$Julian
lDepth<-pollack1$lDepth
QxM<-pollack1$QxM
sstM<-pollack1$sstM
caladoNight<-pollack1$caladoNight
fZoneO<-pollack1$fZoneO
Seafloor<-pollack1$Seafloor
offs1<-pollack1$offs1

RR<-10000 # Number of resamples (reduce this number for testing the code)

bootstrap<-replicate(n=RR,{ 
	
    yStar<-pmax(round(Fit+sample(Pearson)*VarComp,0),0)
    try(mod<-zeroinfl(yStar~offset(log(offs1))+fgear+fcrew+lGRT+fyear+Julian+
    						lDepth+QxM+sstM+caladoNight+fZoneO+Seafloor | lDepth,
    						dist="negbin",link="logit"))
    
    if (exists("mod")) {
    	
    	predict(mod,newdata=newTrend2,type="response")
    	
    } else {rep(NA,times=15)} # If the above model crashes this fills in the gaps
    						  # with NA and the algorithm continues

})

CIs<-t(apply(X=bootstrap,MARGIN=1,FUN=quantile,c(0.025,0.975),na.rm=T))
colnames(CIs)<-c("ciLow","ciUpp")
newTrend2$fit<-predict(zinb1,newdata=newTrend2,type="response")
newTrend2<-cbind(newTrend2,CIs)
newTrend2$Year<-seq(1999,2013,1)

#11.2# Plot of abundance

pdf(file="/Users/jaimeoterovillar/Desktop/YearPredCI.pdf",width=10,height=8)

ggplot(data=newTrend2,aes(x=Year,y=fit))+
	geom_segment(aes(x=Year,y=fit,xend=Year,yend=ciUpp),col="orange",lwd=0.5)+
	geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=fit),col="orange",lwd=0.5)+
	geom_line(lwd=0.5,linetype="dotted",col="blue")+
	geom_point(size=5,col="gray90")+
	geom_point(size=3,col="orange")+
	scale_y_continuous("Standardized Index of Abundance",limits=c(0,0.3))+
	scale_x_continuous("Year",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

dev.off()

#11.3# Plot of biomass

size.data<-pollack1[pollack1$fgear=="VETAS" & pollack1$fZoneO=="1",] # Faltarian fcrew=="2" y Seafloor=="hard" pero no hay esas categorias para 1999 y 2010-13 

size.data$Wtot.M<-size.data$Wtot/size.data$Ntot # Mean fish weight per haul (in kg)

biomass<-newTrend2$fit*tapply(size.data$Wtot.M,size.data$Year,mean,na.rm=T) # Mean body size per year (for the predicted categories)

pollackBiom<-data.frame(cbind(biomass,newTrend2$Year))
colnames(pollackBiom)<-c("Index","Year")

pdf(file="/Users/jaimeoterovillar/Desktop/BiomassTrend.pdf",width=10,height=8)

ggplot(data=pollackBiom,aes(x=Year,y=Index))+
	geom_line(lwd=0.3,linetype="dotted",col="blue")+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("Standardized Index")+
	scale_x_continuous("Year",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

dev.off()

#11.4# Plot of abundance, nominal cpue, and landings

#11.4.1# Calculate nominal cpue and average for the Rias Baixas

head(pollack1)

pollack1$cpue<-(pollack1$Ntot*200)/pollack1$offs1 # Standardize at 200 m2 hour
cpue.RB<-pollack1[pollack1$ZoneO==1,]
cpues<-tapply(cpue.RB$cpue,cpue.RB$Year,mean,na.rm=T)
cpues.L<-tapply(cpue.RB$cpue,cpue.RB$Year,ci95Low)
cpues.U<-tapply(cpue.RB$cpue,cpue.RB$Year,ci95Upp)

cpues.RB<-data.frame(cbind(cpues,rep(NA,15),rep(NA,15),seq(1999,2013,1))) # No le aÃ±ado los intervalos calculados en las lineas anteriores
colnames(cpues.RB)<-c("Index","ciLow","ciUpp","Year")

#11.4.2# Abundance fit

head(newTrend2)
pollackTrends<-newTrend2[,13:16]
colnames(pollackTrends)<-c("Index","ciLow","ciUpp","Year")

#11.4.3# Official landings

catches<-read.csv2(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/2_abadejo landings.csv",header=T,dec=".",sep=",")

catches

catches.RB<-data.frame(cbind(catches$Rias_Baixas[3:17]/1000,rep(NA,15),rep(NA,15),seq(1999,2013,1)))

colnames(catches.RB)<-c("Index","ciLow","ciUpp","Year")

cor(pollackTrends$Index,cpues.RB$Index,method="spearman")
cor(pollackTrends$Index,catches.RB$Index,method="spearman")
cor(cpues.RB$Index,catches.RB$Index,method="spearman")

#11.4.4# Plot

pollackSummary<-rbind(pollackTrends,cpues.RB,catches.RB)
pollackSummary$Trend<-rep(c("Abundance","Nominal CPUE","Official Landings"),each=15)

pdf(file="/Users/jaimeoterovillar/Desktop/PollackTrends.pdf",width=8,height=10)

ggplot(data=pollackSummary,aes(x=Year,y=Index))+
	geom_segment(aes(x=Year,y=Index,xend=Year,yend=ciUpp),col="orange",lwd=0.5)+
	geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=Index),col="orange",lwd=0.5)+
	geom_line(lwd=0.5,linetype="dotted",col="blue")+
	geom_point(size=5,col="gray90")+
	geom_point(size=3,col="orange")+
	facet_grid(Trend ~.,scales="free_y")+
	scale_y_continuous("Standardized Index")+
	scale_x_continuous("Year",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

dev.off()


# ----------------------------- #
#12# Plotting NB model results (just for comparison!!)

#12.1# Continuous variables

z1<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=seq(min(pollack1$lGRT,na.rm=T),
					   			max(pollack1$lGRT,na.rm=T),length=100), # lGRT
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=0,
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probGRT<-data.frame(z1$fit,z1$se.fit,seq(min(pollack1$lGRT),max(pollack1$lGRT),length=100))
colnames(probGRT)<-c("z1","z1SE","GRTSeq")

z2<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=seq(from=1,to=366,by=1), # Julian
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=0,
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDOY<-data.frame(z2$fit,z2$se.fit,seq(from=1,to=366,by=1))
colnames(probDOY)<-c("z2","z2SE","DOYSeq")

z3<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(pollack1$lDepth,na.rm=T),
					   			  max(pollack1$lDepth,na.rm=T),length=100), # lDepth
					   QxM=0,
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDepth<-data.frame(z3$fit,z3$se.fit,seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100))
colnames(probDepth)<-c("z3","z3SE","DepthSeq")

z4<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=seq(min(pollack1$QxM,na.rm=T),
					   			max(pollack1$QxM,na.rm=T),length=100), # Qx
					   sstM=0,
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probQX<-data.frame(z4$fit,z4$se.fit,seq(min(pollack1$QxM),max(pollack1$QxM),length=100))
colnames(probQX)<-c("z4","z4SE","QXSeq")

z5<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=0,
					   sstM=seq(min(pollack1$sstM,na.rm=T),
					   			max(pollack1$sstM,na.rm=T),length=100), # sst
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probSST<-data.frame(z5$fit,z5$se.fit,seq(min(pollack1$sstM),max(pollack1$sstM),length=100))
colnames(probSST)<-c("z5","z5SE","SSTSeq")

z6<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS",
					   fcrew="2",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   QxM=0,
					   sstM=0,
					   caladoNight=seq(min(pollack1$caladoNight,na.rm=T),
					   			   max(pollack1$caladoNight,na.rm=T),length=100), # caladoNight
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probNight<-data.frame(z6$fit,z6$se.fit,seq(min(pollack1$caladoNight),max(pollack1$caladoNight),length=100))
colnames(probNight)<-c("z6","z6SE","NightSeq")

f1<-ggplot(data=probGRT,aes(x=GRTSeq,y=exp(z1)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(z1-1.96*z1SE),ymax=exp(z1+1.96*z1SE)),
		fill="#66CCFF",alpha=0.15)+
	scale_y_continuous("Standardized Abundance",limits=c(0,0.25))+
	scale_x_continuous("ln-GRT",limits=c(-1,4))

f2<-ggplot(data=probDOY,aes(x=DOYSeq,y=exp(z2)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(z2-1.96*z2SE),ymax=exp(z2+1.96*z2SE)),
		fill="#66CCFF",alpha=0.15)+
	scale_y_continuous("",limits=c(0,0.25))+
	scale_x_continuous("Day of the Year",limits=c(1,366),breaks=c(1,100,200,300))

f3<-ggplot(data=probDepth,aes(x=DepthSeq,y=exp(z3)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(z3-1.96*z3SE),ymax=exp(z3+1.96*z3SE)),
		fill="#66CCFF",alpha=0.15)+
	scale_y_continuous("Standardized Abundance",limits=c(0,2.6))+
	scale_x_continuous("ln-Depth (m)",limits=c(0,7))

f4<-ggplot(data=probQX,aes(x=QXSeq,y=exp(z4)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(z4-1.96*z4SE),ymax=exp(z4+1.96*z4SE)),
		fill="#66CCFF",alpha=0.15)+
	scale_y_continuous("",limits=c(0,0.25))+
	scale_x_continuous(expression(paste(italic(-Q[X])," ","(",m^3," ",s^-1," ",km^-1,") Ã—",10^3)),limits=c(-5,2))

f5<-ggplot(data=probSST,aes(x=SSTSeq,y=exp(z5)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(z5-1.96*z5SE),ymax=exp(z5+1.96*z5SE)),
		fill="#66CCFF",alpha=0.15)+
	scale_y_continuous("",limits=c(0,0.25))+
	scale_x_continuous("SST (ÂºC)",limits=c(-4,3))

f6<-ggplot(data=probNight,aes(x=NightSeq,y=exp(z6)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(z6-1.96*z6SE),ymax=exp(z6+1.96*z6SE)),
		fill="#66CCFF",alpha=0.15)+
	scale_y_continuous("",limits=c(0,0.6))+
	scale_x_continuous("Night soak (%)",limits=c(0,1))

pdf(file="/Users/jaimeoterovillar/Desktop/NBprobCont.pdf",width=10,height=6)

multiplot(f1,f3,f2,f4,f6,f5,cols=3)

dev.off()

#12.2# Categorical variables

#12.2.1# Year trend

zYears<-predict(negbin1,type="link",se=T,newdata=newTrend)

abundInd<-zYears$fit
abundIndSE<-zYears$se.fit
years<-seq(1999,2013,1)

abundIndLow<-exp(abundInd-1.96*abundIndSE)
abundIndUpp<-exp(abundInd+1.96*abundIndSE)

pollackAbund<-data.frame(cbind(exp(abundInd),abundIndLow,abundIndUpp,years))
colnames(pollackAbund)<-c("Index","CIlow","CIupp","Year")

pdf(file="/Users/jaimeoterovillar/Desktop/NBprobTrend.pdf",width=10,height=8)

ggplot(data=pollackAbund,aes(x=Year,y=Index))+
	geom_segment(aes(x=Year,y=Index,xend=Year,yend=CIupp),col="orange",lwd=0.5)+
	geom_segment(aes(x=Year,y=CIlow,xend=Year,yend=Index),col="orange",lwd=0.5)+
	geom_line(lwd=0.5,linetype="dotted",col="blue")+
	geom_point(size=5,col="gray90")+
	geom_point(size=3,col="orange")+
	scale_y_continuous("Standardized Index",limits=c())+
	scale_x_continuous("Year",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

dev.off()

#12.2.2# Gear

zGears<-predict(negbin1,type="link",se=T,newdata=newGear)

gearP<-zGears$fit
gearSE<-zGears$se.fit
vals<-c(1,2)
probGear<-data.frame(cbind(gearP,gearSE,vals))

#12.2.3# Crew

zCrew<-predict(negbin1,type="link",se=T,newdata=newCrew)

crewP<-zCrew$fit
crewSE<-zCrew$se.fit
probCrew<-data.frame(cbind(crewP,crewSE,vals))

#12.2.4# ZoneO

zZoneO<-predict(negbin1,type="link",se=T,newdata=newZoneO)

zoneP<-zZoneO$fit
zoneSE<-zZoneO$se.fit
vals<-c(1,2,3)
probZone<-data.frame(cbind(zoneP,zoneSE,vals))

#12.2.5# Seafloor

zSeafloor<-predict(negbin1,type="link",se=T,newdata=newSeafloor)

seafloorP<-zSeafloor$fit
seafloorSE<-zSeafloor$se.fit
probSeafloor<-data.frame(cbind(seafloorP,seafloorSE,vals))

#12.2.6# Plots

f7<-ggplot(data=probGear,aes(x=vals,y=exp(gearP)))+
	geom_segment(aes(x=vals,y=exp(gearP-1.96*gearSE),
		xend=vals,yend=exp(gearP+1.96*gearSE)),col="orange",lwd=0.3)+
	geom_point(size=4.5,col="gray90")+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("Standardized Abundance",limits=c(0,0.25))+
	scale_x_continuous("Gear",breaks=c(1,2),labels=c("Veta","MiÃ±o"),
		limits=c(0.75,2.25))
	
f8<-ggplot(data=probCrew,aes(x=vals,y=exp(crewP)))+
	geom_segment(aes(x=vals,y=exp(crewP-1.96*crewSE),
		xend=vals,yend=exp(crewP+1.96*crewSE)),col="orange",lwd=0.3)+
	geom_point(size=4.5,col="gray90")+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("",limits=c(0,0.22))+
	scale_x_continuous("Crew",breaks=c(1,2),labels=c("1-3","4-6"),
		limits=c(0.75,2.25))
		
f9<-ggplot(data=probZone,aes(x=vals,y=exp(zoneP)))+
	geom_segment(aes(x=vals,y=exp(zoneP-1.96*zoneSE),
		xend=vals,yend=exp(zoneP+1.96*zoneSE)),col="orange",lwd=0.3)+
	geom_point(size=4.5,col="gray90")+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("Standardized Abundance",limits=c(0,0.4))+
	scale_x_continuous("Oceanographic Zone",breaks=c(1,2,3),
		labels=c("RÃ­as Baixas","Arco Ãrtabro","CantÃ¡brico"),limits=c(0.75,3.25))

f10<-ggplot(data=probSeafloor,aes(x=vals,y=exp(seafloorP)))+
	geom_segment(aes(x=vals,y=exp(seafloorP-1.96*seafloorSE),
		xend=vals,yend=exp(seafloorP+1.96*seafloorSE)),col="orange",lwd=0.3)+
	geom_point(size=4.5,col="gray90")+
	geom_point(size=2.5,col="orange")+
	scale_y_continuous("",limits=c(0,0.22))+
	scale_x_continuous("Seafloor type",breaks=c(1,2,3),
		labels=c("Hard","Mixed","Soft"),limits=c(0.75,3.25))

pdf(file="/Users/jaimeoterovillar/Desktop/NBprobFacts.pdf")

multiplot(f7,f9,f8,f10,cols=2)

dev.off()


# ----------------------------- #
#13# Considerations about the potential random effects (multiple boats and
# multiple trips within the same boat which performs multiple hauls)

# vv<-length(unique(pollack1$Idflota)) # Number of sampled vessels
# jj<-length(unique(pollack1$Idjornada)) # Number of sampled trips
# hh<-length(unique(pollack1$Idlance)) # Number of hauls 

# pollack1$fvessel<-factor(pollack1$Idflota) # Sampled vessel
# pollack1$fsample<-factor(pollack1$Idjornada) # Sampled day 

#13.1# Identify single observations

# r1<-table(pollack1$fvessel);hist(r1)
# r2<-table(pollack1$fsample);hist(r2)

# length(which(r1==1)) # Number of vessel with a single observation
# length(which(r2==1)) # Number of trips with a single haul

#13.2# Remove vessels with a single observation

# names(which(r1==1))

# pollackRV<-pollack1[-which(pollack1$fvessel=="-1885739447" |
	# pollack1$fvessel=="-1817503041" | pollack1$fvessel=="-1816091332" |
	# pollack1$fvessel=="-1721819171" | pollack1$fvessel=="-739178201" |
	# pollack1$fvessel=="-327161842" | pollack1$fvessel=="1832" |
	# pollack1$fvessel=="2144" | pollack1$fvessel=="2795" |
	# pollack1$fvessel=="2903" | pollack1$fvessel=="3624" |
	# pollack1$fvessel=="4182" | pollack1$fvessel=="4928" |
	# pollack1$fvessel=="5242" | pollack1$fvessel=="5812" |
	# pollack1$fvessel=="6324" | pollack1$fvessel=="6848" |
	# pollack1$fvessel=="7298" | pollack1$fvessel=="7428" |
	# pollack1$fvessel=="8503" | pollack1$fvessel=="8682" |
	# pollack1$fvessel=="8713" | pollack1$fvessel=="979057547" |
	# pollack1$fvessel=="1066635779" | pollack1$fvessel=="1504209601" |
	# pollack1$fvessel=="1743595046"),]


# ----------------------------- #
#14# Modelling of catches using mixed models (DON'T RUN ALL THIS!!!!)

library(lme4)
library(glmmADMB)

#14.1# With glmer

poissM1<-glmer(Ntot~fgear+
					fcrew+lGRT+
					fyear+poly(Julian,2)+
					lDepth+
					QxM+sstM+
					poly(caladoNight,2)+
					offset(log(offs1))+
					(1|fvessel),
					family="poisson",data=pollack1)

#14.2# With ADMB

poissADMB<-glmmadmb(Ntot~fgear+
						 fcrew+lGRT+
						 fyear+poly(Julian,2)+
						 lDepth+
						 QxM+sstM+
						 poly(caladoNight,2)+
						 offset(log(offs1))+
						 (1|fvessel),
						 family="poisson",data=pollack1)

summary(poissADMB)

negbinADMB<-glmmadmb(Ntot~fgear+
						  fcrew+lGRT+
						  fyear+poly(Julian,2)+
						  lDepth+
						  QxM+sstM+
						  poly(caladoNight,2)+
						  offset(log(offs1))+
						  (1|fvessel),
			 	 	      family="nbinom",data=pollack1)

summary(negbinADMB)

poissADMBzi<-glmmadmb(Ntot~fgear+
						   fcrew+lGRT+
						   fyear+poly(Julian,2)+
						   lDepth+
						   QxM+sstM+
						   poly(caladoNight,2)+
						   offset(log(offs1))+
						   (1|fvessel),zeroInflation=T,
			 	 	       family="poisson",data=pollack1)

summary(poissADMBzi)

negbinADMBzi<-glmmadmb(Ntot~fgear+
						    fcrew+lGRT+
						    fyear+poly(Julian,2)+
						    lDepth+
						    QxM+sstM+
						    poly(caladoNight,2)+
						    offset(log(offs1))+
						    (1|fvessel),zeroInflation=T,
			 	 	        family="nbinom",data=pollack1)

summary(negbinADMBzi)

AIC(poissADMB,negbinADMB,poissADMBzi,negbinADMBzi)

