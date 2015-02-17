##########################################################
# Analysis of Trisopterus luscus catch rates in Galicia#
# using data sampled onboard fishing vessels by the UTPB #
##########################################################

#Date start: April 2014 by Jaime Otero (bulk of code from pollachius analysis)
#Last update: 17-02-2015 by Alex Alonso
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
setwd("D:\\iim-csic\\proyectos\\ICES_CASGASS\\analysis\\pouting")

#to start since the last point
load("pouting.RData")

#1# Load species data
#el código está hecho en base al análisis de pollachius
#adaptar el código a cada especies
pouting<-read.table(file="C:\\Users\\alex\\Dropbox\\casgass\\faneca.txt",header=T,dec=".",sep=",")

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

basic.info

pdf(file="gearsPouting.pdf",width=15,height=8)

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
levels(pouting1$Gear)

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
pouting1$Observer<-factor(pouting1$Observer)

#4.2# Gear size

#4.2.1# Recorded gear size in each trip

gearSize<-read.csv2(file="C:\\Users\\alex\\Dropbox\\casgass\\dimensiones_artes_pesca.csv",header=T,dec=".",sep=";")

head(gearSize)

#we are not going to consider trap size as an offset
#only number of traps per haul
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

gear7<-lm(ZoneO~Gear,data=pouting5)
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

write.table(oceano,"D:\\iim-csic\\proyectos\\ICES_CASGASS\\analysis\\oceano.txt",col.names = TRUE,sep=",")

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

#pouting5<-pouting5[-which(is.na(pouting5$QxM)),]
#no la usaremos por lo que no es necesario retirar los NA

names(pouting5)
pouting6<-pouting5[,-c(30:38,40,43)]
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

pouting6$Gear<-relevel(pouting6$Gear,ref="NASA-PEIXES")
poll.dist$Gear<-relevel(poll.dist$Gear,ref="NASA-PEIXES")
summary(poll.dist)

pdf(file="respPouting.pdf",width=15,height=6)

ggplot(data=poll.dist,aes(x=Count,y=Freq))+facet_wrap(~Gear,nrow=3,scales="free_y")+
	geom_bar(stat="identity")+
	scale_y_continuous("Frequency")+
	scale_x_discrete("Number of pouting caught per haul",breaks=seq(0,1915,25))+
	theme(axis.text.x=element_text(size=7),
		  axis.text.y=element_text(size=8))

dev.off()

par(mar=c(5,5,3,3))
hist(pouting6$Ntot,prob=T,ylab="Probability",xlab="Nº of fished pouting",main="")

hist(log(pouting6$Wtot))

plot(log(Ntot)~log(Wtot),data=pouting6,ylab="Total number",xlab="Total biomass")

#6.1.2# Plots of spatial distribution of effort (m2 of net and fishing hour) 

p1<-ggplot(data=pouting6[pouting6$Gear=="MINOS"|pouting6$Gear=="VETAS",],aes(x=Lon,y=Lat))
p1+geom_point(aes(color=Area))+
	scale_color_gradient(low="blue",high="red")
summary(pouting6[pouting6$Gear=="MINOS"|pouting6$Gear=="VETAS",]$Area)

p2<-ggplot(data=pouting6[pouting6$Gear=="NASA-PEIXES",],aes(x=Lon,y=Lat))
p2+geom_point(aes(color=Pieces))+
  scale_color_gradient(low="blue",high="red")
summary(pouting6[pouting6$Gear=="NASA-PEIXES",]$Pieces)

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
abline(lm(GRT~Pieces,data=pouting6),col="red",lwd=2)
plot(GRT~Pieces,data=pouting6)
abline(lm(GRT~Soak,data=pouting6),col="red",lwd=2)
plot(Crew~Area,data=pouting6)
abline(lm(Crew~Area,data=pouting6),col="red",lwd=2)
plot(Crew~Pieces,data=pouting6)
abline(lm(Crew~Pieces,data=pouting6),col="red",lwd=2)
plot(Crew~Soak,data=pouting6)
abline(lm(Crew~Soak,data=pouting6),col="red",lwd=2)
plot(Area~Soak,data=pouting6)
abline(lm(Area~Soak,data=pouting6),col="red",lwd=2)
plot(Pieces~Soak,data=pouting6)
abline(lm(Pieces~Soak,data=pouting6),col="red",lwd=2)

par(mfrow=c(1,3))
boxplot(Area~fcrew,data=pouting6,notch=T,xlab="Crew",ylab="Area")
abline(lm(Area~fcrew,data=pouting6),col="red",lwd=2)
boxplot(Pieces~fcrew,data=pouting6,notch=T,xlab="Crew",ylab="Pieces")
abline(lm(Pieces~fcrew,data=pouting6),col="red",lwd=2)
boxplot(Soak~fcrew,data=pouting6,notch=T,xlab="Crew",ylab="Soak time")
abline(lm(Soak~fcrew,data=pouting6),col="red",lwd=2)

#6.3# Calculation of nighttime fishing

library(lubridate) # Necesaria para obtener horas y minutos de objetos POSIXct

clock.M <- function (t) {hour(t)*60 + minute(t)} # Calcula el minuto entre 1 y 1440 de un dia

pho<-read.csv2(file="C:\\Users\\alex\\Dropbox\\casgass\\Sunrise_set.csv",header=T,dec=".",sep=";")

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

#we are not work with it

#añadimos una nueva variable indicativa del tamaño promedio de los individuos capturados
#Wtot/Ntot puede ser utilizada como indicativo de cambios de distribución ontogénicos???
#summary(pouting6$Ntot)
#summary(pouting6$Wtot)
#plot(pouting6$Ntot,pouting6$Wtot)
#necesitamos eliminar NA de los datos Wtot para poder calcular esta nueva variable
#pouting6$Wtot<-ifelse(pouting6$Ntot==0,0,pouting6$Wtot)
#par(mfrow=c(1,2))
#dotchart(pouting6$Wtot)
#hist(pouting6$Wtot,breaks=100)
#summary(pouting6$Wtot)
#plot(pouting6$Ntot,pouting6$Wtot)
#calculamos la relación
#pouting6$AvgSize<-ifelse(pouting6$Ntot==0,0,pouting6$Wtot/pouting6$Ntot)
#hist(pouting6$AvgSize,breaks=100)
#summary(pouting6$AvgSize)
#does it make sense 0 values for this variable?
#ggplot(data=pouting6,aes(x=AvgSize,y=Ntot))+
#  geom_point()+
#  facet_wrap(~Gear,scales="free")

# ----------------------------- #
#7# Exploratory steps before modelling

#7.1# Distribution of all potential explanatory variables

#likely outlier in enmalle >500
which(pouting6$Depth>500) # Solo un lance por encima de 200 m (utpb dice que es real; pero no representativo de la pesquería)
#pouting7<-pouting6[pouting6$Depth<500,]
pouting7<-pouting6 #lo dejamos como fuente de falso "0"

pouting7$fgear<-relevel(pouting7$fgear,ref="NASA-PEIXES")
levels(pouting7$fgear)
levels(pouting7$Gear)

d1<-ggplot(data=pouting7,aes(x=Crew))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Crew")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))
d2<-ggplot(data=pouting7,aes(x=GRT))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("GRT")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d3<-ggplot(data=pouting7,aes(x=Area))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous(expression(paste("Gear area"," ","(",m^2,")")))+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d4<-ggplot(data=pouting7,aes(x=Soak/1440))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Soak time (days)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d5<-ggplot(data=pouting7,aes(x=hLarg/60))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Deployment (Local time)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d6<-ggplot(data=pouting7,aes(x=hVir/60))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Retrieval (Local time)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d7<-ggplot(data=pouting7,aes(x=Julian))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Day of the Year")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d8<-ggplot(data=pouting7,aes(x=Lat))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Latitude (ºN)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d9<-ggplot(data=pouting7,aes(x=Lon))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Longitude (ºW)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d10<-ggplot(data=pouting7,aes(x=Depth))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Depth (m)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d11<-ggplot(data=pouting7,aes(x=QxM))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous(expression(paste("Qx"," ","(",m^3,s^-1,km^-1,") Ã—",10^3)))+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d12<-ggplot(data=pouting7,aes(x=sstM))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("SST (ºC)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d13<-ggplot(data=pouting7,aes(x=AvgSize))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Individual size (kg)")+
  scale_y_continuous("")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  theme(legend.position="none")
d14<-ggplot(data=pouting7,aes(x=caladoNight))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Night-time soak (%)")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position="none")

pdf(file="explPouting.pdf",width=20,height=8)

multiplot(d1,d7,d8,d2,d5,d9,d3,d6,d10,d4,d14,d12,cols=4)

dev.off()

ggplot(data=pouting7,aes(x=hLargC))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Corrected deployment time")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position=c(0.8,0.8),legend.text=element_text(size=8))

ggplot(data=pouting7,aes(x=hVirC))+
  geom_density(aes(fill=fgear),alpha=0.5)+
  scale_x_continuous("Corrected retrieval time")+
  scale_y_continuous("Density")+
  scale_fill_manual(name="Gear",values=c("orange","#3182bd","light blue"),labels=c("Nasa","Miño","Veta"))+
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
nasa.pouting<-nasa.pouting[,-c(30:32,36:41,44:53,55,57)]
names(nasa.pouting)

#data base for enmalle
enmalle.pouting<-pouting.dat[pouting.dat$Gear=="MINOS"|pouting.dat$Gear=="VETAS",]
dim(enmalle.pouting)
names(enmalle.pouting)
enmalle.pouting<-enmalle.pouting[,-c(30,31,36:41,44:53,55,57)]
names(enmalle.pouting)

#some exploratory plots
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

#7.3# Adding Offset

# nasa

nasa.pouting$OffSet<-nasa.pouting$Pieces*(nasa.pouting$Soak/60)
summary(nasa.pouting$OffSet);hist(log(nasa.pouting$OffSet),xlab="OffSet",main="")

#enmalle

enmalle.pouting$OffSet<-enmalle.pouting$Area*(enmalle.pouting$Soak/60)
summary(enmalle.pouting$OffSet);hist(log(enmalle.pouting$OffSet),xlab="OffSet",main="")

#7.4# Collinearity

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
vifs1<-c("lGRT","Crew","Year","Julian","Depth","sstM","caladoNight")
vifs2<-c("lGRT","Crew","Year","Julian","lDepth","sstM","caladoNight")
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
nasa.pouting$Observer<-as.factor(nasa.pouting$Observer)
nasa.pouting$Idflota<-as.factor(nasa.pouting$Idflota)
nasa.pouting$Idlance<-as.factor(nasa.pouting$Idlance)
nasa.pouting$Idjornada<-as.factor(nasa.pouting$Idjornada)
#not well ditributed GRT
#so we are going to do 3 main classes to avoid extrange patterns in linear response
# 0-8 / 8-12 / >12
nasa.pouting$fGRT<-as.factor(ifelse(nasa.pouting$GRT<8,"C8",
                                    ifelse(nasa.pouting$GRT<12,"C8-12","C12")))
nasa.pouting$fGRT<-factor(nasa.pouting$fGRT,levels = c("C8", "C8-12", "C12"))
str(nasa.pouting$fGRT)
boxplot(Ntot~fGRT,data=nasa.pouting)

names(enmalle.pouting)
enmalle.pouting$fyear<-as.factor(enmalle.pouting$fyear)
enmalle.pouting$fcrew<-as.factor(enmalle.pouting$fcrew)
enmalle.pouting$fZoneA<-as.factor(enmalle.pouting$fZoneA)
enmalle.pouting$fZoneO<-as.factor(enmalle.pouting$fZoneO)
enmalle.pouting$Observer<-as.factor(enmalle.pouting$Observer)
enmalle.pouting$Idflota<-as.factor(enmalle.pouting$Idflota)
enmalle.pouting$Idlance<-as.factor(enmalle.pouting$Idlance)
enmalle.pouting$Idjornada<-as.factor(enmalle.pouting$Idjornada)
#not well ditributed GRT
#so we are going to do 3 main classes to avoid extrange patterns in linear response
# 0-8 / 8-12 / >12
enmalle.pouting$fGRT<-as.factor(ifelse(enmalle.pouting$GRT<5,"C5",
                                    ifelse(enmalle.pouting$GRT<10,"C5-10",
                                           ifelse(enmalle.pouting$GRT<15,"C10-15",
                                                  ifelse(enmalle.pouting$GRT<20,"C15-20","C>20")))))
enmalle.pouting$fGRT<-factor(enmalle.pouting$fGRT,levels = c("C5", "C5-10", "C10-15","C15-20","C>20"))
str(enmalle.pouting$fGRT)
boxplot(Ntot~fGRT,data=enmalle.pouting)

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

nasa.poiss0<-gam(Ntot~offset(log(OffSet))+
                   fGRT+
                   fyear+s(Julian,k=6,bs="cc")+
                   s(Depth,k=3)+
                   s(sstM,k=3)+
                   s(caladoNight,k=3),
                 family=poisson,data=nasa.pouting)

summary(nasa.poiss0)
anova(nasa.poiss0)
plot(nasa.poiss0,pages=1,all.terms=T,scale=0,shade="T")
sum(residuals(nasa.poiss0,type="pearson")^2)/nasa.poiss0$df.resid # Overdispersion

#8.1.2# GAM NegBin modelling

nasa.negbin0<-gam(Ntot~offset(log(OffSet))+
                    fGRT+
                    fyear+s(Julian,k=6,bs="cc")+
                    s(Depth,k=3)+
                    s(sstM,k=3)+
                    s(caladoNight,k=3),
                  family=nb(),data=nasa.pouting)

summary(nasa.negbin0)
anova(nasa.negbin0)
plot(nasa.negbin0,pages=1,all.terms=T,scale=0,shade="T")
sum(residuals(nasa.negbin0,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin0))+1)) 

#8.1.3# GLM Poisson modelling

nasa.poiss1<-glm(Ntot~offset(log(OffSet))+fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight,
			 	 family=poisson,data=nasa.pouting)

summary(nasa.poiss1)
Anova(nasa.poiss1)
plot(allEffects(nasa.poiss1))
sum(residuals(nasa.poiss1,type="pearson")^2)/nasa.poiss1$df.resid # Overdispersion

#8.1.4# Negative Binomial modelling

#indirect parametrization GLM.NB
nasa.negbin1<-glm.nb(Ntot~offset(log(OffSet))+fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight,
                     data=nasa.pouting)

summary(nasa.negbin1)
Anova(nasa.negbin1)
plot(allEffects(nasa.negbin1))
sum(residuals(nasa.negbin1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin1))+1))

#direct parametrization (Joseph M. HILBE)
#http://cran.r-project.org/web/packages/msme/msme.pdf
#http://cran.r-project.org/web/packages/COUNT/COUNT.pdf
#https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=10&cad=rja&uact=8&ved=0CG4QFjAJ&url=http%3A%2F%2Fwww.cambridge.org%2Fch%2Fdownload_file%2F846106%2F&ei=xofbVOLSA4WBU82HgbAC&usg=AFQjCNEGGL69wK3m3m5Q7K8tN5ZHljTFFQ&sig2=42OpyYD42s7ju9B2OPcdNw

OFFSET<-nasa.pouting$OffSet

# TRADITIONAL NB REGRESSION WITH ALPHA
# do not pay attention to warnings()

nasa.negbin2<-nbinomial(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight,
                        formula2 = ~ 1,
                        data=nasa.pouting,offset=log(OFFSET),family = "nb2",mean.link = "log",scale.link = "inverse_s")
#the same as
#nasa.negbin2<-nbinomial(Ntot~offset(log(OffSet))+fGRT+fyear+poly(Julian,2)+poly(Depth,2)+sstM+caladoNight,
#                        formula2 = ~ 1,
#                        family = "nb2",mean.link = "log",scale.link = "inverse_s",data=nasa.pouting)

summary(nasa.negbin2)
sum(residuals(nasa.negbin2,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin2))+1))
#dispersion statistics diferent to calculated with pearson residuals (teh one provided in the summary)

#TRADITIONAL NB REGRESSION with DISPERSION PARAMETERIZED

nasa.negbin3<-nbinomial(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight,
                        formula2 = ~ fyear+Depth,
                        family = "nb2",mean.link = "log",scale.link = "inverse_s",offset=log(OFFSET),data=nasa.pouting)

summary(nasa.negbin3)
sum(residuals(nasa.negbin3,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin3))+1))

# R GLM.NB-TYPE INVERTED DISPERSON --THETA ; WITH DEFAULTS

nasa.negbin4 <-nbinomial(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight,
                         formula2 = ~ 1,
                         data = nasa.pouting,offset=log(OFFSET),family = "negBinomial",mean.link = "log",scale.link = "log_s")

summary(nasa.negbin4)
sum(residuals(nasa.negbin4,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin4))+1))

# HETEROGENEOUS NB; DISPERSION PARAMETERIZED

nasa.negbin5 <- nbinomial(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight,
                   formula2 = ~ fyear+Depth,
                   data = nasa.pouting,offset=log(OFFSET),family = "negBinomial",mean.link = "log",scale.link = "log_s")

summary(nasa.negbin5)
sum(residuals(nasa.negbin5,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin5))+1))

#GLMMADMB

nasa.pouting$Nasaoffset<-log(nasa.pouting$OffSet)

#negative binomial (log): "nbinom", "nbinom1" (negative binomial "type 1": see Details), "nbinom2" (a synonym for "nbinom")

#nbinom
nasa.negbin6 <- glmmadmb(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight+offset(log(OffSet)),
                         data=nasa.pouting,zeroInflation=FALSE,family="nbinom")

summary(nasa.negbin6)
Anova(nasa.negbin6)
sum(residuals(nasa.negbin6,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin6))+1)) 

#type 1
nasa.negbin7 <- glmmadmb(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight+offset(log(OffSet)),
                         data=nasa.pouting, zeroInflation=FALSE,family="nbinom1")

summary(nasa.negbin7)
Anova(nasa.negbin7)
sum(residuals(nasa.negbin7,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin7))+1)) 

#comparing NG parametrizations

coef(nasa.negbin1)
coef(nasa.negbin2)
coef(nasa.negbin3)
coef(nasa.negbin4)
coef(nasa.negbin5)
coef(nasa.negbin6)
coef(nasa.negbin7)

modelfit(nasa.negbin1)
modelfit(nasa.negbin2)
modelfit(nasa.negbin3)
modelfit(nasa.negbin4)
modelfit(nasa.negbin5)

AICNB<-data.frame(cbind(c("nb1","nb2","nb3","nb4","nb5","nb6","nb7"),
                        as.numeric(c(AIC(nasa.negbin1),nasa.negbin2$aic,nasa.negbin3$aic,nasa.negbin4$aic,nasa.negbin5$aic,AIC(nasa.negbin6),AIC(nasa.negbin7)))))
DispNB<-data.frame(c(sum(residuals(nasa.negbin1,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin1))+1)),
                     sum(residuals(nasa.negbin2,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin2))+1)),
                     sum(residuals(nasa.negbin3,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin3))+1)),
                     sum(residuals(nasa.negbin4,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin4))+1)),
                     sum(residuals(nasa.negbin5,type="deviance")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin5))+1)),
                     sum(residuals(nasa.negbin6,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin6))+1)),
                     sum(residuals(nasa.negbin7,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin7))+1))))

compMNB<-data.frame(cbind(AICNB,DispNB))
colnames(compMNB)<-c("Model","AIC","Dispersion")
compMNB$AIC<-round(as.numeric(levels(compMNB$AIC))[compMNB$AIC],digits = 2)
compMNB

nb1<-data.frame(list(summary(nasa.negbin1)$coefficients)[[1]])
nb1<-data.frame(cbind(nb1,confint(nasa.negbin1)))
nb1$Coef<-rownames(nb1)
colnames(nb1)<-c("Estimate","SE","Z","p","LCL","UCL","Coef")
nb1$Mod<-"glm.nb"
dim(nb1)
nb2<-list(summary(nasa.negbin2)$coefficients)[[1]]
nb2$Coef<-rownames(nb2)
nb2<-nb2[c(1:14),]
nb2$Mod<-"nb2"
dim(nb2)
nb3<-list(summary(nasa.negbin3)$coefficients)[[1]]
nb3$Coef<-rownames(nb3)
nb3$Mod<-"nb2-Mdisp"
nb3<-nb3[c(1:14),]
dim(nb3)
nb4<-list(summary(nasa.negbin4)$coefficients)[[1]]
nb4$Coef<-rownames(nb4)
nb4$Mod<-"negBinomial"
nb4<-nb4[c(1:14),]
dim(nb4)
nb5<-list(summary(nasa.negbin5)$coefficients)[[1]]
nb5$Coef<-rownames(nb5)
nb5$Mod<-"negBinomial-Mdisp"
nb5<-nb5[c(1:14),]
dim(nb5)
nb6<-data.frame(list(summary(nasa.negbin6)$coefficients)[[1]])
nb6<-data.frame(cbind(nb6,confint(nasa.negbin6)))
nb6$Coef<-rownames(nb6)
colnames(nb6)<-c("Estimate","SE","Z","p","LCL","UCL","Coef")
nb6$Mod<-"ADMB-nbinom"
dim(nb6)
nb7<-data.frame(list(summary(nasa.negbin7)$coefficients)[[1]])
nb7<-data.frame(cbind(nb7,confint(nasa.negbin7)))
nb7$Coef<-rownames(nb7)
colnames(nb7)<-c("Estimate","SE","Z","p","LCL","UCL","Coef")
nb7$Mod<-"ADMB-type1"
dim(nb7)

NbComp<-data.frame(rbind(nb1,nb2,nb3,nb4,nb5,nb6,nb7))
str(NbComp)

pdf(file="NBpouting-coef.pdf",width=15,height=8)

limits <- aes(ymax =  UCL, ymin= LCL)
ggplot(data=NbComp,aes(x=Mod,y=Estimate))+
  geom_errorbar(limits, colour="blue", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="blue")+
  geom_hline(aes(yintercept=0),col="red")+
  scale_y_continuous("Estimate")+
  scale_x_discrete("NB Model")+
  facet_wrap(~Coef,scales="free_y", ncol = 5)+
  ggtitle("NegBin Parametrization")+
  theme(plot.title = element_text(size=26, face="bold"),axis.text.x=element_text(angle=-45, hjust = 0))

dev.off()

gdnb1<-data.frame(cbind(residuals(nasa.negbin1,type="deviance"),fitted(nasa.negbin1),nasa.pouting$Idflota))
colnames(gdnb1)<-c("Residuals","Fitted","Boat")
gdnb1$Mod<-"glm.nb"
gdnb1$Boat<-factor(as.numeric(gdnb1$Boat))
dim(gdnb1)

gdnb2<-data.frame(cbind(residuals(nasa.negbin2,type="deviance"),fitted(nasa.negbin2),nasa.pouting$Idflota))
colnames(gdnb2)<-c("Residuals","Fitted","Boat")
gdnb2$Mod<-"nb2"
gdnb2$Boat<-factor(as.numeric(gdnb2$Boat))
dim(gdnb2)

gdnb3<-data.frame(cbind(residuals(nasa.negbin3,type="deviance"),fitted(nasa.negbin3),nasa.pouting$Idflota))
colnames(gdnb3)<-c("Residuals","Fitted","Boat")
gdnb3$Mod<-"nb2-Mdisp"
gdnb3$Boat<-factor(as.numeric(gdnb3$Boat))
dim(gdnb3)

gdnb4<-data.frame(cbind(residuals(nasa.negbin4,type="deviance"),fitted(nasa.negbin4),nasa.pouting$Idflota))
colnames(gdnb4)<-c("Residuals","Fitted","Boat")
gdnb4$Mod<-"negBinomial"
gdnb4$Boat<-factor(as.numeric(gdnb4$Boat))
dim(gdnb4)

gdnb5<-data.frame(cbind(residuals(nasa.negbin5,type="deviance"),fitted(nasa.negbin5),nasa.pouting$Idflota))
colnames(gdnb5)<-c("Residuals","Fitted","Boat")
gdnb5$Mod<-"negBinomial-Mdisp"
gdnb5$Boat<-factor(as.numeric(gdnb5$Boat))
dim(gdnb5)

gdnb6<-data.frame(cbind(residuals(nasa.negbin6,type="pearson"),fitted(nasa.negbin6),nasa.pouting$Idflota))
colnames(gdnb6)<-c("Residuals","Fitted","Boat")
gdnb6$Mod<-"ADMB-nbinom"
gdnb6$Boat<-factor(as.numeric(gdnb6$Boat))
dim(gdnb6)

gdnb7<-data.frame(cbind(residuals(nasa.negbin7,type="pearson"),fitted(nasa.negbin7),nasa.pouting$Idflota))
colnames(gdnb7)<-c("Residuals","Fitted","Boat")
gdnb7$Mod<-"ADMB-type1"
gdnb7$Boat<-factor(as.numeric(gdnb7$Boat))
dim(gdnb7)

gdnb<-data.frame(rbind(gdnb1,gdnb2,gdnb3,gdnb4,gdnb5,gdnb6,gdnb7))
str(gdnb)
gdnb$Lfit<-log(gdnb$Fitted)

pdf(file="NBpouting-histres.pdf",width=12,height=8)

ggplot(gdnb, aes(x=Residuals))+
  geom_histogram(binwidth = 0.2,aes(fill = ..count..))+
  scale_fill_gradient("Count", low = "green", high = "red") +
  geom_vline(aes(xintercept=0),col="black")+
  facet_wrap(~Mod,scales="free", ncol = 4)

dev.off()

ggplot(gdnb, aes(y=Residuals,x=Fitted,colour=Residuals)) + 
  geom_point(size=3)+
  scale_colour_gradient(limits=c(-5, 5),low="green", high="red")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Fitted")+
  facet_wrap(~Mod,scales="free", ncol = 4)

#fitted log-transformed

pdf(file="NBpouting-fitres.pdf",width=12,height=8)

ggplot(gdnb, aes(y=Residuals,x=Lfit,colour=Residuals)) + 
  geom_point(size=3)+
  scale_colour_gradient(limits=c(-4.5, 4.5),low="green", high="red")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("log(Fitted)")+
  facet_wrap(~Mod,scales="free", ncol = 4)

dev.off()

pdf(file="NBpouting-resmixed.pdf",width=12,height=8)

ggplot(gdnb, aes(x=factor(Boat), y=Residuals))+
  geom_boxplot(col="gray50",fill="gray75")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Boat")+
  facet_wrap(~Mod,scales="free", ncol = 4)

dev.off()

#best soluction
#nasa.negbin3 has the best AIC but we prefer the simpler solution based on dispersion statistics

nasa.negbin1
summary(nasa.negbin1)

#8.1.5 Mixed model NB
#http://glmmadmb.r-forge.r-project.org/glmmADMB.html
#https://rpubs.com/bbolker/glmmchapter

nasa.negbinm1 <- glmmadmb(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight+offset(log(OffSet))+(1|Idflota), 
                            data=nasa.pouting,zeroInflation=FALSE,family="nbinom")
summary(nasa.negbinm1)
Anova(nasa.negbinm1)
sum(residuals(nasa.negbinm1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbinm1))+1))
#no mejoramos la sobredispersión

par(mfrow=c(1,3))
hist(residuals(nasa.negbinm1,type="pearson"))
plot(residuals(nasa.negbinm1,type="pearson")~log(fitted(nasa.negbinm1)))
boxplot(residuals(nasa.negbinm1,type="pearson")~nasa.pouting$Idflota)
abline(0,0,col="red")

#8.1.6# ZI
#exploring binomial 

nasa.pouting$binNtot<-ifelse(nasa.pouting$Ntot==0,0,1)

nasa.bin<-glm(binNtot~lGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight,
                 family=binomial,data=nasa.pouting)
summary(nasa.bin)
anova(nasa.bin)
plot(allEffects(nasa.bin))
#no clear patterns

#8.1.6.1# ZIpoisson
nasa.zipoiss1<-zeroinfl(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight+offset(log(OffSet))|1,
                           dist="poisson",link="logit",data=nasa.pouting)

summary(nasa.zipoiss1)
sum(residuals(nasa.zipoiss1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zipoiss1))))

#8.1.6.2# ZINB
nasa.zinb1<-zeroinfl(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight+offset(log(OffSet))|1,
                        dist="negbin",link="logit",data=nasa.pouting)

summary(nasa.zinb1)
sum(residuals(nasa.zinb1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zinb1))+1))

#glmmadmb does not allow to include explanmatory variables in ZI part
nasa.zinbm1 <- glmmadmb(Ntot~fGRT+fyear+poly(Julian,2)+Depth+sstM+caladoNight+offset(log(OffSet))+(1|Idflota), 
                             data=nasa.pouting,zeroInflation=TRUE,family="nbinom")

summary(nasa.zinbm1)
sum(residuals(nasa.zinbm1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zinbm1))+1))
#even worse

#8.1.7# full model comparison
# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

#8.1.7.1# AIC/BIC

AICtab(nasa.poiss0,nasa.poiss1,nasa.negbin0,nasa.negbin1,nasa.negbinm1,nasa.zipoiss1,nasa.zinb1,nasa.zinbm1)

nasa.aic=AIC(nasa.poiss0,nasa.poiss1,nasa.negbin0,nasa.negbin1,nasa.negbinm1,nasa.zipoiss1,nasa.zinb1,nasa.zinbm1,k=2)
nasa.aic[order(nasa.bic$AIC), ]
nasa.bic=AIC(nasa.poiss0,nasa.poiss1,nasa.negbin0,nasa.negbin1,nasa.negbinm1,nasa.zipoiss1,nasa.zinb1,nasa.zinbm1,k=log(dim(nasa.pouting)[1]))
nasa.bic[order(nasa.bic$AIC), ]

#8.1.7.2# overdispersion parameter
#If the dispersion parameter is significantly greater than one, indicating overdispersion (variance greater than the mean),
#then the scale parameter should be used to adjust the variance.  
#Failing to account for the overdispersion can result in inflated test statistics. 
#However, when the dispersion parameter is less than one, then the test statistics become more conservative,
#which is not considered as much of a problem.

nasa.disp<-data.frame(cbind(c("nasa.poiss0","nasa.poiss1","nasa.negbin0","nasa.negbin1","nasa.negbinm1","nasa.zipoiss1","nasa.zinb1","nasa.zinbm1"),
                      c(sum(residuals(nasa.poiss0,type="pearson")^2)/nasa.poiss0$df.resid,
                        sum(residuals(nasa.poiss1,type="pearson")^2)/nasa.poiss1$df.resid,
                        sum(residuals(nasa.negbin0,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin0))+1)),
                        sum(residuals(nasa.negbin1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbin1))+1)),
                        sum(residuals(nasa.negbinm1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.negbinm1))+1)),
                        sum(residuals(nasa.zipoiss1,type="pearson")^2)/nasa.zipoiss1$df.resid,
                        sum(residuals(nasa.zinb1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zinb1))+1)),
                        sum(residuals(nasa.zinbm1,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.zinbm1))+1)))))
colnames(nasa.disp)<-c("Model","Dispersion")
nasa.disp$Dispersion<-round(as.numeric(levels(nasa.disp$Dispersion))[nasa.disp$Dispersion],digits = 3)
nasa.disp[order(nasa.disp$Dispersion),]
#sobredispersión un tanto elevada

#8.1.7.3# loglikelihood

fm.nasa<-list("nasa.poiss0"=nasa.poiss0,
               "nasa.poiss1"=nasa.poiss1,
               "nasa.negbin0"=nasa.negbin0,
               "nasa.negbin1"=nasa.negbin1,
               "nasa.negbinm1"=nasa.negbinm1,
               "nasa.zipoiss1"=nasa.zipoiss1,
               "nasa.zinb1"=nasa.zinb1,
               "nasa.zinbm1"=nasa.zinbm1)

logliks.nasa<-data.frame(cbind(c("nasa.poiss0","nasa.poiss1","nasa.negbin0","nasa.negbin1","nasa.negbinm1","nasa.zipoiss1","nasa.zinb1","nasa.zinbm1"),
                    logLik=sapply(fm.nasa,function (x) round(logLik(x),digits=0)),
                     Df=sapply(fm.nasa,function (x) attr(logLik(x),"df"))))
colnames(logliks.nasa)<-c("Model","logLik","DF")
logliks.nasa$logLik<-round(as.numeric(levels(logliks.nasa$logLik))[logliks.nasa$logLik],digits = 0)
logliks.nasa$DF<-round(as.numeric(levels(logliks.nasa$DF))[logliks.nasa$DF],digits = 2)
logliks.nasa[order(-logliks.nasa$logLik),]

#8.1.7.4# Vuong test to compare NB vs ZINB
# do not run with random models

vuong(nasa.zinb1,nasa.negbin1)
#mejor NB

#8.1.7.5# Predicting probabilities

#no mixed models neither gams included in this comparison
phat.pois.nasa<-predprob(nasa.poiss1)
phat.pois.mn.nasa<-apply(phat.pois.nasa,2,mean)
phat.nb.nasa<-predprob(nasa.negbin1)
phat.nb.mn.nasa<-apply(phat.nb.nasa,2,mean)
phat.zipoiss.nasa<-predprob(nasa.zipoiss1)
phat.zipoiss.mn.nasa<-apply(phat.zipoiss.nasa,2,mean)
phat.zinb.nasa<-predprob(nasa.zinb1)
phat.zinb.mn.nasa<-apply(phat.zinb.nasa,2,mean)

with(nasa.pouting,{
  hist(Ntot,prob=TRUE,col="grey",breaks=seq(-0.5,246,2),xlab="",main="")
  lines(x=seq(0,245,1),y=phat.zipoiss.mn.nasa,type="l",lwd=2,col="green")
  lines(x=seq(0,245,1),y=phat.zinb.mn.nasa,type="l",lwd=2,col="yellow")
  lines(x=seq(0,245,1),y=phat.pois.mn.nasa,type="l",lwd=2,col="red")
  lines(x=seq(0,245,1),y=phat.nb.mn.nasa,type="l",lwd=2,col="blue")
})

#8.1.7.6# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

nasa.aic
nasa.bic
logliks.nasa
logliks.nasa2<-rbind(logLik=sapply(fm.nasa,function (x) round(logLik(x),digits=0)),
                     Df=sapply(fm.nasa,function (x) attr(logLik(x),"df")))

zeroes.nasa<-c("Obs"=sum(nasa.pouting$Ntot<1),
                     "nasa.poiss0"=sum(dpois(0,fitted(nasa.poiss0))),
                     "nasa.poiss1"=sum(dpois(0,fitted(nasa.poiss1))),
                     "nasa.negbin0"="NA",
                     "nasa.negbin1"=sum(dnbinom(0,mu=fitted(nasa.negbin1),size=nasa.negbin1$theta)),
                     "nasa.negbinm1"="NA",
                     "nasa.zipoiss1"=sum(predict(nasa.zipoiss1,type="prob")[,1]),
                     "nasa.zinb1"=sum(predict(nasa.zinb1,type="prob")[,1]),
                     "nasa.zinbm1"="NA")

modelsComp.nasa<-data.frame(cbind(logliks.nasa2[c(2,4,6,8,10,12,14,16)],nasa.bic[[2]],nasa.aic[[2]],logliks.nasa2[c(1,3,5,7,9,11,13,15)],zeroes.nasa[2:9]))
colnames(modelsComp.nasa)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes.nasa[[1]],")"))
modelsComp.nasa$Df<-round(as.numeric(levels(modelsComp.nasa$Df))[modelsComp.nasa$Df],digits = 2)
modelsComp.nasa$BIC<-round(as.numeric(levels(modelsComp.nasa$BIC))[modelsComp.nasa$BIC],digits = 2)
modelsComp.nasa$AIC<-round(as.numeric(levels(modelsComp.nasa$AIC))[modelsComp.nasa$AIC],digits = 2)
modelsComp.nasa$logLik<-round(as.numeric(levels(modelsComp.nasa$logLik))[modelsComp.nasa$logLik],digits = 0)
modelsComp.nasa$Zeroes<-round(as.numeric(levels(modelsComp.nasa$Zeroes))[modelsComp.nasa$Zeroes],digits = 4)
modelsComp.nasa$deltaAIC<-round(modelsComp.nasa$AIC-min(modelsComp.nasa$AIC),2)
modelsComp.nasa$deltaBIC<-round(modelsComp.nasa$BIC-min(modelsComp.nasa$BIC),2)
modelsComp.nasa$ModLikelihood<-round(exp(-modelsComp.nasa$deltaAIC/2),2)
modelsComp.nasa$AICweight<-round(modelsComp.nasa$ModLikelihood/sum(modelsComp.nasa$ModLikelihood),2)
modelsComp.nasa<-modelsComp.nasa[order(modelsComp.nasa$deltaAIC),]
modelsComp.nasa

#8.2#ENMALLE modelling

#8.2.1# GAM Poisson modelling

names(enmalle.pouting)

enmalle.poiss0<-gam(Ntot~offset(log(OffSet))+
                      fGRT*Gear+
                      fyear*fZoneO+s(Julian,k=6,bs="cc")+
                      s(lDepth,by=Gear,k=3)+
                      s(sstM,k=3)+
                      s(caladoNight,by=Gear,k=3)+
                      Seafloor,
                    family=poisson,data=enmalle.pouting)

summary(enmalle.poiss0)
anova(enmalle.poiss0)
plot(enmalle.poiss0,pages=1,all.terms=T,scale=0,shade="T")
sum(residuals(enmalle.poiss0,type="pearson")^2)/enmalle.poiss0$df.resid # Overdispersion

#8.2.2# GAM NegBin modelling
enmalle.negbin0<-gam(Ntot~offset(log(OffSet))+
                      fGRT*Gear+
                      fyear*fZoneO+
                      s(Julian,k=6,bs="cc")+
                      s(lDepth,by=Gear,k=3)+
                      s(sstM,k=3)+
                      s(caladoNight,by=Gear,k=3)+
                      Seafloor,
                     family=nb(),data=enmalle.pouting)

summary(enmalle.negbin0)
anova(enmalle.negbin0)
plot(enmalle.negbin0,pages=1,all.terms=T,scale=0,shade="T")
sum(residuals(enmalle.negbin0,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbin0))+1)) # Overdispersion

#8.2.3# GLM Poisson modelling

enmalle.poiss1<-glm(Ntot~offset(log(OffSet))+
                      Gear*fGRT+
                      fyear*fZoneO+
                      poly(Julian,2)+
                      lDepth*Gear+
                      sstM+
                      caladoNight*Gear+
                      Seafloor,
                 family=poisson,data=enmalle.pouting)

summary(enmalle.poiss1)
Anova(enmalle.poiss1)
plot(allEffects(enmalle.poiss1))
sum(residuals(enmalle.poiss1,type="pearson")^2)/enmalle.poiss1$df.resid # Overdispersion

#8.2.4# Negative Binomial modelling

enmalle.negbin1<-glm.nb(Ntot~offset(log(OffSet))+
                          Gear*fGRT+
                          fyear*fZoneO+
                          poly(Julian,2)+
                          lDepth*Gear+
                          sstM+
                          caladoNight*Gear+
                          Seafloor,
                        data=enmalle.pouting)

summary(enmalle.negbin1)
Anova(enmalle.negbin1)
plot(allEffects(enmalle.negbin1))
sum(residuals(enmalle.negbin1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbin1))+1)) 

#8.2.5# Mixed model NB

enmalle.negbinm1 <- glmmadmb(Ntot~Gear*fGRT+
                                 fyear*fZoneO+
                                 poly(Julian,2)+
                                 lDepth*Gear+
                                 sstM+
                                 caladoNight*Gear+
                                 Seafloor+
                                 offset(log(OffSet))+
                                 (1|Idflota),
                               data=enmalle.pouting,
                               zeroInflation=FALSE,
                               family="nbinom")
summary(enmalle.negbinm1)
Anova(enmalle.negbinm1)
sum(residuals(enmalle.negbinm1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbinm1))+1)) 

par(mfrow=c(1,3))
hist(residuals(enmalle.negbinm1,type="pearson"))
plot(residuals(enmalle.negbinm1,type="pearson")~log(fitted(enmalle.negbinm1)))
boxplot(residuals(enmalle.negbinm1,type="pearson")~enmalle.pouting$Idflota)
abline(0,0,col="red")

#8.2.6# ZI
#exploring binomial 
enmalle.pouting$binNtot<-ifelse(enmalle.pouting$Ntot==0,0,1)

enmalle.bin<-glm(binNtot~Gear+fGRT+fyear+fZoneO+poly(Julian,2)+lDepth+sstM+caladoNight+Seafloor,
                 family=binomial,data=enmalle.pouting)
summary(enmalle.bin)
Anova(enmalle.bin)
plot(allEffects(enmalle.bin))

# ZIPoiss
#zeros intercept
enmalle.zipoiss1<-zeroinfl(Ntot~offset(log(OffSet))+
                             Gear*fGRT+
                             fyear*fZoneO+
                             poly(Julian,2)+
                             poly(lDepth,2)*Gear+
                             sstM+
                             caladoNight*Gear+
                             Seafloor|1,
                           dist="poisson",link="logit",data=enmalle.pouting)

summary(enmalle.zipoiss1)
sum(residuals(enmalle.zipoiss1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zipoiss1))))

#zeros modeled
enmalle.zipoiss2<-zeroinfl(Ntot~offset(log(OffSet))+
                             Gear*fGRT+
                             fyear*fZoneO+
                             poly(Julian,2)+
                             lDepth*Gear+
                             sstM+
                             caladoNight*Gear+
                             Seafloor|lDepth,
                           dist="poisson",link="logit",data=enmalle.pouting)

summary(enmalle.zipoiss2)
sum(residuals(enmalle.zipoiss2,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zipoiss2))))

#ZINB
#zeros intercept
enmalle.zinb1<-zeroinfl(Ntot~offset(log(OffSet))+
                          Gear*fGRT+
                          fyear*fZoneO+
                          poly(Julian,2)+
                          lDepth*Gear+
                          sstM+
                          caladoNight*Gear+
                          Seafloor|1,
                        dist="negbin",link="logit",data=enmalle.pouting)

summary(enmalle.zinb1)
sum(residuals(enmalle.zinb1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zinb1))+1))

#zeros modeled
enmalle.zinb2<-zeroinfl(Ntot~offset(log(OffSet))+
                          Gear*fGRT+
                          fyear*fZoneO+
                          poly(Julian,2)+
                          lDepth*Gear+
                          sstM+
                          caladoNight*Gear+
                          Seafloor|lDepth,
                        dist="negbin",link="logit",data=enmalle.pouting)

summary(enmalle.zinb2)
sum(residuals(enmalle.zinb2,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zinb2))+1))

#Mixed ZINB
#we decided not to include random effects
#glmmadmb does not allow to include explanmatory variables in ZI part
#enmalle.zinbm1 <- glmmadmb(Ntot~Gear*lGRT+
#                               fyear*fZoneO+
#                               poly(Julian,2)+
#                               lDepth*Gear+
#                               sstM+
#                               Gear*caladoNight+
#                               Seafloor+
#                               offset(log(OffSet))+
#                               (1|Idflota), 
#                             data=enmalle.pouting,
#                             zeroInflation=TRUE, 
#                             family="nbinom")

#Parameters were estimated, but not standard errors were not

#8.2.7# full model comparison
# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

#8.2.7.1# AIC/BIC

AICtab(enmalle.poiss0,enmalle.poiss1,
       enmalle.negbin0,enmalle.negbin1,enmalle.negbinm1,
       enmalle.zipoiss1,enmalle.zipoiss2,enmalle.zinb1,enmalle.zinb2)

enmalle.aic=AIC(enmalle.poiss0,enmalle.poiss1,
             enmalle.negbin0,enmalle.negbin1,enmalle.negbinm1,
             enmalle.zipoiss1,enmalle.zipoiss2,enmalle.zinb1,enmalle.zinb2,k=2)
enmalle.aic[order(enmalle.aic$AIC), ]
enmalle.bic=AIC(enmalle.poiss0,enmalle.poiss1,
             enmalle.negbin0,enmalle.negbin1,enmalle.negbinm1,
             enmalle.zipoiss1,enmalle.zipoiss2,enmalle.zinb1,enmalle.zinb2,k=log(dim(enmalle.pouting)[1]))
enmalle.bic[order(enmalle.bic$AIC), ]

#8.2.7.2# overdispersion parameter

enmalle.disp<-data.frame(cbind(c("enmalle.poiss0","enmalle.poiss1",
                              "enmalle.negbin0","enmalle.negbin1","enmalle.negbinm1",
                              "enmalle.zipoiss1","enmalle.zipoiss2","enmalle.zinb1","enmalle.zinb2"),
                            c(sum(residuals(enmalle.poiss0,type="pearson")^2)/enmalle.poiss0$df.resid,
                              sum(residuals(enmalle.poiss1,type="pearson")^2)/enmalle.poiss1$df.resid,
                              sum(residuals(enmalle.negbin0,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbin0))+1)),
                              sum(residuals(enmalle.negbin1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbin1))+1)),
                              sum(residuals(enmalle.negbinm1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.negbinm1))+1)),
                              sum(residuals(enmalle.zipoiss1,type="pearson")^2)/enmalle.zipoiss1$df.resid,
                              sum(residuals(enmalle.zipoiss2,type="pearson")^2)/enmalle.zipoiss2$df.resid,
                              sum(residuals(enmalle.zinb1,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zinb1))+1)),
                              sum(residuals(enmalle.zinb2,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.zinb2))+1)))))
colnames(enmalle.disp)<-c("Model","Dispersion")
enmalle.disp$Dispersion<-round(as.numeric(levels(enmalle.disp$Dispersion))[enmalle.disp$Dispersion],digits = 3)
enmalle.disp[order(enmalle.disp$Dispersion),]

#8.2.7.3# loglikelihood

fm.enmalle<-list("enmalle.poiss0"=enmalle.poiss0,
              "enmalle.poiss1"=enmalle.poiss1,
              "enmalle.negbin0"=enmalle.negbin0,
              "enmalle.negbin1"=enmalle.negbin1,
              "enmalle.negbinm1"=enmalle.negbinm1,
              "enmalle.zipoiss1"=enmalle.zipoiss1,
              "enmalle.zipoiss2"=enmalle.zipoiss2,
              "enmalle.zinb1"=enmalle.zinb1,
              "enmalle.zinb2"=enmalle.zinb2)

logliks.enmalle<-data.frame(cbind(c("enmalle.poiss0","enmalle.poiss1",
                                 "enmalle.negbin0","enmalle.negbin1","enmalle.negbinm1",
                                 "enmalle.zipoiss1","enmalle.zipoiss2","enmalle.zinb1","enmalle.zinb2"),
                               logLik=sapply(fm.enmalle,function (x) round(logLik(x),digits=0)),
                               Df=sapply(fm.enmalle,function (x) attr(logLik(x),"df"))))
colnames(logliks.enmalle)<-c("Model","logLik","DF")
logliks.enmalle$logLik<-round(as.numeric(levels(logliks.enmalle$logLik))[logliks.enmalle$logLik],digits = 0)
logliks.enmalle$DF<-round(as.numeric(levels(logliks.enmalle$DF))[logliks.enmalle$DF],digits = 2)
logliks.enmalle[order(-logliks.enmalle$logLik),]

#8.2.7.4# Vuong test to compare NB vs ZINB
# do not run with random models

vuong(enmalle.zinb1,enmalle.negbin1)
vuong(enmalle.zinb2,enmalle.negbin1)
vuong(enmalle.zinb1,enmalle.zinb2)

#8.2.7.5# Predicting probabilities

#no mixed models neither gams included in this comparison
phat.pois.enmalle<-predprob(enmalle.poiss1)
phat.pois.mn.enmalle<-apply(phat.pois.enmalle,2,mean)
phat.nb.enmalle<-predprob(enmalle.negbin1)
phat.nb.mn.enmalle<-apply(phat.nb.enmalle,2,mean)
phat.zipoiss1.enmalle<-predprob(enmalle.zipoiss1)
phat.zipoiss1.mn.enmalle<-apply(phat.zipoiss1.enmalle,2,mean)
phat.zipoiss2.enmalle<-predprob(enmalle.zipoiss2)
phat.zipoiss2.mn.enmalle<-apply(phat.zipoiss2.enmalle,2,mean)
phat.zinb1.enmalle<-predprob(enmalle.zinb1)
phat.zinb1.mn.enmalle<-apply(phat.zinb1.enmalle,2,mean)
phat.zinb2.enmalle<-predprob(enmalle.zinb2)
phat.zinb2.mn.enmalle<-apply(phat.zinb2.enmalle,2,mean)

with(enmalle.pouting,{
  hist(Ntot,prob=TRUE,col="grey",breaks=seq(-0.5,245,1),xlab="",main="")
  lines(x=seq(0,245,1),y=phat.pois.mn.enmalle,type="l",lwd=2,col="red")
  lines(x=seq(0,245,1),y=phat.nb.mn.enmalle,type="l",lwd=2,col="blue")
  lines(x=seq(0,245,1),y=phat.zipoiss1.mn.enmalle,type="l",lwd=2,col="green")
  lines(x=seq(0,245,1),y=phat.zinb1.mn.enmalle,type="l",lwd=2,col="yellow")
  lines(x=seq(0,245,1),y=phat.zipoiss2.mn.enmalle,type="l",lwd=2,col="grey")
  lines(x=seq(0,245,1),y=phat.zinb2.mn.enmalle,type="l",lwd=2,col="orange")
  
}) #???

#8.2.7.6# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

enmalle.aic
enmalle.bic
logliks.enmalle
logliks.enmalle2<-rbind(logLik=sapply(fm.enmalle,function (x) round(logLik(x),digits=0)),
                     Df=sapply(fm.enmalle,function (x) attr(logLik(x),"df")))

zeroes.enmalle<-c("Obs"=sum(enmalle.pouting$Ntot<1),
               "enmalle.poiss0"=sum(dpois(0,fitted(enmalle.poiss0))),
               "enmalle.poiss1"=sum(dpois(0,fitted(enmalle.poiss1))),
               "enmalle.negbin0"="NA",
               "enmalle.negbin1"=sum(dnbinom(0,mu=fitted(enmalle.negbin1),size=enmalle.negbin1$theta)),
               "enmalle.negbinm1"="NA",
               "enmalle.zipoiss1"=sum(predict(enmalle.zipoiss1,type="prob")[,1]),
               "enmalle.zipoiss2"=sum(predict(enmalle.zipoiss2,type="prob")[,1]),
               "enmalle.zinb1"=sum(predict(enmalle.zinb1,type="prob")[,1]),
               "enmalle.zinb2"=sum(predict(enmalle.zinb2,type="prob")[,1]))

modelsComp.enmalle<-data.frame(cbind(logliks.enmalle2[c(2,4,6,8,10,12,14,16,18)],enmalle.bic[[2]],enmalle.aic[[2]],logliks.enmalle2[c(1,3,5,7,9,11,13,15,17)],zeroes.enmalle[2:10]))
colnames(modelsComp.enmalle)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes.enmalle[[1]],")"))
modelsComp.enmalle$Df<-round(as.numeric(levels(modelsComp.enmalle$Df))[modelsComp.enmalle$Df],digits = 2)
modelsComp.enmalle$BIC<-round(as.numeric(levels(modelsComp.enmalle$BIC))[modelsComp.enmalle$BIC],digits = 2)
modelsComp.enmalle$AIC<-round(as.numeric(levels(modelsComp.enmalle$AIC))[modelsComp.enmalle$AIC],digits = 2)
modelsComp.enmalle$logLik<-round(as.numeric(levels(modelsComp.enmalle$logLik))[modelsComp.enmalle$logLik],digits = 0)
modelsComp.enmalle$Zeroes<-round(as.numeric(levels(modelsComp.enmalle$Zeroes))[modelsComp.enmalle$Zeroes],digits = 4)
modelsComp.enmalle$deltaAIC<-round(modelsComp.enmalle$AIC-min(modelsComp.enmalle$AIC),2)
modelsComp.enmalle$deltaBIC<-round(modelsComp.enmalle$BIC-min(modelsComp.enmalle$BIC),2)
modelsComp.enmalle$ModLikelihood<-round(exp(-modelsComp.enmalle$deltaAIC/2),2)
modelsComp.enmalle$AICweight<-round(modelsComp.enmalle$ModLikelihood/sum(modelsComp.enmalle$ModLikelihood),2)
modelsComp.enmalle<-modelsComp.enmalle[order(modelsComp.enmalle$deltaAIC),]
modelsComp.enmalle

#final candidate models
modelsComp.nasa
modelsComp.enmalle

#NO random effects
summary(nasa.negbin1)
Anova(nasa.negbin1)
summary(enmalle.zinb2)
Anova(enmalle.zinb2) #it does not work

#8.3# Selection of predictors

#8.3.1# nasa model selection
nasa.NB1<-nasa.negbin1
summary(nasa.NB1)
Anova(nasa.NB1)

#sstM
nasa.NB2<-glm.nb(Ntot~offset(log(OffSet))+fGRT+fyear+poly(Julian,2)+Depth+caladoNight,
                    data=nasa.pouting)
summary(nasa.NB2)
Anova(nasa.NB2)
lrtest(nasa.NB1,nasa.NB2)
AICtab(nasa.NB1,nasa.NB2)

#Depth
nasa.NB3<-glm.nb(Ntot~offset(log(OffSet))+fGRT+fyear+poly(Julian,2)+caladoNight,
                    data=nasa.pouting)
summary(nasa.NB3)
Anova(nasa.NB3)
lrtest(nasa.NB2,nasa.NB3)
AICtab(nasa.NB2,nasa.NB3)

AICtab(nasa.NB1,nasa.NB2,nasa.NB3)
BICtab(nasa.NB1,nasa.NB2,nasa.NB3)

coefplot2(nasa.NB2)
sum(residuals(nasa.NB2,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.NB2))+1))

#comparing models
AICnasa<-AIC(nasa.NB1,nasa.NB2,nasa.NB3,k=2) # AIC
BICnasa<-AIC(nasa.NB1,nasa.NB2,nasa.NB3,k=log(dim(nasa.pouting)[1])) # BIC

FMnasa<-list("M1"=nasa.NB1,"M2"=nasa.NB2,"M3"=nasa.NB3)
Logliksnasa<-rbind(logLik=sapply(FMnasa,function (x) round(logLik(x),digits=0)),
                     Df=sapply(FMnasa,function (x) attr(logLik(x),"df")))

ZEROnasa<-round(c("Obs"=sum(nasa.pouting$Ntot<1),
                  "M1"=sum(dnbinom(0,mu=fitted(nasa.NB1),size=nasa.NB1$theta)),
                  "M2"=sum(dnbinom(0,mu=fitted(nasa.NB2),size=nasa.NB2$theta)),
                  "M3"=sum(dnbinom(0,mu=fitted(nasa.NB3),size=nasa.NB3$theta))))

MCompnasa<-data.frame(cbind(Logliksnasa[c(2,4,6)],BICnasa[[2]],AICnasa[[2]],Logliksnasa[c(1,3,5)],ZEROnasa[2:4]))
colnames(MCompnasa)<-c("Df","BIC","AIC","logLik",paste("Zeroes (Obs=",zeroes.nasa[[1]],")"))
MCompnasa$deltaAIC<-round(MCompnasa$AIC-min(MCompnasa$AIC),2)
MCompnasa$deltaBIC<-round(MCompnasa$BIC-min(MCompnasa$BIC),2)
MCompnasa$ModLikelihood<-round(exp(-MCompnasa$deltaAIC/2),2)
MCompnasa$AICweight<-round(MCompnasa$ModLikelihood/sum(MCompnasa$ModLikelihood),2)
MCompnasa

summary(nasa.NB2)
Anova(nasa.NB2)
plot(allEffects(nasa.NB2))
sum(residuals(nasa.NB2,type="pearson")^2)/(nrow(nasa.pouting)-(length(coef(nasa.NB2))+1))
#slightly overdispersed

#8.3.2# enmalle model selection
enmalle.ZINB1 <- enmalle.zinb2
summary(enmalle.ZINB1)

#year:zone
enmalle.ZINB2 <- zeroinfl(Ntot~offset(log(OffSet))+
                            Gear*fGRT+
                            fyear+fZoneO+
                            poly(Julian,2)+
                            lDepth*Gear+
                            sstM+
                            caladoNight*Gear+
                            Seafloor|lDepth,
                          dist="negbin",link="logit",data=enmalle.pouting)
summary(enmalle.ZINB2)
lrtest(enmalle.ZINB1,enmalle.ZINB2)
AICtab(enmalle.ZINB1,enmalle.ZINB2)

#sstM, but reintroducing interaction
enmalle.ZINB3 <- zeroinfl(Ntot~offset(log(OffSet))+
                            Gear*fGRT+
                            fyear*fZoneO+
                            poly(Julian,2)+
                            lDepth*Gear+
                            caladoNight*Gear+
                            Seafloor|lDepth,
                          dist="negbin",link="logit",data=enmalle.pouting)
summary(enmalle.ZINB3)
lrtest(enmalle.ZINB1,enmalle.ZINB3)
AICtab(enmalle.ZINB1,enmalle.ZINB3)

#poly(Julian,2)
enmalle.ZINB4 <- zeroinfl(Ntot~offset(log(OffSet))+
                            Gear*fGRT+
                            fyear*fZoneO+
                            Julian+
                            lDepth*Gear+
                            caladoNight*Gear+
                            Seafloor|lDepth,
                          dist="negbin",link="logit",data=enmalle.pouting)
summary(enmalle.ZINB4)
lrtest(enmalle.ZINB3,enmalle.ZINB4)
AICtab(enmalle.ZINB3,enmalle.ZINB4)

#Seafloor
enmalle.ZINB5 <- zeroinfl(Ntot~offset(log(OffSet))+
                            Gear*fGRT+
                            fyear*fZoneO+
                            Julian+
                            lDepth*Gear+
                            caladoNight*Gear|lDepth,
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

summary(enmalle.ZINB4)
sum(residuals(enmalle.ZINB4,type="pearson")^2)/(nrow(enmalle.pouting)-(length(coef(enmalle.ZINB4))+1))
#slightly overdispersed

#final model selection
MCompnasa
MCompenmalle

summary(nasa.NB2)
summary(enmalle.ZINB4)

#8.4# Basic model checking

#8.4.1# nasa model checking
nasa.pouting$residNB<-residuals(nasa.NB2,type="pearson")
nasa.pouting$fitNB<-fitted(nasa.NB2,type="response")
nasa.pouting$LfitNB<-log(nasa.pouting$fitNB)

#8.4.1.1# Normality (not really relevant) and heterogeneity

hn<-ggplot(nasa.pouting, aes(x=residNB))+
  geom_histogram(binwidth = 0.2,aes(fill = ..count..))+
  scale_fill_gradient("Count", low = "green", high = "red") +
  geom_vline(aes(xintercept=0),col="black")+
  scale_x_continuous("Pearson residuals")
frn<-ggplot(data=nasa.pouting,aes(x=fitNB,y=residNB,colour=residNB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red")+
  scale_y_continuous("Pearson residuals")+
  scale_x_continuous("Fitted values")+
  geom_hline(aes(yintercept=0),col="black")+
  ggtitle("Residual vs. Fitted") # some lack of fit at low and high values
#we are no able to perfectly catch extreme values variation, too much overdispersion

pdf(file="NBassumption_nasapouting.pdf",width=10,height=10)

multiplot(hn,frn)

dev.off()

#8.4.1.2# Residuals vs predictors

rpn1<-ggplot(nasa.pouting, aes(x=factor(Idflota), y=residNB))+
  geom_boxplot(col="gray50",aes(fill = GRT))+
  scale_fill_gradient(low = "green", high = "red") +
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Boat")
rpn2<-ggplot(nasa.pouting, aes(x=fyear, y=residNB))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Year")+
  facet_wrap(~fZoneO,scales="free")
rpn3<-ggplot(nasa.pouting, aes(x=fyear, y=residNB))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Year")+
  facet_wrap(~Seafloor,scales="free")
rpn4<-ggplot(nasa.pouting, aes(x=Julian, y=residNB,colour=fitNB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Julian")
rpn5<-ggplot(nasa.pouting, aes(x=Depth, y=residNB,colour=fitNB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Depth(m)")
rpn6<-ggplot(nasa.pouting, aes(x=caladoNight, y=residNB,colour=fitNB))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("% Night time")

pdf(file="NBresid-expl_nasapouting.pdf",width=10,height=16)

multiplot(rpn2,rpn4,rpn1,rpn3,rpn5,rpn6,cols=2)

dev.off()

#8.4.1.3# Spatial patterns

ggplot(data=nasa.pouting,aes(x=Lon,y=Lat))+
  geom_point(aes(colour=residNB))+
  scale_colour_gradient(low="green",high="red")

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

pdf(file="NBresid-nasapouting.pdf",width=10,height=10)

multiplot(hn,var1.nasa,frn,var2.nasa,cols=2)

dev.off()

#8.4.2# enmalle model checking
enmalle.pouting$residZI<-residuals(enmalle.ZINB4,type="pearson")
enmalle.pouting$fitZI<-fitted(enmalle.ZINB4,type="response")
enmalle.pouting$LfitZI<-log(enmalle.pouting$fitZI)

#8.4.2.1# Normality (not really relevant) and heterogeneity

he<-ggplot(enmalle.pouting, aes(x=residZI))+
  geom_histogram(binwidth = 0.2,aes(fill = ..count..))+
  scale_fill_gradient("Count", low = "green", high = "red") +
  geom_vline(aes(xintercept=0),col="black")+
  scale_x_continuous("Pearson residuals") #very bad normality stuff
fre<-ggplot(data=enmalle.pouting,aes(x=LfitZI,y=residZI,colour=residZI))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red")+
  scale_y_continuous("Pearson residuals")+
  scale_x_continuous("Fitted values")+
  geom_hline(aes(yintercept=0),col="black")+
  ggtitle("Residual vs. Fitted") # some lack of fit at low and high values
#we are no able to perfectly catch extreme values variation, too much overdispersion

pdf(file="ZIassumption_enmallepouting.pdf",width=10,height=10)

multiplot(he,fre)

dev.off()

#8.4.2.2# Residuals vs predictors

rpe1<-ggplot(enmalle.pouting, aes(x=factor(Idflota), y=residZI))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Boat")+
  facet_wrap(~Gear,scales="free")+
  theme(axis.text.x = element_blank(), axis.text.x = element_blank())
rpe2<-ggplot(enmalle.pouting, aes(x=fyear, y=residZI))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Year",breaks=c("2000","2005","2010"))+
  facet_grid(fZoneO~Gear,scales="free")
rpe3<-ggplot(enmalle.pouting, aes(x=fyear, y=residZI))+
  geom_boxplot(col="gray50",fill = "gray70")+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_discrete("Year")+
  facet_grid(Seafloor~Gear,scales="free")
rpe4<-ggplot(enmalle.pouting, aes(x=Julian, y=residZI,colour=fitZI))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Julian")+
  facet_wrap(~Gear,scales="free")
rpe5<-ggplot(enmalle.pouting, aes(x=Depth, y=residZI,colour=fitZI))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("Depth(m)")+
  facet_wrap(~Gear,scales="free")
rpe6<-ggplot(enmalle.pouting, aes(x=caladoNight, y=residZI,colour=fitZI))+
  geom_point(size=3)+
  scale_colour_gradient(low="green", high="red",expression("CPUE"))+
  geom_hline(aes(yintercept=0),col="black")+
  scale_y_continuous("Residuals")+
  scale_x_continuous("% Night time")+
  facet_wrap(~Gear,scales="free")

pdf(file="ZIresid-expl_enmallepouting.pdf",width=10,height=16)

multiplot(rpe2,rpe4,rpe6,rpe3,rpe5,rpe1,cols=2)

dev.off()

#8.4.2.3# Spatial patterns

ggplot(data=enmalle.pouting,aes(x=Lon,y=Lat))+
  geom_point(aes(colour=residZI))+
  scale_colour_gradient(low="green",high="red")+
  facet_wrap(~Gear,scales="free")

enmalle.pouting.Spat<-data.frame(enmalle.pouting$residZI,enmalle.pouting$Lon,enmalle.pouting$Lat)
coordinates(enmalle.pouting.Spat)<- ~enmalle.pouting.Lon+enmalle.pouting.Lat
vario1.enmalle<-variogram(enmalle.pouting.residZI~1,data=enmalle.pouting.Spat)
plot(vario1.enmalle,pch=16,col=1,cex=1.5) # No spatial autocorrelation using the default
# classical method of moments estimate of the variogram (this assumes normality)
vario2.enmalle<-variogram(enmalle.pouting.residZI~1,data=enmalle.pouting.Spat,cressie=TRUE)
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

pdf(file="ZIresid-enmallepouting.pdf",width=10,height=10)

multiplot(he,var1.enmalle,fre,var2.enmalle,cols=2)

dev.off()

#8.5# Model coefficients interpretation (care with ln-transformation of Depth)??

summary(nasa.NB2)
summary(enmalle.ZINB4)

exp(coef(nasa.NB2))

exp(coef(enmalle.ZINB4)[c(1:63)]) #count part
exp(coef(enmalle.ZINB4)[64]) #zero part

# ----------------------------- #
#9# Bootstrapping the optimal model coefficients (zero & count parts)

#9.1# Function
#(add starting values to the model if needed!!)
dput(round(coef(nasa.NB2),4))
length(coef(nasa.NB2))

dput(round(coef(enmalle.ZINB4,"count"),4))
dput(round(coef(enmalle.ZINB4,"zero"),4))
length(coef(enmalle.ZINB4))

#http://www.ats.ucla.edu/stat/r/dae/zinbreg.htm

#9.1.1#nasa
#nasa.NB2<-glm.nb(Ntot~offset(log(OffSet))+fGRT+fyear+poly(Julian,2)+poly(Depth,2)+caladoNight,
#data=nasa.pouting)

boot.nb.nasa<-function (data,indices) {
  
  data<-data[indices,]
  try(mod<-glm.nb(Ntot~offset(log(OffSet))+fGRT+fyear+poly(Julian,2)+Depth+caladoNight,
	                data=data))
						   			  
	if (exists("mod")) { 
						   			 
	  as.vector(t(coef(summary(mod))[,1:2]))
	
	} else {rep(NA,times=28)}
  # If the above model crashes this fills in the gaps
  # with NA and the algorithm continues
	
}

#9.1.2#enmalle
#enmalle.ZINB4 <- zeroinfl(Ntot~offset(log(OffSet))+
#Gear*fGRT+fyear*fZoneO+Julian+poly(lDepth,2)*Gear+caladoNight*Gear+Seafloor|1,
#dist="negbin",link="logit",data=enmalle.pouting)

boot.zinb.enmalle<-function (data,indices) {
  
  data<-data[indices,]
  try(mod<-zeroinfl(Ntot~offset(log(OffSet))+
                      Gear*fGRT+
                      fyear*fZoneO+
                      Julian+
                      lDepth*Gear+
                      caladoNight*Gear+
                      Seafloor|lDepth,
                    dist="negbin",link="logit",data=data))
  
  if (exists("mod")) { 
    
    as.vector(t(do.call(rbind,coef(summary(mod)))[,1:2]))
    
  } else {rep(NA,times=128)}
  
}

#9.2# Coefficients (obtain CI of estimates excluding SE and theta)

#9.2.1# nasa
RR<-1000 # Number of resamples (reduce this number for testing the code)
nasa.nb.boot.out<-boot(data=nasa.pouting,statistic=boot.nb.nasa,R=RR)
nasa.nb.boot.out  # Basic output
plot(nasa.nb.boot.out,index=1)
#Example of histogram and qqplot for a given component

nasa.nb.boot.out2<-as.data.frame(nasa.nb.boot.out$t[,seq(1,25,2)]) # Coefficients of interest from the boot object matrix
colnames(nasa.nb.boot.out2)<-names(coef(nasa.NB2))
head(nasa.nb.boot.out2)
write.table(x=nasa.nb.boot.out2,file="nasabootCoefs.txt",row.names=F)

nasa.parmsPERC<-t(sapply(seq(1,25,2),function (i) {
  out<-boot.ci(nasa.nb.boot.out,index=i,type=c("perc"))
  with(out,c(Estimate=t0,pLow=percent[4],pUpp=percent[5]))
})) # Obtain intervals calculated using the bootstrap percentile method
row.names(nasa.parmsPERC)<-names(coef(nasa.NB2))
head(nasa.parmsPERC)
write.table(x=as.data.frame(nasa.parmsPERC),file="nasabootPERC.txt",row.names=T)

#it takes too long
#nasa.parmsBCA<-t(sapply(seq(1,27,2),function (i) {
#  out<-boot.ci(nasa.nb.boot.out,index=i,type=c("bca"))
#  with(out,c(Estimate=t0,bcaLow=bca[4],bcaUpp=bca[5]))
#})) # Obtain intervals calculated using the bootstrap bias-corrected accelerated method
#row.names(nasa.parmsBCA)<-names(coef(nasa.NB2))
#head(nasa.parmsBCA)
#write.table(x=as.data.frame(nasa.parmsBCA),file="nasabootBCA.txt",row.names=T)
#ci.nasa<-data.frame(cbind(nasa.parmsPERC,nasa.parmsBCA,confint(nasa.NB2)))
#ci.nasa<-ci.nasa[,c(1,2,3,5,6,7,8)]
#colnames(ci.nasa)<-c("Estimate","pLow","pUpp","bcaLow","bcaUpp","2.5%","97.%")

#we only present perc CI (BCA just in case reviewers required them)
ci.nasa<-data.frame(cbind(nasa.parmsPERC,confint(nasa.NB2)))
colnames(ci.nasa)<-c("Estimate","pLow","pUpp","2.5%","97.%")
str(ci.nasa)
ci.nasa
write.table(ci.nasa,file="ci_nasa.txt",row.names=T)

#9.2.2# enmalle
RR<-1000 # Number of resamples (reduce this number for testing the code)
enmalle.zinb.boot.out<-boot(data=enmalle.pouting,statistic=boot.zinb.enmalle,R=RR)
enmalle.zinb.boot.out  # Basic output
plot(enmalle.zinb.boot.out,index=1)

enmalle.zinb.boot.out2<-as.data.frame(enmalle.zinb.boot.out$t[,c(seq(1,121,2),125,127)]) # Coefficients of interest from the boot object matrix
colnames(enmalle.zinb.boot.out2)<-names(coef(enmalle.ZINB4))
head(enmalle.zinb.boot.out2)
write.table(x=enmalle.zinb.boot.out2,file="enmallebootCoefs.txt",row.names=F)

enmalle.parmsPERC<-t(sapply(c(seq(1,121,2),125,127),function (i) {
  out<-boot.ci(enmalle.zinb.boot.out,index=i,type=c("perc"))
  with(out,c(Estimate=t0,pLow=percent[4],pUpp=percent[5]))
})) # Obtain intervals calculated using the bootstrap percentile method
row.names(enmalle.parmsPERC)<-names(coef(enmalle.ZINB4))
head(enmalle.parmsPERC)
write.table(x=as.data.frame(enmalle.parmsPERC),file="enmallebootPERC.txt",row.names=T)

#it takes too long
#enmalle.parmsBCA<-t(sapply(seq(1,127,2),function (i) {
#  out<-boot.ci(enmalle.zinb.boot.out,index=i,type=c("bca"))
#  with(out,c(Estimate=t0,bcaLow=bca[4],bcaUpp=bca[5]))
#})) # Obtain intervals calculated using the bootstrap bias-corrected accelerated method
#row.names(enmalle.parmsBCA)<-names(coef(enmalle.ZINB4))
#head(enmalle.parmsBCA)
#write.table(x=as.data.frame(enmalle.parmsBCA),file="enmallebootBCA.txt",row.names=T)
#ci.enmalle<-data.frame(cbind(enmalle.parmsPERC,enmalle.parmsBCA,confint(enmalle.ZINB4)))
#ci.enmalle<-ci.enmalle[,c(1,2,3,5,6,7,8)]
#colnames(ci.enmalle)<-c("Estimate","pLow","pUpp","bcaLow","bcaUpp","2.5%","97.%")

#we only present perc CI (BCA just in case reviewers required them)
ci.enmalle<-data.frame(cbind(enmalle.parmsPERC,confint(enmalle.ZINB4)))
colnames(ci.enmalle)<-c("Estimate","pLow","pUpp","2.5%","97.%")
str(ci.enmalle)
ci.enmalle
write.table(ci.enmalle,file="ci_enmalle.txt",row.names=T)

#9.4# Histograms of all components

pdf(file="hist_bootnasa.pdf",width=10,height=12)

par(mfrow=c(5,3))
for (i in 1:13) {
  hist(nasa.nb.boot.out2[,i],breaks=50,col="light blue",main="",
      xlab=names(coef(nasa.NB2))[i])
  abline(v=coef(nasa.NB2)[i],col="red",lwd=2)
  }

dev.off()

pdf(file="hist_bootenmalle.pdf",width=20,height=20)

par(mfrow=c(8,8))
for (i in 1:63) {
  hist(enmalle.zinb.boot.out2[,i],breaks=50,col="light blue",main="",
       xlab=names(coef(enmalle.ZINB4))[i])
  abline(v=coef(enmalle.ZINB3)[i],col="red",lwd=2)
}

dev.off()
#check again this bootstrap

# ----------------------------- #
#10# Sketching results for the optimal model

#10.1# plotting nasa model
summary(coef(nasa.NB2))

#10.1.1# Continuous variables

rugscalN<-data.frame(rugscalN=unique(nasa.pouting$caladoNight))
rugsDoY<-data.frame(rugsDoY=unique(nasa.pouting$Julian))
rugsDepth<-data.frame(rugsDepth=unique(nasa.pouting$Depth))

nasa.precaladoNight<-predict(nasa.NB2,type="response",
                 newdata=data.frame(fGRT="C8-12",
                                    fyear="2006",
                                    Julian=183,
                                    Depth=mean(nasa.pouting$Depth,na.rm=T),
                                    caladoNight=seq(min(nasa.pouting$caladoNight),max(nasa.pouting$caladoNight),length=100), #caladonight
                                    OffSet=12))

nasa.precaladoNight<-data.frame(nasa.precaladoNight,seq(min(nasa.pouting$caladoNight),max(nasa.pouting$caladoNight),length=100))
colnames(nasa.precaladoNight)<-c("CalNPre","CalNSeq")
head(nasa.precaladoNight)
summary(nasa.precaladoNight)

nasa.calN<-ggplot(data=nasa.precaladoNight,aes(x=CalNSeq,y=CalNPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(2.5,55))+
  scale_x_continuous("% Night operation",limits=c(0,1))+
  #((55-2.5)*0.01)+2.5
  geom_segment(data=rugscalN,aes(x=rugscalN,xend=rugscalN,y=2.5,yend=3.025),stat="identity",lwd=0.5,col="gray50")
nasa.calN

nasa.preDepth<-predict(nasa.NB2,type="response",
                             newdata=data.frame(fGRT="C8-12",
                                                fyear="2006",
                                                Julian=183,
                                                Depth=seq(min(nasa.pouting$Depth),max(nasa.pouting$Depth),length=100), 
                                                caladoNight=mean(nasa.pouting$caladoNight,na.rm=T),
                                                OffSet=12))

nasa.preDepth<-data.frame(nasa.preDepth,seq(min(nasa.pouting$Depth),max(nasa.pouting$Depth),length=100))
colnames(nasa.preDepth)<-c("DepthPre","DepthSeq")
head(nasa.preDepth)
summary(nasa.preDepth)

nasa.Depth<-ggplot(data=nasa.preDepth,aes(x=DepthSeq,y=DepthPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(25,50))+
  scale_x_continuous("Depth (m)",limits=c(min(nasa.pouting$Depth),max(nasa.pouting$Depth)))+
  #((50-25)*0.01)+25
  geom_segment(data=rugsDepth,aes(x=rugsDepth,xend=rugsDepth,y=25,yend=25.25),stat="identity",lwd=0.5,col="gray50")
nasa.Depth

nasa.preJulian<-predict(nasa.NB2,type="response",
                       newdata=data.frame(fGRT="C8-12",
                                          fyear="2006",
                                          Julian=seq(1,365,1),
                                          Depth=mean(nasa.pouting$Depth,na.rm=T), 
                                          caladoNight=mean(nasa.pouting$caladoNight,na.rm=T),
                                          OffSet=12))

nasa.preJulian<-data.frame(nasa.preJulian,seq(1,365,1))
colnames(nasa.preJulian)<-c("JulianPre","JulianSeq")
head(nasa.preJulian)
summary(nasa.preJulian)

nasa.Julian<-ggplot(data=nasa.preJulian,aes(x=JulianSeq,y=JulianPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(35,145))+
  scale_x_continuous("Day of the year",limits=c(1,365))+
  #((145-35)*0.01)+35
  geom_segment(data=rugsDoY,aes(x=rugsDoY,xend=rugsDoY,y=35,yend=36.1),stat="identity",lwd=0.5,col="gray50")
nasa.Julian

#10.1.2# Categorical variables

nasa.preGRT<-data.frame(fGRT=c("C8","C8-12","C12"),
                     fyear="2006",
                     Julian=rep(183,times=3),
                     Depth=rep(mean(nasa.pouting$Depth,na.rm=T),3),
                     caladoNight=rep(mean(nasa.pouting$caladoNight,na.rm=T),times=3),
                     OffSet=rep(12,times=3))
nasa.preGRT<-predict(nasa.NB2,newdata=nasa.preGRT,type="response")

nasa.preGRT<-data.frame(cbind(nasa.preGRT,c("C8","C8-12","C12")))
colnames(nasa.preGRT)<-c("GRTPre","GRTSeq")
nasa.preGRT$GRTPre<-round(as.numeric(levels(nasa.preGRT$GRTPre))[nasa.preGRT$GRTPre],digits = 2)
nasa.preGRT$GRTSeq<-factor(nasa.preGRT$GRTSeq,levels = c("C8", "C8-12", "C12"))
str(nasa.preGRT)
summary(nasa.preGRT)

nasa.GRT<-ggplot(nasa.preGRT, aes(x = factor(GRTSeq), y = GRTPre))+
  geom_bar(stat = "identity",fill="blue",col="grey")+
  scale_y_continuous("Standardized Index",limits=c(0,60))+
  scale_x_discrete("GRT (tons)",labels=c("<8","8-12 tons",">12"))
nasa.GRT

pdf(file="NBfitvar-nasapouting.pdf",width=10,height=10)
multiplot(nasa.GRT,nasa.Depth,nasa.Julian,nasa.calN,cols=2)
dev.off()

nasa.preYear<-data.frame(fGRT=rep("C8-12",times=7),
                     fyear=c("2002","2005","2006","2007","2008","2009","2012"),
                     Julian=rep(183,times=7),
                     Depth=rep(mean(nasa.pouting$Depth,na.rm=T),7),
                     caladoNight=rep(mean(nasa.pouting$caladoNight,na.rm=T),times=7),
                     OffSet=rep(12,times=7))

nasa.preYear<-predict(nasa.NB2,newdata=nasa.preYear,type="response")

nasa.preYear<-data.frame(cbind(nasa.preYear,c(2002,2005,2006,2007,2008,2009,2012)))
colnames(nasa.preYear)<-c("YearPre","YearSeq")
str(nasa.preYear)
summary(nasa.preYear)

pdf(file="NBTrend-nasapouting.pdf",width=8,height=8)
ggplot(data=nasa.preYear,aes(x=YearSeq,y=YearPre))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index",limits=c(10,50))+
  scale_x_continuous("Year",limits=c(1999,2013))
dev.off()

nasa.Year<-ggplot(data=nasa.preYear,aes(x=YearSeq,y=YearPre))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index",limits=c(65,165))+
  scale_x_continuous("Year",limits=c(1999,2013))

nasapouting.abund<-lm(YearPre~YearSeq,data=nasa.preYear)
summary(nasapouting.abund) #Significant negative trend!
par(mfrow=c(2,2))
plot(nasapouting.abund)

#10.2# plotting enmalle model
summary(enmalle.ZINB4)

rugscalN2<-data.frame(rugscalN=unique(enmalle.pouting[enmalle.pouting$Gear=="MINOS",]$caladoNight))
rugscalN3<-data.frame(rugscalN=unique(enmalle.pouting[enmalle.pouting$Gear=="VETAS",]$caladoNight))
rugsDoY2<-data.frame(rugsDoY=unique(enmalle.pouting[enmalle.pouting$Gear=="MINOS",]$Julian))
rugsDoY3<-data.frame(rugsDoY=unique(enmalle.pouting[enmalle.pouting$Gear=="VETAS",]$Julian))
rugsDepth2<-data.frame(rugsDepth=unique(enmalle.pouting[enmalle.pouting$Gear=="MINOS",]$lDepth))
  rugsDepth2<-data.frame(exp(rugsDepth2))
rugsDepth3<-data.frame(rugsDepth=unique(enmalle.pouting[enmalle.pouting$Gear=="VETAS",]$lDepth))
  rugsDepth3<-data.frame(exp(rugsDepth3))

#10.2.1# ZI part

rugsDepth4<-data.frame(rugsDepth=unique(enmalle.pouting$lDepth))
  rugsDepth4<-data.frame(exp(rugsDepth4))

enmalle.prezero<-predict(enmalle.ZINB4,type="zero",
            newdata=data.frame(Gear="VETAS",
                               fGRT="C5-10",
                               fyear="2006",
                               Julian=183, #Julian
                               lDepth=seq(min(enmalle.pouting$lDepth),max(enmalle.pouting$lDepth),length.out = 100),
                               caladoNight=mean(enmalle.pouting$caladoNight,na.rm=T),
                               fZoneO="1",
                               Seafloor="hard",
                               OffSet=200)) # Estandarizacion a num por 200 m2 por h (200 m2 es el tamaÃ±o legal aproximado de un paÃ±o: 50 m de largo por 4 de alto)

enmalle.prezero<-data.frame(enmalle.prezero,
                            seq(min(enmalle.pouting$lDepth),max(enmalle.pouting$lDepth),length.out = 100))
colnames(enmalle.prezero)<-c("ZeroPre","DepthSeq")
enmalle.prezero$Depth<-exp(enmalle.prezero$DepthSeq)
str(enmalle.prezero)
summary(enmalle.prezero)
enmalle.prezero

pdf(file="enmalle_preZIzero.pdf")

ggplot(data=enmalle.prezero,aes(x=Depth,y=ZeroPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Probability of false zeros",limits=c(-0.02,1))+
  scale_x_continuous("ln-Depth (m)",limits=c(0,605))+
  geom_segment(data=rugsDepth4,aes(x=rugsDepth,xend=rugsDepth,
                                  y=-0.02,yend=-0.01),stat="identity",lwd=0.1,col="gray50")

dev.off()

#10.2.2# Continuous variables

enmalle.preJulian<-predict(enmalle.ZINB4,type="count",
                      newdata=data.frame(Gear="VETAS",
                                         fGRT="C5-10",
                                         fyear="2006",
                                         Julian=seq(1,365,1), #Julian
                                         lDepth=mean(enmalle.pouting$lDepth,na.rm=T),
                                         caladoNight=mean(enmalle.pouting$caladoNight,na.rm=T),
                                         fZoneO="1",
                                         Seafloor="hard",
                                         OffSet=200))

enmalle.preJulian<-data.frame(enmalle.preJulian,seq(1,365,1))
colnames(enmalle.preJulian)<-c("JulianPre","JulianSeq")
str(enmalle.preJulian)
summary(enmalle.preJulian)

enmalle.Julian<-ggplot(data=enmalle.preJulian,aes(x=JulianSeq,y=JulianPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance",limits=c(1,1.75))+
  scale_x_continuous("Day of the year",limits=c(1,365))+
  #((1.75-1)*0.01)+1
  geom_segment(data=rugsDoY2,aes(x=rugsDoY,xend=rugsDoY,y=1,yend=1.0075),stat="identity",lwd=0.5,col="gray50")
enmalle.Julian

enmalle.preDepth<-predict(enmalle.ZINB4,type="count",
                           newdata=data.frame(Gear=c(rep("VETAS",100),rep("MINOS",100)),
                                              fGRT="C5-10",
                                              fyear="2006",
                                              Julian=183,
                                              lDepth=rep(seq(min(enmalle.pouting$lDepth),max(enmalle.pouting$lDepth),length=100),2),
                                              caladoNight=mean(enmalle.pouting$caladoNight,na.rm=T),
                                              fZoneO="1",
                                              Seafloor="hard",
                                              OffSet=200))

enmalle.preDepth<-data.frame(enmalle.preDepth,
                             rep(seq(min(enmalle.pouting$lDepth),max(enmalle.pouting$lDepth),length=100),2),
                             c(rep("VETAS",100),rep("MINOS",100)))
colnames(enmalle.preDepth)<-c("DepthPre","DepthSeq","Gear")
enmalle.preDepth$Depth<-exp(enmalle.preDepth$DepthSeq)
str(enmalle.preDepth)
summary(enmalle.preDepth)

enmalle.Depth1<-ggplot(data=enmalle.preDepth[enmalle.preDepth$Gear=="MINOS",],aes(x=Depth,y=DepthPre))+
  geom_line(colour="blue",lwd=1,linetype=3)+
  scale_y_continuous("Standardized Abundance")+
  scale_x_continuous("Depth (m)",limits=c(1,605))+
  geom_segment(data=rugsDepth2,aes(x=rugsDepth,xend=rugsDepth,y=0,yend=0.001),stat="identity",lwd=0.5,col="gray50")
enmalle.Depth1

enmalle.Depth2<-ggplot(data=enmalle.preDepth[enmalle.preDepth$Gear=="VETAS",],aes(x=Depth,y=DepthPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance")+
  scale_x_continuous("Depth (m)",limits=c(1,605))+
  geom_segment(data=rugsDepth3,aes(x=rugsDepth,xend=rugsDepth,y=0,yend=0.03),stat="identity",lwd=0.5,col="gray50")
enmalle.Depth2

enmalle.precaladoNight<-predict(enmalle.ZINB4,type="count",
                          newdata=data.frame(Gear=c(rep("VETAS",100),rep("MINOS",100)),
                                             fGRT="C5-10",
                                             fyear="2006",
                                             Julian=183,
                                             lDepth=mean(enmalle.pouting$lDepth,na.rm=T),
                                             caladoNight=rep(seq(min(enmalle.pouting$caladoNight),max(enmalle.pouting$caladoNight),length=100),2),
                                             fZoneO="1",
                                             Seafloor="hard",
                                             OffSet=200))

enmalle.precaladoNight<-data.frame(enmalle.precaladoNight,
                             rep(seq(min(enmalle.pouting$caladoNight),max(enmalle.pouting$caladoNight),length=100),2),
                             c(rep("VETAS",100),rep("MINOS",100)))
colnames(enmalle.precaladoNight)<-c("CalNPre","CalNSeq","Gear")
str(enmalle.precaladoNight)
summary(enmalle.precaladoNight)

enmalle.calN1<-ggplot(data=enmalle.precaladoNight[enmalle.precaladoNight$Gear=="MINOS",],aes(x=CalNSeq,y=CalNPre))+
  geom_line(colour="blue",lwd=1,linetype=3)+
  scale_y_continuous("Standardized Abundance")+
  scale_x_continuous("Depth (m)",limits=c(0,1))+
  geom_segment(data=rugscalN2,aes(x=rugscalN,xend=rugscalN,y=0,yend=0.0005),stat="identity",lwd=0.5,col="gray50")
enmalle.calN1

enmalle.calN2<-ggplot(data=enmalle.precaladoNight[enmalle.precaladoNight$Gear=="VETAS",],aes(x=CalNSeq,y=CalNPre))+
  geom_line(colour="blue",lwd=1)+
  scale_y_continuous("Standardized Abundance")+
  scale_x_continuous("% Night",limits=c(0,1))+
  geom_segment(data=rugscalN3,aes(x=rugscalN,xend=rugscalN,y=0,yend=0.035),stat="identity",lwd=0.5,col="gray50")
enmalle.calN2

#10.2.3# Categorical variables

enmalle.preGRT<-data.frame(Gear=c(rep("MINOS",5),rep("VETAS",5)),
                             fGRT=rep(c("C5","C5-10","C10-15","C15-20","C>20"),2),
                             fyear=rep("2006",10),
                             Julian=rep(183,times=10),
                             lDepth=rep(mean(enmalle.pouting$lDepth,na.rm=T),times=10),
                             caladoNight=rep(mean(enmalle.pouting$caladoNight,na.rm=T),times=10),
                             fZoneO=rep("1",times=10),
                             Seafloor=rep("hard",times=10),
                             OffSet=rep(200,times=10))

enmalle.preGRT<-predict(enmalle.ZINB4,newdata=enmalle.preGRT,type="count")

enmalle.preGRT<-data.frame(cbind(enmalle.preGRT),
                           rep(c("C5","C5-10","C10-15","C15-20","C>20"),2),
                           c(rep("MINOS",5),rep("VETAS",5)))
colnames(enmalle.preGRT)<-c("GRTPre","GRTSeq","Gear")
str(enmalle.preGRT)
enmalle.preGRT

enmalle.GRT1<-ggplot(enmalle.preGRT[enmalle.preGRT$Gear=="MINOS",], aes(x = factor(GRTSeq), y = GRTPre))+
  geom_bar(stat = "identity",fill="grey",col="blue",linetype=3,lwd=1)+
  scale_y_continuous("Standardized Index")+
  scale_x_discrete("GRT (tons)",labels=c("<5","5-10","10-15","15-20",">20"))
enmalle.GRT1

enmalle.GRT2<-ggplot(enmalle.preGRT[enmalle.preGRT$Gear=="VETAS",], aes(x = factor(GRTSeq), y = GRTPre))+
  geom_bar(stat = "identity",fill="grey",col="blue",lwd=1)+
  scale_y_continuous("Standardized Index")+
  scale_x_discrete("GRT (tons)",labels=c("<5","5-10","10-15","15-20",">20"))
enmalle.GRT2

enmalle.GRT<-ggplot(enmalle.preGRT, aes(x = factor(GRTSeq), y = GRTPre))+
  geom_bar(stat = "identity",fill="blue",col="black")+
  scale_y_continuous("Standardized Index")+
  scale_x_discrete("GRT (tons)",labels=c("<5","5-10","10-15","15-20",">20"))+
  facet_wrap(~Gear,scales="free")
enmalle.GRT

enmalle.preSeafloor<-data.frame(Gear=rep("VETAS",3),
                           fGRT=rep("C5-10",3),
                           fyear=rep("2006",3),
                           Julian=rep(183,times=3),
                           lDepth=rep(mean(enmalle.pouting$lDepth,na.rm=T),times=3),
                           caladoNight=rep(mean(enmalle.pouting$caladoNight,na.rm=T),times=3),
                           fZoneO=rep("1",times=3),
                           Seafloor=c("hard","mixed","soft"),
                           OffSet=rep(200,times=3))

enmalle.preSeafloor<-predict(enmalle.ZINB4,newdata=enmalle.preSeafloor,type="count")

enmalle.preSeafloor<-data.frame(cbind(enmalle.preSeafloor),c("hard","mixed","soft"))
colnames(enmalle.preSeafloor)<-c("SeafloorPre","SeafloorSeq")
str(enmalle.preSeafloor)
enmalle.preSeafloor

enmalle.Seafloor<-ggplot(enmalle.preSeafloor, aes(x = factor(SeafloorSeq), y = SeafloorPre))+
  geom_bar(stat = "identity",fill="blue",col="grey")+
  scale_y_continuous("Standardized Index",limits=c(0,2.5))+
  scale_x_discrete("Sea-floor type",labels=c("Hard","Mixed","Soft"))
enmalle.Seafloor

pdf(file="ZINBfitvar-enmallepouting.pdf",width=8,height=12)
multiplot(enmalle.Depth1,enmalle.calN1,enmalle.GRT1,enmalle.Julian,enmalle.Depth2,enmalle.calN2,enmalle.GRT2,enmalle.Seafloor,cols=2)
dev.off()

enmalle.preYear<-data.frame(Gear=rep("VETAS",45),
                                fGRT=rep("C5-10",45),
                                fyear=rep(c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013"),3),
                                Julian=rep(183,times=45),
                                lDepth=rep(mean(enmalle.pouting$lDepth,na.rm=T),times=45),
                                caladoNight=rep(mean(enmalle.pouting$caladoNight,na.rm=T),times=45),
                                fZoneO=c(rep("1",times=15),rep("2",times=15),rep("3",times=15)),
                                Seafloor=rep("hard",45),
                                OffSet=rep(200,times=3))

enmalle.preYear<-predict(enmalle.ZINB4,newdata=enmalle.preYear,type="count")

enmalle.preYear<-data.frame(cbind(enmalle.preYear,
                                  rep(c(1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013),3),
                                  c(rep("Rias Baixas",times=15),rep("Golfo Artabro",times=15),rep("Cantabrico",times=15))))
colnames(enmalle.preYear)<-c("YearPre","YearSeq","Zones")
enmalle.preYear$YearPre<-round(as.numeric(levels(enmalle.preYear$YearPre))[enmalle.preYear$YearPre],digits = 2)
enmalle.preYear$YearSeq<-round(as.numeric(levels(enmalle.preYear$YearSeq))[enmalle.preYear$YearSeq],digits = 0)
enmalle.preYear$Zones<-factor(enmalle.preYear$Zones,levels=c("Rias Baixas","Golfo Artabro","Cantabrico"))
str(enmalle.preYear)
summary(enmalle.preYear)
enmalle.preYear

pdf(file="ZINBTrend-enmallepouting.pdf",width=12,height=8)
ggplot(data=enmalle.preYear,aes(x=YearSeq,y=YearPre))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index",limits=c(0,3.5))+
  scale_x_continuous("Year",breaks=c(2000,2005,2010),labels=c("2000","2005","2010"))+
  facet_wrap(~Zones)
dev.off()

enmallepouting.abund<-lm(YearPre~YearSeq*Zones,data=enmalle.preYear)
summary(enmallepouting.abund) 
anova(enmallepouting.abund) #No Significant trend!
par(mfrow=c(2,2))
plot(enmallepouting.abund)

#11# Plot of abundance, nominal cpue, and landings

#11.1# Calculate nominal cpue and average for the Rias Baixas

#nasa
head(nasa.pouting)

nasa.pouting$cpue<-(nasa.pouting$Ntot*12)/nasa.pouting$OffSet # Standardize at 50 pieces
cpues.nasa<-tapply(nasa.pouting$cpue,nasa.pouting$Year,mean,na.rm=T)
cpues.L.nasa<-tapply(nasa.pouting$cpue,nasa.pouting$Year,ci95Low)
cpues.U.nasa<-tapply(nasa.pouting$cpue,nasa.pouting$Year,ci95Upp)

cpues.nasa.pouting<-data.frame(cbind(cpues.nasa,cpues.L.nasa,cpues.U.nasa,c(2002,2005,2006,2007,2008,2009,2012)))
colnames(cpues.nasa.pouting)<-c("Index","ciLow","ciUpp","Year")
str(cpues.nasa.pouting)

limits <- aes(ymax = ciUpp, ymin= ciLow)
ggplot(data=cpues.nasa.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_errorbar(limits, colour="blue", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))
#do not match CPU with predicted values

#enmalle
head(enmalle.pouting)
enmalle.pouting

enmalle.pouting$cpue<-(enmalle.pouting$Ntot*200)/enmalle.pouting$OffSet # Standardize at 200 m2 paño
enmalle.pouting2<-enmalle.pouting[enmalle.pouting$Gear=="VETAS",]


enmalleRB<-enmalle.pouting2[enmalle.pouting2$ZoneO==1,]
enmalleAR<-enmalle.pouting2[enmalle.pouting2$ZoneO==2,]
enmalleCN<-enmalle.pouting2[enmalle.pouting2$ZoneO==3,]

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
cpues.enmalle.AR$Zones<-rep("Golfo Artabro",14)
cpues.enmalle.AR$Year<-as.numeric(seq(2000,2013,1))
colnames(cpues.enmalle.AR)<-c("Index","ciLow","ciUpp","Zones","Year")
str(cpues.enmalle.AR)

cpues.enmalle.CN<-tapply(enmalleCN$cpue,enmalleCN$Year,mean,na.rm=T)
cpues.L.enmalle.CN<-tapply(enmalleCN$cpue,enmalleCN$Year,ci95Low)
cpues.U.enmalle.CN<-tapply(enmalleCN$cpue,enmalleCN$Year,ci95Upp)
cpues.enmalle.CN<-data.frame(cbind(cpues.enmalle.CN,cpues.L.enmalle.CN,cpues.U.enmalle.CN))
cpues.enmalle.CN$Zones<-rep("Cantabrico",14)
cpues.enmalle.CN$Year<-as.numeric(seq(2000,2013,1))
colnames(cpues.enmalle.CN)<-c("Index","ciLow","ciUpp","Zones","Year")
str(cpues.enmalle.CN)

cpues.enmalle.pouting<-data.frame(rbind(cpues.enmalle.RB,cpues.enmalle.AR,cpues.enmalle.CN))
str(cpues.enmalle.pouting)
cpues.enmalle.pouting
cpues.enmalle.pouting$Zones <- ordered(cpues.enmalle.pouting$Zones,levels = c("Rias Baixas", "Golfo Artabro", "Cantabrico"))

limits <- aes(ymax = ciUpp, ymin= ciLow)
ggplot(data=cpues.enmalle.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_errorbar(limits, colour="blue", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
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
landings.pouting<-data.frame(cbind(zones,c(rep("Rias Baixas",18),rep("Golfo Artabro",18),rep("Cantabrico",18))))
colnames(landings.pouting)<-c("Year","Index","Zones")
str(landings.pouting)
head(landings.pouting)
landings.pouting$Zones <- ordered(landings.pouting$Zones,levels = c("Rias Baixas", "Golfo Artabro", "Cantabrico"))

ggplot(data=landings.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)

#11.3# Comparing abundance index, nominal cpue, landings

head(nasa.preYear)
head(enmalle.preYear)

t1a<-ggplot(data=nasa.preYear,aes(x=YearSeq,y=YearPre))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  ggtitle("Standardize abundance index (trap)")+
  theme(plot.title = element_text(size=26, face="bold"))
t1a

t1b<-ggplot(data=enmalle.preYear,aes(x=YearSeq,y=YearPre))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  scale_y_continuous("Standardized Index")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones)+
  ggtitle("Standardize abundance index (gillnet)")+
  theme(plot.title = element_text(size=26, face="bold"))
t1b

t2a<-ggplot(data=cpues.nasa.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  #geom_errorbar(limits, colour="blue", width=0)+
  scale_y_continuous("Nominal CPUE (nº/50pieces*h)")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  ggtitle("Nominal CPUE (trap)")+
  theme(plot.title = element_text(size=26, face="bold"))
t2a

t2b<-ggplot(data=cpues.enmalle.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
  #geom_errorbar(limits, colour="blue", width=0)+
  scale_y_continuous("Nominal CPUE (nº/200m2*h)")+
  scale_x_continuous("Year")+
  facet_wrap(~Zones)+  ggtitle("Nominal CPUE (gillnet)")+
  theme(plot.title = element_text(size=26, face="bold"))
t2b

t3<-ggplot(data=landings.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="blue")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3.5,col="blue")+
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

pdf(file="pouting_trends2.pdf",width=12,height=10)

multiplot(t1b,t2b,t3,cols=1)

dev.off()

#correlations trends

names(enmalle.preYear)
dim(enmalle.preYear)
names(cpues.enmalle.pouting)
dim(cpues.enmalle.pouting)
names(landings.pouting)
dim(landings.pouting)

landings.pouting2<-landings.pouting[landings.pouting$Year>1998&landings.pouting$Year<2014,]
dim(landings.pouting2)

#Rias Baixas
cor(enmalle.preYear[enmalle.preYear$Zones=="Rias Baixas",]$YearPre,
    cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Rias Baixas",]$Index,method="spearman")
cor(enmalle.preYear[enmalle.preYear$Zones=="Rias Baixas",]$YearPre,
    landings.pouting2[landings.pouting2$Zones=="Rias Baixas",]$Index,method="spearman")
cor(cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Rias Baixas",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Rias Baixas",]$Index,method="spearman")

#Artabro
cor(enmalle.preYear[enmalle.preYear$Zones=="Golfo Artabro"&enmalle.preYear$YearSeq>1999,]$YearPre,
    cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Golfo Artabro",]$Index,method="spearman")
cor(enmalle.preYear[enmalle.preYear$Zones=="Golfo Artabro",]$YearPre,
    landings.pouting2[landings.pouting2$Zones=="Golfo Artabro",]$Index,method="spearman")
cor(cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Golfo Artabro",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Golfo Artabro"&landings.pouting2$Year>1999,]$Index,method="spearman")

#Cantabrico
cor(enmalle.preYear[enmalle.preYear$Zones=="Cantabrico"&enmalle.preYear$YearSeq>1999,]$YearPre,
    cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Cantabrico",]$Index,method="spearman")
cor(enmalle.preYear[enmalle.preYear$Zones=="Cantabrico",]$YearPre,
    landings.pouting2[landings.pouting2$Zones=="Cantabrico",]$Index,method="spearman")
cor(cpues.enmalle.pouting[cpues.enmalle.pouting$Zones=="Cantabrico",]$Index,
    landings.pouting2[landings.pouting2$Zones=="Cantabrico"&landings.pouting2$Year>1999,]$Index,method="spearman")

#ridiculous correlation with oficial landings :(

# ----------------------------- #
#12# Plotting ZINB model predicted means for YEAR with Bootstrapping
# of predictions (95% CI) by means of shuffling residuals (Thierry Onkelinx code)

#12.1# Bootstrap

#nasa

newTrend<-data.frame(fGRT=rep("C8-12",times=7),
                     fyear=c("2002","2005","2006","2007","2008","2009","2012"),
                     Julian=rep(183,times=7),
                     Depth=rep(mean(nasa.pouting$Depth,na.rm=T),7),
                     caladoNight=rep(mean(nasa.pouting$caladoNight,na.rm=T),times=7),
                     OffSet=rep(12,times=7))
head(newTrend)

Fit<-predict(nasa.NB2,type="response")

Pearson<-residuals(nasa.NB2,type="pearson") # Pearson residuals
VarComp<-residuals(nasa.NB2,type="response")/Pearson # Raw residuals/Pearson residuals

Gear<-nasa.pouting$Gear
fGRT<-nasa.pouting$fGRT
fyear<-nasa.pouting$fyear
Julian<-nasa.pouting$Julian
Depth<-nasa.pouting$Depth
caladoNight<-nasa.pouting$caladoNight
OffSet<-nasa.pouting$OffSet

#bootstraping residuals
RR<-1000 # Number of resamples (reduce this number for testing the code)

nasa.bootstrap<-replicate(n=RR,{ 
  
  yStar<-pmax(round(Fit+sample(Pearson)*VarComp,0),0)
  try(mod<-glm.nb(yStar~offset(log(OffSet))+fGRT+fyear+poly(Julian,2)+Depth+caladoNight))
  
  if (exists("mod")) {
    
    predict(mod,newdata=newTrend,type="response")
    
  } else {rep(NA,times=15)} # If the above model crashes this fills in the gaps
  # with NA and the algorithm continues
  
})

CIs.nasa<-t(apply(X=nasa.bootstrap,MARGIN=1,FUN=quantile,c(0.025,0.975),na.rm=T))
colnames(CIs.nasa)<-c("ciLow","ciUpp")
newTrend$fit<-predict(nasa.NB2,newdata=newTrend,type="response")
newTrend<-cbind(newTrend,CIs.nasa)
newTrend$Year<-c(2002,2005,2006,2007,2008,2009,2012)
names(newTrend)
newTrend

#enmalle

newTrend2<-data.frame(Gear=rep("VETAS",45),
                     fGRT=rep("C5-10",45),
                     fyear=rep(c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013"),3),
                     Julian=rep(183,times=45),
                     lDepth=rep(mean(enmalle.pouting$lDepth,na.rm=T),times=45),
                     caladoNight=rep(mean(enmalle.pouting$caladoNight,na.rm=T),times=45),
                     fZoneO=c(rep("1",times=15),rep("2",times=15),rep("3",times=15)),
                     Seafloor=rep("hard",45),
                     OffSet=rep(200,times=3))
head(newTrend2)

Fit<-predict(enmalle.ZINB4,type="response")

Pearson<-residuals(enmalle.ZINB4,type="pearson") # Pearson residuals
VarComp<-residuals(enmalle.ZINB4,type="response")/Pearson # Raw residuals/Pearson residuals

Gear<-enmalle.pouting$Gear
fGRT<-enmalle.pouting$fGRT
fyear<-enmalle.pouting$fyear
Julian<-enmalle.pouting$Julian
lDepth<-enmalle.pouting$lDepth
caladoNight<-enmalle.pouting$caladoNight
fZoneO<-enmalle.pouting$fZoneO
Seafloor<-enmalle.pouting$Seafloor
OffSet<-enmalle.pouting$OffSet

#bootstraping residuals
RR<-1000 # Number of resamples (reduce this number for testing the code)

enmalle.bootstrap<-replicate(n=RR,{ 
	
    yStar<-pmax(round(Fit+sample(Pearson)*VarComp,0),0)
    try(mod<-zeroinfl(yStar~offset(log(OffSet))+
                        Gear*fGRT+
                        fyear*fZoneO+
                        Julian+
                        lDepth*Gear+
                        caladoNight*Gear+
                        Seafloor|lDepth,
                      dist="negbin",link="logit"))
    
    if (exists("mod")) {
    	
    	predict(mod,newdata=newTrend2,type="response")
    	
    } else {rep(NA,times=15)} # If the above model crashes this fills in the gaps
    						  # with NA and the algorithm continues

})

CIs.enmalle<-t(apply(X=enmalle.bootstrap,MARGIN=1,FUN=quantile,c(0.025,0.975),na.rm=T))
colnames(CIs.enmalle)<-c("ciLow","ciUpp")
newTrend2$fit<-predict(enmalle.ZINB4,newdata=newTrend2,type="response")
newTrend2<-cbind(newTrend2,CIs.enmalle)
newTrend2$Year<-seq(1999,2013,1)
names(newTrend2)
newTrend2

#12.2# Plot of abundance

#nasa

pdf(file="pouting_nasaYearPredCI.pdf",width=10,height=8)

limits <- aes(ymax = ciUpp, ymin= ciLow)
ggplot(data=newTrend,aes(x=Year,y=fit))+
  #geom_segment(aes(x=Year,y=fit,xend=Year,yend=ciUpp),col="black",lwd=0.5)+
  #geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=fit),col="black",lwd=0.5)+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_errorbar(limits, colour="black", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Standardized Index of Abundance")+
  scale_x_continuous("Year",limits=c(1999,2013),breaks=c(2000,2005,2010),labels=c(2000,2005,2010))

dev.off()

nasa.boot1<-ggplot(data=newTrend,aes(x=Year,y=fit))+
  #geom_segment(aes(x=Year,y=fit,xend=Year,yend=ciUpp),col="black",lwd=0.5)+
  #geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=fit),col="black",lwd=0.5)+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_errorbar(limits, colour="black", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Standardized Index of Abundance")+
  scale_x_continuous("Year",limits=c(1999,2013),breaks=c(2000,2005,2010),labels=c(2000,2005,2010))+
  ggtitle("Standardized abundance index")+
  theme(plot.title = element_text(size=26, face="bold"))
nasa.boot1

cpues.nasa.pouting2<-cpues.nasa.pouting
min(cpues.nasa.pouting2$ciLow)
#cpues.enmalle.pouting2$ciLow<-ifelse(cpues.enmalle.pouting2$ciLow<0,0,cpues.enmalle.pouting2$ciLow)
limits <- aes(ymax = ciUpp, ymin= ciLow)
nasa.cpue1<-ggplot(data=cpues.nasa.pouting2,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_errorbar(limits, colour="black", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Nominal CPUE (nº/50traps*h)")+
  scale_x_continuous("Year",limits=c(1999,2013),breaks=c(2000,2005,2010),labels=c(2000,2005,2010))+
  ggtitle("Nominal CPUE")+
  theme(plot.title = element_text(size=26, face="bold"))
nasa.cpue1

landings.pouting2<-read.table(file="faneca_landings.txt",header=T,dec=".")
landings.pouting2<-landings.pouting2[,c(1,2)]

nasa.landing1<-ggplot(data=landings.pouting2,aes(x=year,y=totland))+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Total landings (kg)")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  ggtitle("Oficial landings")+
  theme(plot.title = element_text(size=26, face="bold"))
nasa.landing1

multiplot(nasa.boot1,nasa.cpue1,nasa.landing1,cols=1)

pdf(file="pouting-nasa_trends.pdf",width=15,height=20)

multiplot(nasa.boot1,nasa.cpue1,nasa.landing1,cols=1)

dev.off()

#enmalle

levels(newTrend2$fZoneO)<-c("1"="Rias Baixas", "2"="Golfo Artabro","3"="Cantabrico")
levels(newTrend2$fZoneO)

pdf(file="pouting_enmalleYearPredCI.pdf",width=10,height=8)

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
  facet_wrap(~fZoneO)

dev.off()

enmalle.boot1<-ggplot(data=newTrend2,aes(x=Year,y=fit))+
  #geom_segment(aes(x=Year,y=fit,xend=Year,yend=ciUpp),col="black",lwd=0.5)+
  #geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=fit),col="black",lwd=0.5)+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_errorbar(limits, colour="black", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Standardized Index of Abundance")+
  scale_x_continuous("Year")+
  facet_wrap(~fZoneO,scales="free_y")+
  ggtitle("Standardized abundance index")+
  theme(plot.title = element_text(size=26, face="bold"))
enmalle.boot1

cpues.enmalle.pouting2<-cpues.enmalle.pouting
min(cpues.enmalle.pouting2$ciLow)
cpues.enmalle.pouting2$ciLow<-ifelse(cpues.enmalle.pouting2$ciLow<0,0,cpues.enmalle.pouting2$ciLow)
limits <- aes(ymax = ciUpp, ymin= ciLow)
enmalle.cpue1<-ggplot(data=cpues.enmalle.pouting2,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_errorbar(limits, colour="black", width=0)+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Nominal CPUE (nº/200m2*h)")+
  scale_x_continuous("Year")+
  facet_wrap(~Zones,scales="free_y")+
  ggtitle("Nominal CPUE")+
  theme(plot.title = element_text(size=26, face="bold"))
enmalle.cpue1

enmalle.landing1<-ggplot(data=landings.pouting,aes(x=Year,y=Index))+
  geom_line(lwd=0.5,linetype="dotted",col="black")+
  geom_point(size=5,col="gray90")+
  geom_point(size=3,col="black")+
  scale_y_continuous("Total landings (kg)")+
  scale_x_continuous("Year",limits=c(1999,2013))+
  facet_wrap(~Zones,scales="free_y")+
  ggtitle("Oficial landings")+
  theme(plot.title = element_text(size=26, face="bold"))
enmalle.landing1

multiplot(enmalle.boot1,enmalle.cpue1,enmalle.landing1,cols=1)

pdf(file="pouting-enmalle_trends.pdf",width=15,height=20)

multiplot(enmalle.boot1,enmalle.cpue1,enmalle.landing1,cols=1)

dev.off()
