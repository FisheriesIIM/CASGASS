##########################################################
# Analysis of Pollachius pollachius catches in Galicia   #
# using data sampled onboard fishing vessels by the UTPB #
##########################################################

library(pscl)
library(mgcv)
library(MASS)
library(effects)
library(ggplot2)
library(boot)
library(lmtest)
library(sp)
library(gstat)
quartz.options(dpi=75)

source(file="/Users/jaimeoterovillar/Documents/Bibliografia/Estadistica/R/Funciones/HighStatLib.r")
source(file="/Users/jaimeoterovillar/Documents/Bibliografia/Estadistica/R/Funciones/multiple_ggplot.r")
source(file="/Users/jaimeoterovillar/Documents/Bibliografia/Estadistica/R/Funciones/CI_mean.r")


# ----------------------------- #
#1# Load pollack data

pollack<-read.table(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/Data/Catches/Working files/2_abadejo.txt",header=T,dec=".",sep=",")

head(pollack)
dim(pollack)
length(unique(pollack$Idflota)) # Number of vessels surveyed
sum(pollack$Ntot,na.rm=T) # Number of individuals
length(unique(pollack$Gear)) # Number of gears


# ----------------------------- #
#2# Evaluate fishing per gear, distribution of response and other elementary information

aba.tot<-sum(pollack$Ntot,na.rm=T) # Total number of caught pollack
abadejos<-as.numeric(tapply(pollack$Ntot,pollack$Gear,sum,na.rm=T)) # Total number per gear
gear.hauls<-tapply(pollack$Ntot,pollack$Gear,length) # Hauls per gear
fish<-ifelse(pollack$Ntot==0,0,1) # Variable to classify 0/1 hauls
zero.fish<-gear.hauls-tapply(fish,pollack$Gear,sum,na.rm=T) # Zeros per gear
gear.zeros<-round((zero.fish*100)/gear.hauls,2) # Percentage of zeros per gear

basic.info<-data.frame(cbind(gear.hauls,gear.zeros,(abadejos*100)/aba.tot))
colnames(basic.info)<-c("hauls","zeros","catch")
basic.info$gear<-c("Boliche","Bou de Man","Bou de Vara","Liña-Cordel","Miño",
	"Nasa Choco","Nasa Nécora-Camarón","Nasa Peixes","Palangrillo",
	"Racú","Raeira","Rasco","Rastro Camarón","Trasmallo","Trueiro","Veta")

pdf(file="/Users/jaimeoterovillar/Desktop/gearsPollack.pdf",width=15,height=8)

ggplot(data=basic.info,aes(x=gear,y=catch))+
	geom_bar(stat="identity")+
	scale_y_continuous(limits=c(0,100),"Frequency")+
	scale_x_discrete("Gear")+
	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=14),
		  axis.text.y=element_text(size=14))+
	geom_text(aes(label=paste("[",hauls,",",zeros,"]"),vjust=-0.5),size=3)

dev.off()

par(mar=c(5,5,3,3))
plot(table(pollack$Ntot),ylab="Frequency",xlab="Number of pollack caught per haul")


# ----------------------------- #
#3# Select data set that only include the 2 main gears

pollack1<-pollack[pollack$Gear=="MINOS" | pollack$Gear=="VETAS",]

head(pollack1)
dim(pollack1)


# ----------------------------- #
#4# Add variables and recode as factors when needed

#4.1# Variables already in the data (remove NAs)

summary(pollack1)

pollack1<-pollack1[-which(is.na(pollack1$Lat)),]
pollack1<-pollack1[-which(is.na(pollack1$Depth)),]
pollack1<-pollack1[-which(is.na(pollack1$Seafloor)),]

str(pollack1)

pollack1$fyear<-factor(pollack1$Year)
pollack1$fgear<-factor(pollack1$Gear)
pollack1$fcrew<-factor(ifelse(pollack1$Crew <= 3,1,2)) # Crew as factor
pollack1$fZoneA<-factor(pollack1$ZoneA)
pollack1$fZoneO<-factor(pollack1$ZoneO)

#4.2# Gear size

#4.2.1# Recorded gear size in each trip

gearSize<-read.csv2(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/Data/Catches/Original files/dimensiones_artes_pesca.csv",header=T,dec=",",sep=",")

head(gearSize)

gearSize1<-gearSize[gearSize$ARTE=="MINOS" | gearSize$ARTE=="VETAS",] # Only Miños y Vetas
gearSize1<-gearSize1[,c(1:3,5)]

head(gearSize1)
dim(gearSize1)

summary(gearSize1$longitud);hist(gearSize1$longitud) # Looks OK
summary(gearSize1$altura);hist(gearSize1$altura) # Strange values

gearSize1$altura<-ifelse(gearSize1$altura==2500,2.5,ifelse(gearSize1$altura < 5,gearSize1$altura,gearSize1$altura/100)) # Remove mistakes and scale in m

pollack1<-merge(pollack1,gearSize1,by.x=c("Idlance","Gear"),by.y=c("Idlance","ARTE")) # Merge datasets

head(pollack1)

ids<-unique(pollack1$Idflota) # Fill in NAs with vessel-average data
long<-length(ids)
lista<-as.list(1:long)
names(lista)<-ids

for (i in 1:long) {lista[[i]]<-pollack1[pollack1$Idflota==ids[i],]}
for (i in 1:long) {lista[[i]]$longitud2<-ifelse(sum(lista[[i]]$longitud,na.rm=T) > 0,mean(lista[[i]]$longitud,na.rm=T),NA)}
for (i in 1:long) {lista[[i]]$altura2<-ifelse(sum(lista[[i]]$altura,na.rm=T) > 0,mean(lista[[i]]$altura,na.rm=T),NA)}

pollack1<-do.call(rbind,lista)

head(pollack1)

pollack1$longitud3<-ifelse(is.na(pollack1$longitud),pollack1$longitud2,
	pollack1$longitud)
pollack1$altura3<-ifelse(is.na(pollack1$altura),pollack1$altura2,pollack1$altura)

pollack1$longitud4<-ifelse(is.na(pollack1$longitud3),
	mean(pollack1$longitud3,na.rm=T),pollack1$longitud3) # Fill in with fleet average data for those vessels without any information
pollack1$altura4<-ifelse(is.na(pollack1$altura3),
	mean(pollack1$altura3,na.rm=T),pollack1$altura3) 

pollack1$Area<-(pollack1$Pieces*pollack1$longitud4*pollack1$altura4) # MIÑOS & VETAS Area (m2) 

#4.2.2# Common value (according to law) for all MIÑOS and VETAS

# pollack1$Area<-ifelse(pollack1$Gear=="MINOS",pollack1$Pieces*49*4.46,
#					  pollack1$Pieces*52.3*3.67)/1000000 # MIÑOS & VETAS Area (km2) 

#4.2.3# Differences between Vetas and Miños

head(pollack1)

gear1<-lm(log(Area)~Gear,data=pollack1)
summary(gear1)
plot(allEffects(gear1))

gear2<-lm(log(Depth)~Gear,data=pollack1)
summary(gear2)
plot(allEffects(gear2))

gear3<-lm(log(Soak)~Gear,data=pollack1)
summary(gear3)
plot(allEffects(gear3))


# ----------------------------- #
#5# Add environmental data (At this time I run models without AR structure)

#5.1# Load upwelling and remove seasonal cycle (Esta variable no la usaremos
# por las dificultades a la hora de interpretarla...)

qx<-read.csv2(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Upwelling.csv",header=T,dec=",",sep=",")

head(qx)
qx$QX<-qx$QX/1000 # Change scale

summary(qx$QX)
IntQua<-3*(summary(qx$QX)[[5]]-summary(qx$QX)[[2]])

hlow<-summary(qx$QX)[[2]]-IntQua
hupp<-summary(qx$QX)[[5]]+IntQua

qx$QX<-ifelse(qx$QX < hlow,hlow,qx$QX)
qx$QX<-ifelse(qx$QX > hupp,hupp,qx$QX)

g1<-gam(QX~s(DoY,k=6,bs="cc")+s(Timer,k=4,bs="cr"),data=qx,na.action=na.exclude)
summary(g1)
#plot(g1,pages=1)

qx$qxAno<-residuals(g1)

#5.2# Load sst and remove seasonal cycle for each zone

sst<-read.table(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/Data/Oceanography/NOAA_OISSTv2_daily/sstGal.selection.txt",header=T,dec=".",sep=" ")

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

win.Up<-6 # 6 days for upwelling (the previous weeks)
win.Te<-6 # Days for sst: 3 (the previous 3 days), 6 (the previous week), 12 (the two previous weeks)

#5.4.2# Calculate environmental covariate (from the window day to the preceeding day of the catch)

pollack1$QxM<-as.vector(rep(NA,dim(pollack1)[1]))

for (i in 1:dim(pollack1)[1]) {pollack1$QxM[i]<-mean(oceano$qxAno[(which(oceano$Year==pollack1$Year[i] &
	oceano$DoY==pollack1$Julian[i])-win.Up):(which(oceano$Year==pollack1$Year[i] &
	oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)}

pollack1$sstM<-as.vector(rep(NA,dim(pollack1)[1]))

for (i in 1:dim(pollack1)[1]) {pollack1$sstM[i]<-
	
	if (pollack1$ZoneA[i]=="1") {mean(oceano$sstAnoZ1[(which(oceano$Year==pollack1$Year[i] &
		oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
		oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)} 
	
	else {if (pollack1$ZoneA[i]=="2") {mean(oceano$sstAnoZ2[(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)} 
		
	else {if (pollack1$ZoneA[i]=="3") {mean(oceano$sstAnoZ3[(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)} 
			
	else {if (pollack1$ZoneA[i]=="4") {mean(oceano$sstAnoZ4[(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)} 
				
	else {if (pollack1$ZoneA[i]=="5") {mean(oceano$sstAnoZ5[(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)} 
					
	else {if (pollack1$ZoneA[i]=="6") {mean(oceano$sstAnoZ6[(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)} 
						
	else {if (pollack1$ZoneA[i]=="7") {mean(oceano$sstAnoZ7[(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)} 
							
	else {if (pollack1$ZoneA[i]=="8") {mean(oceano$sstAnoZ8[(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
			oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)} 
								
	else {mean(oceano$sstAnoZ9[(which(oceano$Year==pollack1$Year[i] &
		oceano$DoY==pollack1$Julian[i])-win.Te):(which(oceano$Year==pollack1$Year[i] &
		oceano$DoY==pollack1$Julian[i])-1)],na.rm=T)}
								
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

head(pollack1)
summary(pollack1)

#pollack1<-pollack1[-which(is.na(pollack1$QxM)),] # Remove Qx NAs (only if we were to use this variable)

which(pollack1$Depth > 200) # Solo un lance por encima de 200 m (603 m de profundidad): la utpb dice que es real, pero que es de un fulano que va a rape y no representa el artesanal con que lo borro

pollack1<-pollack1[-5561,]

names(pollack1)

pollack1<-pollack1[,-c(30:37)]
dim(pollack1)

rep.lances<-tapply(pollack1$Ntot,pollack1$Idlance,length)
which(rep.lances > 1)
rep.lances[2530]
pollack1[pollack1$Idlance=="5069",] # Miño y veta unidos en el mismo lance
rep.lances[2552]
pollack1[pollack1$Idlance=="5143",] # Miño y veta unidos en el mismo lance

#6.1# Response variable

#6.1.1# Percentage of zeroes and other stuff

(length(which(pollack1$Ntot==0))*100)/nrow(pollack1) # Percentage of zeros
dim(pollack1)[1] # Number of hauls
length(unique(pollack1$Idlance)) # 2 Hauls repetidos aclarado anteriormente
length(unique(pollack1$Idflota)) # Number of vessels
length(unique(pollack1$Idjornada)) # Fishing trips

poll.dist<-as.data.frame(table(pollack1$Ntot))
colnames(poll.dist)<-c("Count","Freq")

pdf(file="/Users/jaimeoterovillar/Desktop/respPollack.pdf",width=15,height=6)

ggplot(data=poll.dist,aes(x=Count,y=Freq))+
	geom_bar(stat="identity")+
	scale_y_continuous(limits=c(0,4100),"Frequency")+
	scale_x_discrete("Number of pollack caught per haul")+
	theme(axis.text.x=element_text(size=7),
		  axis.text.y=element_text(size=8))

dev.off()

par(mar=c(5,5,3,3))
hist(pollack1$Ntot,prob=T,breaks=seq(-0.5,247.5,1),ylab="Probability",xlab="Nº of fished pollack",main="")

#6.1.2# Plots of spatial distribution of effort (m2 of net and fishing hour) 

ggplot(data=pollack1,aes(x=Lon,y=Lat))+geom_point(aes(color=Area*Soak))

#6.1.3# Distribution of cases across explanatory variables

table(pollack1$Year) # Few data in 1999
table(pollack1$Year,pollack1$Month)
table(pollack1$fcrew)
table(pollack1$fgear)
table(pollack1$ZoneA)
table(pollack1$ZoneO)
table(pollack1$Year,pollack1$ZoneO)

#6.1.4# Vessels per year to compare with Despachos (vetas+miños)

v99<-length(unique(pollack1$Idflota[pollack1$Year=="1999"]))
v00<-length(unique(pollack1$Idflota[pollack1$Year=="2000"]))
v01<-length(unique(pollack1$Idflota[pollack1$Year=="2001"]))
v02<-length(unique(pollack1$Idflota[pollack1$Year=="2002"]))
v03<-length(unique(pollack1$Idflota[pollack1$Year=="2003"]))
v04<-length(unique(pollack1$Idflota[pollack1$Year=="2004"]))
v05<-length(unique(pollack1$Idflota[pollack1$Year=="2005"]))
v06<-length(unique(pollack1$Idflota[pollack1$Year=="2006"]))
v07<-length(unique(pollack1$Idflota[pollack1$Year=="2007"]))
v08<-length(unique(pollack1$Idflota[pollack1$Year=="2008"]))
v09<-length(unique(pollack1$Idflota[pollack1$Year=="2009"]))
v10<-length(unique(pollack1$Idflota[pollack1$Year=="2010"]))
v11<-length(unique(pollack1$Idflota[pollack1$Year=="2011"]))
v12<-length(unique(pollack1$Idflota[pollack1$Year=="2012"]))
v13<-length(unique(pollack1$Idflota[pollack1$Year=="2013"]))

vessels<-c(v03,v04,v05,v06,v07,v08,v09,v10,v11,v12,v13)
despachos<-c(929,1029,1105,1062,1156,1170,1208,1054,979,963,934)

represent1<-(vessels*100)/despachos

#6.1.5# Hauls per year

h99<-length(unique(pollack1$Idlance[pollack1$Year=="1999"]))
h00<-length(unique(pollack1$Idlance[pollack1$Year=="2000"]))
h01<-length(unique(pollack1$Idlance[pollack1$Year=="2001"]))
h02<-length(unique(pollack1$Idlance[pollack1$Year=="2002"]))
h03<-length(unique(pollack1$Idlance[pollack1$Year=="2003"]))
h04<-length(unique(pollack1$Idlance[pollack1$Year=="2004"]))
h05<-length(unique(pollack1$Idlance[pollack1$Year=="2005"]))
h06<-length(unique(pollack1$Idlance[pollack1$Year=="2006"]))
h07<-length(unique(pollack1$Idlance[pollack1$Year=="2007"]))
h08<-length(unique(pollack1$Idlance[pollack1$Year=="2008"]))
h09<-length(unique(pollack1$Idlance[pollack1$Year=="2009"]))
h10<-length(unique(pollack1$Idlance[pollack1$Year=="2010"]))
h11<-length(unique(pollack1$Idlance[pollack1$Year=="2011"]))
h12<-length(unique(pollack1$Idlance[pollack1$Year=="2012"]))
h13<-length(unique(pollack1$Idlance[pollack1$Year=="2013"]))

hauls<-c(h99,h00,h01,h02,h03,h04,h05,h06,h07,h08,h09,h10,h11,h12,h13)

#6.1.6# Fishes per year

f99<-sum(pollack1$Ntot[pollack1$Year=="1999"],na.rm=T)
f00<-sum(pollack1$Ntot[pollack1$Year=="2000"],na.rm=T)
f01<-sum(pollack1$Ntot[pollack1$Year=="2001"],na.rm=T)
f02<-sum(pollack1$Ntot[pollack1$Year=="2002"],na.rm=T)
f03<-sum(pollack1$Ntot[pollack1$Year=="2003"],na.rm=T)
f04<-sum(pollack1$Ntot[pollack1$Year=="2004"],na.rm=T)
f05<-sum(pollack1$Ntot[pollack1$Year=="2005"],na.rm=T)
f06<-sum(pollack1$Ntot[pollack1$Year=="2006"],na.rm=T)
f07<-sum(pollack1$Ntot[pollack1$Year=="2007"],na.rm=T)
f08<-sum(pollack1$Ntot[pollack1$Year=="2008"],na.rm=T)
f09<-sum(pollack1$Ntot[pollack1$Year=="2009"],na.rm=T)
f10<-sum(pollack1$Ntot[pollack1$Year=="2010"],na.rm=T)
f11<-sum(pollack1$Ntot[pollack1$Year=="2011"],na.rm=T)
f12<-sum(pollack1$Ntot[pollack1$Year=="2012"],na.rm=T)
f13<-sum(pollack1$Ntot[pollack1$Year=="2013"],na.rm=T)

fishes<-c(f99,f00,f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13)

#6.1.7# Days per year to compare with Jornadas despachadas (vetas+miños)

j99<-length(unique(pollack1$Idjornada[pollack1$Year=="1999"]))
j00<-length(unique(pollack1$Idjornada[pollack1$Year=="2000"]))
j01<-length(unique(pollack1$Idjornada[pollack1$Year=="2001"]))
j02<-length(unique(pollack1$Idjornada[pollack1$Year=="2002"]))
j03<-length(unique(pollack1$Idjornada[pollack1$Year=="2003"]))
j04<-length(unique(pollack1$Idjornada[pollack1$Year=="2004"]))
j05<-length(unique(pollack1$Idjornada[pollack1$Year=="2005"]))
j06<-length(unique(pollack1$Idjornada[pollack1$Year=="2006"]))
j07<-length(unique(pollack1$Idjornada[pollack1$Year=="2007"]))
j08<-length(unique(pollack1$Idjornada[pollack1$Year=="2008"]))
j09<-length(unique(pollack1$Idjornada[pollack1$Year=="2009"]))
j10<-length(unique(pollack1$Idjornada[pollack1$Year=="2010"]))
j11<-length(unique(pollack1$Idjornada[pollack1$Year=="2011"]))
j12<-length(unique(pollack1$Idjornada[pollack1$Year=="2012"]))
j13<-length(unique(pollack1$Idjornada[pollack1$Year=="2013"]))

daysD<-c(j03,j04,j05,j06,j07,j08,j09,j10,j11,j12,j13)
jornadas<-c(91481,113782,120188,114164,117539,116286,117552,100601,91303,90902,89674)

represent2<-(daysD*100)/jornadas

#6.2# Set up the offset

#6.2.1# Exploration of effort variables

head(pollack1)

par(mfrow=c(2,3))
plot(GRT~Crew,data=pollack1)
plot(GRT~Area,data=pollack1)
plot(GRT~Soak,data=pollack1)
plot(Crew~Area,data=pollack1)
plot(Crew~Soak,data=pollack1)
plot(Area~Soak,data=pollack1)
quartz()
par(mfrow=c(1,2))
boxplot(Area~fcrew,data=pollack1,notch=T,xlab="Crew",ylab="Area")
boxplot(Soak~fcrew,data=pollack1,notch=T,xlab="Crew",ylab="Soak time")

#6.2.2# Two possible offsets

pollack1$offs1<-pollack1$Area*(pollack1$Soak/60)
hist(log(pollack1$offs1),xlab="Offset.1",main="")

# Offset (num pollack per m2 per h fishing). Considero que al calcular el area
# de la red en funcion del tamaño y numero de paños que tiene cada arte
# las artes se igualan y el esfuerzo en cuanto a diferencia entre artes
# queda ya contemplado, pero no la eficiencia de cada arte. Para ello
# podemos incluir arte como factor. En cuanto a GRT y crew, al tratarse de
# artes pasivas no deberian influir mucho, asi que se pueden incluir como
# explicativas en el modelo (factores o continuas) 

# pollack1$offs2<-pollack1$Area
# hist(log(pollack1$offs2),xlab="Offset.2",main="")

# Offset secundario. Consideramos solo el tamaño de la red y usaríamos el tiempo
# de calado como variable explicativa

#6.3# Calculation of nighttime fishing

library(lubridate) # Necesaria para obtener horas y minutos de objetos POSIXct

clock.M <- function (t) {hour(t)*60 + minute(t)} # Calcula el minuto entre 1 y 1440 de un dia

pho<-read.csv2(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Sunrise_set.csv",header=T,dec=".",sep=",")

head(pho)
pho<-pho[,c(1:3,7:12)] # Nos quedamos con las columnas que interesan

#6.3.1# Convertir a POSIXct y obtener informacion sobre las fechas y horas

head(pollack1)

pollack1$Deployment<-as.POSIXct(pollack1$Deployment,format="%Y-%m-%d %H:%M:%S") 
pollack1$Retrieval<-as.POSIXct(pollack1$Retrieval,format="%Y-%m-%d %H:%M:%S")
pho$Rise<-as.POSIXct(pho$Rise,format="%d/%m/%Y %H:%M:%S")
pho$Set<-as.POSIXct(pho$Set,format="%d/%m/%Y %H:%M:%S")

pollack1$yLarg<-year(pollack1$Deployment) # Obtener año de largada
pollack1$mLarg<-month(pollack1$Deployment) # Obtener mes de largada
pollack1$dLarg<-day(pollack1$Deployment) # Obtener dia de largada
pollack1$yVir<-year(pollack1$Retrieval) # Obtener año de virada
pollack1$mVir<-month(pollack1$Retrieval) # Obtener mes de virada
pollack1$dVir<-day(pollack1$Retrieval) # Obtener dia de virada

pollack1$hLarg<-clock.M(pollack1$Deployment) # Calcula minuto largada
pollack1$hVir<-clock.M(pollack1$Retrieval) # Calcula minuto virada
pho$hRise<-clock.M(pho$Rise) # Calcula minuto Sunrise
pho$hSet<-clock.M(pho$Set) # Calcula minuto Sunset

#6.3.2# Obtener informacion cruzando los datos de pesca con fotoperiodo

# Obtener dia unico para la largada de entre todos los dias posibles desde
# el 1/1/1999 (dia 1) hasta el 31/12/2013 (dia 5479) 

pollack1$tLarg<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$tLarg[i]<-pho$Timer[which(pollack1$yLarg[i]==pho$Year & pollack1$mLarg[i]==pho$Month & pollack1$dLarg[i]==pho$Day)]} 

# Obtener dia unico para la virada de entre todos los dias posibles desde
# el 1/1/1999 (dia 1) hasta el 31/12/2013 (dia 5479) 

pollack1$tVir<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$tVir[i]<-pho$Timer[which(pollack1$yVir[i]==pho$Year & pollack1$mVir[i]==pho$Month & pollack1$dVir[i]==pho$Day)]}  

# Obtener minuto Sunrise el dia de largada

pollack1$hRiseL<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$hRiseL[i]<-pho$hRise[which(pollack1$tLarg[i]==pho$Timer)]} 

# Obtener minuto Sunset el dia de largada

pollack1$hSetL<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$hSetL[i]<-pho$hSet[which(pollack1$tLarg[i]==pho$Timer)]} 

# Obtener minuto Sunrise el dia de virada

pollack1$hRiseV<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$hRiseV[i]<-pho$hRise[which(pollack1$tVir[i]==pho$Timer)]} 

# Obtener minuto Sunset el dia de virada

pollack1$hSetV<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$hSetV[i]<-pho$hSet[which(pollack1$tVir[i]==pho$Timer)]} 

# Obtener minutos nocturnos transcurridos entre dias en los casos en que
# pasó mas de un dia entre largada y virada

pollack1$minNigTot<-as.vector(rep(NA,dim(pollack1)[1]))

for (i in 1:dim(pollack1)[1]) {pollack1$minNigTot[i]<-if (pollack1$tVir[i]-pollack1$tLarg[i]<=1) {0}
	else {sum(pho$Night[which(pho$Timer==pollack1$tLarg[i]+1) : which(pho$Timer== pollack1$tVir[i])])}
	} 

#6.3.3# Calcular minutos nocturnos

#6.3.3.1# Minutos nocturnos si la largada y la virada tienen lugar el mismo dia (6 combinaciones posibles) 

pollack1$minNig1<-as.vector(rep(NA,dim(pollack1)[1]))

for (i in 1:dim(pollack1)[1]) {pollack1$minNig1[i]<-if (pollack1$tVir[i]-pollack1$tLarg[i]==0) {
		
		# Largada y Virada ocurren antes del Sunrise
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] <= pollack1$hRiseL[i],pollack1$hVir[i]-pollack1$hLarg[i],
		
		# Largada ocurre antes del Sunrise y Virada entre Sunrise y Sunset 
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] > pollack1$hRiseL[i] & pollack1$hVir[i] <= pollack1$hSetL[i],pollack1$hRiseL[i]-pollack1$hLarg[i], 
		
		# Largada ocurre antes del Sunrise y Virada despues del Sunset
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] > pollack1$hSetL[i],(pollack1$hRiseL[i]-pollack1$hLarg[i])+(pollack1$hVir[i]-pollack1$hSetL[i]),
		
		# Largada ocurre entre Sunrise y Sunset y Virada ocurre entre Sunrise y Sunset 
		ifelse(pollack1$hLarg[i] >= pollack1$hRiseL[i] & pollack1$hLarg[i] <= pollack1$hSetL & pollack1$hVir[i] >= pollack1$hRiseL & pollack1$hVir[i] <= pollack1$hSetL[i],0,
		
		# Largada ocurre entre Sunrise y Sunset y Virada ocurre despues del Sunset
		ifelse(pollack1$hLarg[i] >= pollack1$hRiseL[i] & pollack1$hLarg[i] <= pollack1$hSetL & pollack1$hVir[i] > pollack1$hSetL[i],pollack1$hVir[i]-pollack1$hSetL[i],
		
		# Largada y Virada ocurren despues del Sunset
		pollack1$hVir[i]-pollack1$hLarg[i])))))
		
	} else {0}

} 

#6.3.3.2# Minutos nocturnos si la virada tiene lugar al dia siguiente de la largada (9 combinaciones posibles)

pollack1$minNig2<-as.vector(rep(NA,dim(pollack1)[1]))

for (i in 1:dim(pollack1)[1]) {pollack1$minNig2[i]<-if (pollack1$tVir[i]-pollack1$tLarg[i]==1) {
		
		# Largada ocurre antes del Sunrise y Virada ocurre antes del Sunrise del dia siguiente 
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] <= pollack1$hRiseV[i],(pollack1$hRiseL[i]-pollack1$hLarg[i])+((1440-pollack1$hSetL[i])+pollack1$hVir[i]),
		
		# Largada ocurre antes del Sunrise y Virada ocurre entre el Sunrise y el Sunset del dia siguiente
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] > pollack1$hRiseV[i] & pollack1$hVir[i] <= pollack1$hSetV[i],(pollack1$hRiseL[i]-pollack1$hLarg[i])+((1440-pollack1$hSetL[i])+pollack1$hRiseV[i]),
		
		# Largada ocurre antes del Sunrise y Virada ocurre despues del Sunset del dia siguiente
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] > pollack1$hSetV[i],(pollack1$hRiseL[i]-pollack1$hLarg[i])+(1440-pollack1$hSetL[i])+pollack1$hRiseV[i]+(pollack1$hVir[i]-pollack1$hSetV[i]),
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre antes del Sunrise del dia siguiente 
		ifelse(pollack1$hLarg[i] >= pollack1$hRiseL[i] & pollack1$hLarg[i] <= pollack1$hSetL[i] & pollack1$hVir[i] <= pollack1$hRiseV[i],(1440-pollack1$hSetL[i])+pollack1$hVir[i],
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre entre el Sunrise y el Sunset del dia siguiente 
		ifelse(pollack1$hLarg[i] >= pollack1$hRiseL[i] & pollack1$hLarg[i] <= pollack1$hSetL[i] & pollack1$hVir[i] >= pollack1$hRiseV[i] & pollack1$hVir[i] <= pollack1$hSetV[i],(1440-pollack1$hSetL[i])+pollack1$hRiseV[i],
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre despues del Sunset del dia siguiente
		ifelse(pollack1$hLarg[i] >= pollack1$hRiseL[i] & pollack1$hLarg[i] <= pollack1$hSetL[i] & pollack1$hVir[i] > pollack1$hSetV[i],(1440-pollack1$hSetL[i])+pollack1$hRiseV[i]+(pollack1$hVir[i]-pollack1$hSetV[i]),
		
		# Largada ocurre despues del Sunset y virada ocurre antes del Sunrise del dia siguiente
		ifelse(pollack1$hLarg[i] > pollack1$hSetL[i] & pollack1$hVir[i] <= pollack1$hRiseV[i],(1440-pollack1$hLarg[i])+pollack1$hVir[i],
		
		# Largada ocurre despues del Sunset y virada ocurre entre el Sunrise y el Sunset del dia siguiente
		ifelse(pollack1$hLarg[i] > pollack1$hSetL[i] & pollack1$hVir[i] >= pollack1$hRiseV[i] & pollack1$hVir[i] <= pollack1$hSetV[i],(1440-pollack1$hLarg[i])+pollack1$hRiseV[i],
		
		# Largada ocurre despues del Sunset y virada ocurre despues del Sunset del dia siguiente
		(1440-pollack1$hLarg[i])+pollack1$hRiseV[i]+(pollack1$hVir[i]-pollack1$hSetV[i])))))))))
		
	} else {0}

}

#6.3.3.3# Minutos nocturnos si entre la largada y la virada pasa mas de un dia (las mismas 9 combinaciones posibles)

pollack1$minNig3<-as.vector(rep(NA,dim(pollack1)[1]))

for (i in 1:dim(pollack1)[1]) {pollack1$minNig3[i]<-if (pollack1$tVir[i]-pollack1$tLarg[i] > 1) {
		
		# Largada ocurre antes del Sunrise y Virada ocurre antes del Sunrise varios dias despues
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] <= pollack1$hRiseV[i],(pollack1$hRiseL[i]-pollack1$hLarg[i])+(pollack1$minNigTot[i]-(pollack1$hRiseV[i]-pollack1$hVir[i])),
		
		# Largada ocurre antes del Sunrise y Virada ocurre entre el Sunrise y el Sunset varios dias despues
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] > pollack1$hRiseV[i] & pollack1$hVir[i] <= pollack1$hSetV[i],(pollack1$hRiseL[i]-pollack1$hLarg[i])+pollack1$minNigTot[i],
		
		# Largada ocurre antes del Sunrise y Virada ocurre despues del Sunset varios dias despues
		ifelse(pollack1$hLarg[i] < pollack1$hRiseL[i] & pollack1$hVir[i] > pollack1$hSetV[i],(pollack1$hRiseL[i]-pollack1$hLarg[i])+pollack1$minNigTot[i]+(pollack1$hVir[i]-pollack1$hSetV[i]),
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre antes del Sunrise varios dias despues
		ifelse(pollack1$hLarg[i] >= pollack1$hRiseL[i] & pollack1$hLarg[i] <= pollack1$hSetL[i] & pollack1$hVir[i] <= pollack1$hRiseV[i],pollack1$minNigTot[i]-(pollack1$hRiseV[i]-pollack1$hVir[i]),
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre entre el Sunrise y el Sunset varios dias despues
		ifelse(pollack1$hLarg[i] >= pollack1$hRiseL[i] & pollack1$hLarg[i] <= pollack1$hSetL[i] & pollack1$hVir[i] >= pollack1$hRiseV[i] & pollack1$hVir[i] <= pollack1$hSetV[i],pollack1$minNigTot[i],
		
		# Largada ocurre entre el Sunrise y el Sunset y Virada ocurre despues del Sunset varios dias despues
		ifelse(pollack1$hLarg[i] >= pollack1$hRiseL[i] & pollack1$hLarg[i] <= pollack1$hSetL[i] & pollack1$hVir[i] > pollack1$hSetV[i],pollack1$minNigTot[i]+(pollack1$hVir[i]-pollack1$hSetV[i]),
		
		# Largada ocurre despues del Sunset y virada ocurre antes del Sunrise varios dias despues
		ifelse(pollack1$hLarg[i] > pollack1$hSetL[i] & pollack1$hVir[i] <= pollack1$hRiseV[i],pollack1$minNigTot[i]-((pollack1$hLarg[i]-pollack1$hSetL[i])-(pollack1$hRiseV[i]-pollack1$hVir[i])),
		
		# Largada ocurre despues del Sunset y virada ocurre entre el Sunrise y el Sunset varios dias despues
		ifelse(pollack1$hLarg[i] > pollack1$hSetL[i] & pollack1$hVir[i] >= pollack1$hRiseV[i] & pollack1$hVir[i] <= pollack1$hSetV[i],pollack1$minNigTot[i]-(pollack1$hLarg[i]-pollack1$hSetL[i]),
		
		# Largada ocurre despues del Sunset y virada ocurre despues del Sunset varios dias despues
		(pollack1$minNigTot[i]-(pollack1$hLarg[i]-pollack1$hSetL[i]))+(pollack1$hVir[i]-pollack1$hSetV[i])))))))))
	
	} else {0}

}

#6.3.4# Nueva variable '% minutos nocturnos'

#6.3.4.1# Check that soak from UTPB is correct

# pollack1$Soaktime<-as.numeric(difftime(pollack1$Retrieval,pollack1$Deployment,units="mins"))

# plot(pollack1$Soak~Soaktime)
# soakTest<-ifelse(pollack1$Soak==Soaktime,0,1)
# sum(soakTest) # Los dos valores son iguales asi que usamos el suyo
# nn<-which(soakTest==1,)
# pollack1[nn,]

#6.3.4.2# Calculate proportion of night time

pollack1$caladoNight<-as.vector(rep(NA,dim(pollack1)[1])) 

for (i in 1:dim(pollack1)[1]) {pollack1$caladoNight[i]<-if (pollack1$tVir[i]==pollack1$tLarg[i]) {
	
	# % minutos nocturnos si largada y virada ocurren el mismo dia
	round(pollack1$minNig1[i]/pollack1$Soak[i],2)
	
	} else {if (pollack1$tVir[i]-pollack1$tLarg[i]==1) {
		
		# % minutos nocturnos si virada ocurre al dia siguiente de largada
		round(pollack1$minNig2[i]/pollack1$Soak[i],2)
		
		} else {
			
			# % minutos nocturnos si virada ocurre varios dias despues de largada
			round(pollack1$minNig3[i]/pollack1$Soak[i],2)
		}
	}
}

summary(pollack1$caladoNight)
hist(pollack1$caladoNight)

#6.3.5# Standarizar hora de largada segun fotoperiodo

#6.3.5.1# Obtener minuto Sunrise el dia siguiente de largada

names(pollack1)

pollack1$hRiseL1<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$hRiseL1[i]<-pho$hRise[which(pollack1$tLarg[i]+1==pho$Timer)]} 

#6.3.5.2# Corregir hora de largada segun fotoperiodo: si la hora de largada ocurre antes del
# sunset restamos hora de largada-hora sunrise (con lo que largadas antes de sunrise
# tendran valores negativos, y largadas despues de sunrise hasta sunset tedran valores
# positivos. Si la largada ocurre despues del sunset los valores seran muy negativos
# por que los contaremos respecto del sunrise del dia siguiente)

pollack1$hLargC<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$hLargC[i]<-if (pollack1$hLarg[i] <= pollack1$hSetL[i]) {
		
		pollack1$hLarg[i] - pollack1$hRiseL[i]
		
		} else {
			
			((1440 - pollack1$hLarg[i]) + pollack1$hRiseL1[i])*-1
			
		}
	}

#6.3.6# Standarizar hora de virada segun fotoperiodo

#6.3.6.1# Obtener minuto Sunrise el dia siguiente de largada

names(pollack1)

pollack1$hSetV1<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$hSetV1[i]<-pho$hSet[which(pollack1$tVir[i]-1==pho$Timer)]} 

#6.3.6.2# Corregir hora de virada segun fotoperiodo: si la hora de virada ocurre antes del
# sunrise restamos hora de virada-hora sunset del dia anterior (con lo que viradas mas
# lejos del sunset tienen valores muy negativos). Si la virada ocurre despues del sunrise
# los valores seran positivos

pollack1$hVirC<-as.vector(rep(NA,dim(pollack1)[1]))
for (i in 1:dim(pollack1)[1]) {pollack1$hVirC[i]<-if (pollack1$hVir[i] <= pollack1$hRiseV[i]) {
		
		((1440 - pollack1$hSetV1[i]) + pollack1$hVir[i])*-1
		
		} else {
			
			pollack1$hVir[i] - pollack1$hRiseV[i]
			
		}
	}

#6.3.7# Uso de factores (vale 0 si la proporcion de calado nocturno es
# mayor que 75%, si esta entre el 25 y 75% vale 1, y si es menor o igual que
# el 25%, esto es, el lance es practicamente todo diurno, vale 2)

# pollack1$Period<-factor(ifelse(pollack1$caladoNight <= 0.25,2,
#	ifelse(pollack1$caladoNight > 0.25 & pollack1$caladoNight <= 0.75,1,0)))

#6.4# Distribution of all potential explanatory variables

d1<-ggplot(data=pollack1,aes(x=Crew))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("Crew")+
	scale_y_continuous("Density")+
	scale_fill_manual(name="Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d2<-ggplot(data=pollack1,aes(x=GRT))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("GRT")+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d3<-ggplot(data=pollack1,aes(x=Area))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous(expression(paste("Gear area"," ","(",m^2,")")))+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d4<-ggplot(data=pollack1,aes(x=Soak/1440))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("Soak time (days)")+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d5<-ggplot(data=pollack1,aes(x=hLarg/60))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("Deployment (Local time)")+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d6<-ggplot(data=pollack1,aes(x=hVir/60))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("Retrieval (Local time)")+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d7<-ggplot(data=pollack1,aes(x=Julian))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("Day of the Year")+
	scale_y_continuous("Density")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d8<-ggplot(data=pollack1,aes(x=Lat))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("Latitude (ºN)")+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d9<-ggplot(data=pollack1,aes(x=Lon))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("Longitude (ºW)")+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d10<-ggplot(data=pollack1,aes(x=Depth))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("Depth (m)")+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"))+
	theme(legend.position="none")
d11<-ggplot(data=pollack1,aes(x=sstM))+
	geom_density(aes(fill=fgear),alpha=0.5)+
	scale_x_continuous("SST (ºC)")+
	scale_y_continuous("")+
	scale_fill_manual("Gear",values=c("#3182bd","orange"),
		labels=c("Miño","Veta"))+
	guides(fill=guide_legend(override.aes=list(colour=NULL)))+
	theme(legend.position=c(1.8,0.5),
		  legend.text=element_text(size=12))

pdf(file="/Users/jaimeoterovillar/Desktop/explPollack.pdf",width=20,height=8)

multiplot(d1,d7,d2,d8,d3,d9,d4,d10,d5,d11,d6,cols=6)

dev.off()

ggplot(data=pollack1,aes(x=caladoNight))+
	geom_density(fill="gray40",alpha=0.5,col="gray50")+
	scale_x_continuous("Nighttime soak (%)")+
	scale_y_continuous("Density")

ggplot(data=pollack1,aes(x=hLargC))+
	geom_density(fill="gray40",alpha=0.5,col="gray50")+
	scale_x_continuous("Corrected deployment time")+
	scale_y_continuous("Density")

ggplot(data=pollack1,aes(x=hVirC))+
	geom_density(fill="gray40",alpha=0.5,col="gray50")+
	scale_x_continuous("Corrected retrieval time")+
	scale_y_continuous("Density")

summary(pollack1[,c("GRT","Area","hLarg","hVir","Depth","Soak")])


# ----------------------------- #
#7# Exploratory steps before modelling

names(pollack1)
pollack1<-pollack1[,-c(31,34:51,53,55)]
names(pollack1)
summary(pollack1) # No NAs

#7.1 Simple relationships

par(mfcol=c(2,5))
plot(log(Ntot+1)~log(GRT),data=pollack1)
boxplot(log(Ntot+1)~fcrew,data=pollack1)
boxplot(log(Ntot+1)~fyear,data=pollack1)
plot(log(Ntot+1)~Julian,data=pollack1)
plot(log(Ntot+1)~Depth,data=pollack1)
plot(log(Ntot+1)~sstM,data=pollack1)
boxplot(log(Ntot+1)~ZoneA,data=pollack1)
boxplot(log(Ntot+1)~ZoneO,data=pollack1)
plot(log(Ntot+1)~Soak,data=pollack1)

par(mfrow=c(1,3))
hist(pollack1$GRT)
hist(pollack1$Depth)
hist(pollack1$sstM)

#7.2# Collinearity

head(pollack1)

# pollack1$fZoneI<-factor(ifelse(pollack1$ZoneO==1,1,2))

pollack1$lGRT<-log(pollack1$GRT)
pollack1$lDepth<-log(pollack1$Depth)
#pollack1$JulianC<-pollack1$Julian-183

vifs1<-c("lGRT","Crew","Year","Julian","lDepth","sstM","caladoNight")
corvif(pollack1[,vifs1]) # High correlation between Crew & GRT. Keep GRT because is more informative

#7.3# Further exploration of response

names(pollack1)

hist(log(pollack1$Wtot))
quartz()
plot(log(Ntot+1)~log(Wtot),data=pollack1,ylab="Total number",xlab="Total biomass")

size1<-lm(I(Wtot/Ntot)~fgear,data=pollack1)
summary(size1)
plot(allEffects(size1))

pollack1$cpue<-(pollack1$Ntot*200)/pollack1$offs1 # Nominal cpue standardized at 200 m2 hour

dotchart(pollack1$cpue) # Maybe an outlier? La utpb dice que es real!!
which(pollack1$cpue > 5)
pollack1[5466,]

#7.4# Eploration of cpue relationships

plot(cpue~lDepth,data=pollack1[pollack1$fgear=="VETAS",])
plot(cpue~lDepth,data=pollack1[pollack1$fgear=="MINOS",])
plot(cpue~caladoNight,data=pollack1[pollack1$fgear=="VETAS",])
plot(cpue~caladoNight,data=pollack1[pollack1$fgear=="MINOS",])
boxplot(cpue~fZoneO,data=pollack1)


# ----------------------------- #
#8# Modelling of catches using non-mixed models

#8.1# GAM Poisson modelling 

poiss0<-gam(Ntot~offset(log(offs1))+
				 fgear+fZoneO+Seafloor+
				 s(lGRT,k=3)+
			 	 s(Year,k=3)+s(Julian,k=6,bs="cc")+
			 	 s(lDepth,k=3)+
			 	 s(sstM,k=3)+
			 	 s(caladoNight,k=3),
			 	 family=poisson,data=pollack1)

summary(poiss0)
plot(poiss0,pages=1,all.terms=T,scale=0)

#8.2# GLM Poisson modelling

poiss1<-glm(Ntot~offset(log(offs1))+
				 fgear+
				 lGRT+
			 	 fyear+poly(Julian,2)+
			 	 lDepth+
			 	 sstM+
			 	 caladoNight+
			 	 fZoneO+Seafloor+
			 	 fgear:lDepth+fgear:lGRT+fgear:caladoNight+
			 	 fZoneO:fyear,
			 	 family=poisson,data=pollack1)

summary(poiss1)
sum(residuals(poiss1,type="pearson")^2)/poiss1$df.resid # Overdispersion
plot(allEffects(poiss1))

#8.3# Negative Binomial modelling

negbin1<-glm.nb(Ntot~offset(log(offs1))+
					 fgear+
					 lGRT+
			 		 fyear+poly(Julian,2)+
			 		 lDepth+
			 		 sstM+
			 		 caladoNight+
			 		 fZoneO+Seafloor+
			 		 fgear:lDepth+fgear:lGRT+fgear:caladoNight+
			 	 	 fZoneO:fyear,
			 		 data=pollack1)

summary(negbin1)
sum(residuals(negbin1,type="pearson")^2)/negbin1$df.resid
hist(residuals(negbin1,type="deviance"))
plot(residuals(negbin1,type="deviance")~fitted(negbin1))
plot(allEffects(negbin1))

#8.4# Hurdle modelling (assumption: two step model)

hurdle1<-hurdle(Ntot~offset(log(offs1))+
					 fgear+
					 lGRT+
					 fyear+poly(Julian,2)+
					 lDepth+
					 sstM+
					 caladoNight+
					 fZoneO+Seafloor+
					 fgear:lDepth+fgear:lGRT+fgear:caladoNight+
			 	 	 fZoneO:fyear |
					 lDepth,
					 dist="poisson",link="logit",data=pollack1)

summary(hurdle1)

hurdle2<-hurdle(Ntot~offset(log(offs1))+
					 fgear+
					 lGRT+
					 fyear+poly(Julian,2)+
					 lDepth+
					 sstM+
					 caladoNight+
					 fZoneO+Seafloor+
					 fgear:lDepth+fgear:lGRT+fgear:caladoNight+
			 	 	 fZoneO:fyear |
					 lDepth,
					 dist="negbin",link="logit",data=pollack1)

summary(hurdle2)

#8.5# ZIP modelling (assumption: excess zeroes related to depth)

zipoiss1<-zeroinfl(Ntot~offset(log(offs1))+
				   		fgear+
				   		lGRT+
				   		fyear+poly(Julian,2)+
				   		lDepth+
				   		sstM+
				   		caladoNight+
				   		fZoneO+Seafloor+
				   		fgear:lDepth+fgear:lGRT+fgear:caladoNight+
			 	 		fZoneO:fyear |
				   		lDepth,
				   		dist="poisson",link="logit",data=pollack1)

summary(zipoiss1)
sum(residuals(zipoiss1,type="pearson")^2)/(nrow(pollack1)-(length(coef(zipoiss1)))) # Overdispersed

#8.6# ZINB

zinb1<-zeroinfl(Ntot~offset(log(offs1))+
				   	 fgear+
				   	 lGRT+
				   	 fyear+poly(Julian,2)+
				   	 lDepth+
				   	 sstM+
				   	 caladoNight+
				   	 fZoneO+Seafloor+
				   	 fgear:lDepth+fgear:lGRT+fgear:caladoNight+
			 	 	 fZoneO:fyear |
				   	 lDepth,
				   	 dist="negbin",link="logit",data=pollack1)

summary(zinb1)
sum(residuals(zinb1,type="pearson")^2)/(nrow(pollack1)-(length(coef(zinb1))+1)) # Slightly overdispersed

#8.7# Selection of predictors

zinb2<-zeroinfl(Ntot~offset(log(offs1))+
				   	 fgear+
				   	 lGRT+
				   	 fyear+poly(Julian,2)+
				   	 lDepth+
				   	 caladoNight+
				   	 fZoneO+Seafloor+
				   	 fgear:lDepth+fgear:lGRT+fgear:caladoNight+
			 	 	 fZoneO:fyear |
				   	 lDepth,
				   	 dist="negbin",link="logit",data=pollack1)

lrtest(zinb1,zinb2) # Removing of sstM

zinb1<-zinb2 # Rename final model
summary(zinb1)
sum(residuals(zinb1,type="pearson")^2)/(nrow(pollack1)-(length(coef(zinb1))+1))

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

par(mfcol=c(2,4))
boxplot(residZINB~fgear,data=pollack1)
boxplot(residZINB~fyear,data=pollack1)
boxplot(residZINB~fZoneO,data=pollack1)
boxplot(residZINB~Seafloor,data=pollack1)
plot(residZINB~lGRT,data=pollack1)
plot(residZINB~Julian,data=pollack1)
plot(residZINB~lDepth,data=pollack1)
plot(residZINB~caladoNight,data=pollack1)

#8.8.3# Spatial patterns

ggplot(data=pollack1,aes(x=Lon,y=Lat))+
	geom_point(aes(colour=residZINB))+
	scale_colour_gradient2(low="blue",mid="white",high="red",midpoint=15)

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
	scale_y_continuous("Sample variogram",limits=c(0,3))+
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

#8.9.1# Update competing models (remove sstM)

poiss1<-update(poiss1,.~.-sstM)
negbin1<-update(negbin1,.~.-sstM)
hurdle1<-hurdle(Ntot~offset(log(offs1))+fgear+lGRT+fyear+poly(Julian,2)+
					 lDepth+caladoNight+fZoneO+Seafloor+fgear:lDepth+
					 fgear:lGRT+fgear:caladoNight+fZoneO:fyear |
					 lDepth,dist="poisson",link="logit",data=pollack1)
hurdle2<-hurdle(Ntot~offset(log(offs1))+fgear+lGRT+fyear+poly(Julian,2)+
					 lDepth+caladoNight+fZoneO+Seafloor+fgear:lDepth+
					 fgear:lGRT+fgear:caladoNight+fZoneO:fyear |
					 lDepth,dist="negbin",link="logit",data=pollack1)
zipoiss1<-zeroinfl(Ntot~offset(log(offs1))+fgear+lGRT+fyear+poly(Julian,2)+
					    lDepth+caladoNight+fZoneO+Seafloor+fgear:lDepth+
					    fgear:lGRT+fgear:caladoNight+fZoneO:fyear |
					    lDepth,dist="poisson",link="logit",data=pollack1)					 
#8.9.2# Comparison of models (just some tools, more could be added, see Potts & Elith 2006)

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

#8.9.3# Vuong test to compare NB vs ZINB

vuong(negbin1,zinb1) # zinb model is sligthly better than negbin model

#8.9.4# Predicting probabilities (Poisson and Negative Binomial)

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
# NOT YET FINISHED!!!

exp(coef(zinb1)[c(57,58)]) # Zero-inflated part
# Baseline odds of no catching a pollack is 9.327983e-06; the odds is
# increased by one unit increase in Depth by 1.015645e+01 (for a one-unit
# change in Depth there is a ~10% increase in the odds of non-catching a pollack)

exp(coef(zinb1)[c(1,20)]) # Count part
# The baseline number of pollack catch is 1.442720e-06. A unit increase in
# Depth decrease 0.48 times the expected catch of pollack (for a one-unit
# change in Depth there is a ~52% decrease in the expected catch of pollack)
# Fishing in soft seafloor decreases the expected catch by ~73% 100*(2.709858e-01-1)...


# ----------------------------- #
#9# Bootstrapping the optimal model coefficients (zero & count parts)

#9.1# Function (add starting values to the model if needed!!)

#dput(round(coef(zinb1,"count"),4))
#dput(round(coef(zinb1,"zero"),4))

boot.zinb<-function (data,indices) {
	
	data<-data[indices,] 
	try(mod<-zeroinfl(Ntot~offset(log(offs1))+fgear+lGRT+
						   fyear+poly(Julian,2)+lDepth+
						   caladoNight+fZoneO+Seafloor+fgear:lDepth+
						   fgear:lGRT+fgear:caladoNight+fZoneO:fyear |
						   lDepth,data=data,dist="negbin",link="logit"))
						   			  
	if (exists("mod")) { 
						   			 
	as.vector(t(do.call(rbind,coef(summary(mod)))[,1])) # Select 'estimate' column
	
	} else {rep(NA,times=59)} # If the above model crashes this fills in the gaps
    						  # with 59 NA (number of parameters of the model)
    						  # and the algorithm continues
		
}

#9.2# Coefficients (obtain CI of estimates excluding SE and theta)

#library(snow)

RR<-1000 # Number of resamples (reduce this number for testing the code)
zinb.boot.out<-boot(data=pollack1,statistic=boot.zinb,R=RR) # If no problems with snow add ,parallel="snow",ncpus=4
# zinb.boot.out  # Basic output
# plot(zinb.boot.out,index=1) # Example of histogram and qqplot for a given component

zinb.boot.out2<-as.data.frame(zinb.boot.out$t[,c(1:56,58:59)]) # Coefficients of interest from the boot object matrix
colnames(zinb.boot.out2)<-names(coef(zinb1))
write.table(x=zinb.boot.out2,file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootCoefs.txt",row.names=F)

parmsPERC<-t(sapply(c(1:56,58:59),function (i) {
	out<-boot.ci(zinb.boot.out,index=i,type=c("perc"))
	with(out,c(Estimate=t0,pLow=percent[4],pUpp=percent[5]))
})) # Obtain intervals calculated using the bootstrap percentile method

row.names(parmsPERC)<-names(coef(zinb1))
write.table(x=as.data.frame(parmsPERC),file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootPERC.txt",row.names=T)

parmsBCA<-t(sapply(c(1:56,58:59),function (i) {
	out<-boot.ci(zinb.boot.out,index=i,type=c("bca"))
	with(out,c(Estimate=t0,bcaLow=bca[4],bcaUpp=bca[5]))
})) # Obtain intervals calculated using the bootstrap bias-corrected accelerated method

row.names(parmsBCA)<-names(coef(zinb1))
write.table(x=as.data.frame(parmsBCA),file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootBCA.txt",row.names=T)

#9.3# Simple histograms of all components

zinb.boot.out2<-read.table(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootCoefs.txt",header=T,sep=" ",dec=".")

pdf(file="/Users/jaimeoterovillar/Desktop/bootCoefs1.pdf",width=22,height=14)

par(mfrow=c(6,10))
par(mar=c(4,4,4,4))
par(oma=c(1,5,0,0))
for (i in 1:58) {
  hist(zinb.boot.out2[,i],col="light blue",main="",ylab="",
      xlab=names(coef(zinb1))[i],cex.axis=1.2,cex.lab=1.4)
  abline(v=coef(zinb1)[i],lwd=3,col="orange")
  }

mtext("Frequency",side=2,cex=2,outer=T,line=2)

dev.off()

#9.4# Same figure using ggplot of all components (NOT FINISHED!!)

# zinb.boot.out3<-stack(zinb.boot.out2)

# dat.coefs<-as.data.frame(coef(zinb1))
# colnames(dat.coefs)<-"intcpt"
# dat.coefs$ind<-as.character(names(coef(zinb1)))

# pdf(file="/Users/jaimeoterovillar/Desktop/bootCoefs2.pdf",width=20,height=10)

# ggplot(data=zinb.boot.out3,aes(x=values))+
	# geom_histogram(binwidth=0.5,aes(fill=..count..))+
	# geom_vline(data=dat.coefs,aes(xintercept=intcpt,colour="orange"),size=1)+
	# facet_wrap(~ind,scales="free",nrow=6,ncol=10)+
	# scale_y_continuous("Count",limits=c(0,10))+
	# scale_x_continuous("Values")+
	# theme(legend.position=c(1,0),
		  # legend.justification=c(1.4,0),
		  # legend.text=element_text(size=10))+
	# scale_fill_gradient("Count",
		# guide=guide_legend(direction="horizontal",title.position="top",
							# title.hjust=0.5,title.vjust=0.5,
							# label.position="bottom",
							# label.hjust=0.5,label.vjust=0.5,
							# keywidth=1.5,keyheight=1.5))
	

# dev.off()


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
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(pollack1$lDepth),
					   			  max(pollack1$lDepth),length=100),
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200)) # Estandarizacion a num por 200 m2 por h (200 m2 es el tamaño legal aproximado de un paño: 50 m de largo por 4 de alto)

prob0<-data.frame(l1,seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100))
colnames(prob0)<-c("l1","DepthSeq")

rugsDepth<-data.frame(rugsDepth=unique(pollack1$lDepth))

pdf(file="/Users/jaimeoterovillar/Desktop/ZIzero.pdf")

ggplot(data=prob0,aes(x=DepthSeq,y=l1))+
	geom_line(colour="black",lwd=1)+
	scale_y_continuous("Probability of false zeros",limits=c(-0.02,1))+
	scale_x_continuous("ln-Depth (m)",limits=c(0,6))+
	geom_segment(data=rugsDepth,aes(x=rugsDepth,xend=rugsDepth,
		y=-0.02,yend=-0.01),stat="identity",lwd=0.1,col="gray50")

dev.off()

#10.2# Plotting count (NB) curves from ZINB (this coding could be shortend.
# However, I wrote it this way to keep track on what I'm doing!!)

rugsGRT<-data.frame(rugsGRT=unique(pollack1$lGRT))
rugsDoY<-data.frame(rugsDoY=unique(pollack1$Julian))
rugsNight<-data.frame(rugsNight=unique(pollack1$caladoNight))

#10.2.1# Continuous variables

z1<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=seq(from=1,to=366,by=1), # Julian
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDOY<-data.frame(z1,seq(from=1,to=366,by=1))
colnames(probDOY)<-c("Pred","DOYSeq")

z2<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS", # VETAS
					   lGRT=seq(min(pollack1$lGRT,na.rm=T),
					   			max(pollack1$lGRT,na.rm=T),length=100), # lGRT
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probGRT.V<-data.frame(z2,seq(min(pollack1$lGRT),max(pollack1$lGRT),length=100))
colnames(probGRT.V)<-c("Pred","GRTSeq")

z3<-predict(zinb1,type="count",
	newdata=data.frame(fgear="MINOS", # MIÑOS
					   lGRT=seq(min(pollack1$lGRT,na.rm=T),
					   			max(pollack1$lGRT,na.rm=T),length=100), # lGRT
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probGRT.M<-data.frame(z3,seq(min(pollack1$lGRT),max(pollack1$lGRT),length=100))
colnames(probGRT.M)<-c("Pred","GRTSeq")

z4<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS", # VETAS
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(pollack1$lDepth,na.rm=T),
					   			  max(pollack1$lDepth,na.rm=T),length=100), # lDepth
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDepth.V<-data.frame(z4,seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100))
colnames(probDepth.V)<-c("Pred","DepthSeq")

z5<-predict(zinb1,type="count",
	newdata=data.frame(fgear="MINOS", # MIÑOS
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(pollack1$lDepth,na.rm=T),
					   			  max(pollack1$lDepth,na.rm=T),length=100), # lDepth
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDepth.M<-data.frame(z5,seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100))
colnames(probDepth.M)<-c("Pred","DepthSeq")

z6<-predict(zinb1,type="count",
	newdata=data.frame(fgear="VETAS", # VETAS
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=seq(min(pollack1$caladoNight,na.rm=T),
					   			   max(pollack1$caladoNight,na.rm=T),length=100), # caladoNight
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probNight.V<-data.frame(z6,seq(min(pollack1$caladoNight),max(pollack1$caladoNight),length=100))
colnames(probNight.V)<-c("Pred","NightSeq")

z7<-predict(zinb1,type="count",
	newdata=data.frame(fgear="MINOS", # MIÑOS
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=seq(min(pollack1$caladoNight,na.rm=T),
					   			   max(pollack1$caladoNight,na.rm=T),length=100), # caladoNight
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probNight.M<-data.frame(z7,seq(min(pollack1$caladoNight),max(pollack1$caladoNight),length=100))
colnames(probNight.M)<-c("Pred","NightSeq")

f1<-ggplot(data=probDOY,aes(x=DOYSeq,y=Pred))+
	geom_line(colour="black",lwd=0.8)+
	scale_y_continuous("Standardized Abundance",limits=c(-0.01,0.3))+
	scale_x_continuous("Day of the Year",limits=c(1,366),breaks=c(1,100,200,300))+
	geom_segment(data=rugsDoY,aes(x=rugsDoY,xend=rugsDoY,
		y=-0.01,yend=0),stat="identity",lwd=0.1,col="gray50")

f2<-ggplot(data=probGRT.V,aes(x=GRTSeq,y=Pred))+
	geom_line(colour="black",lwd=0.8)+
	scale_y_continuous("Standardized Abundance",limits=c(-0.02,0.8))+
	scale_x_continuous("ln-GRT",limits=c(-1,4))+
	geom_segment(data=rugsGRT,aes(x=rugsGRT,xend=rugsGRT,
		y=-0.02,yend=0),stat="identity",lwd=0.1,col="gray50")

f3<-ggplot(data=probGRT.M,aes(x=GRTSeq,y=Pred))+
	geom_line(colour="black",lwd=0.8,linetype="twodash")+
	scale_y_continuous("",limits=c(-0.0001,0.005))+
	scale_x_continuous("ln-GRT",limits=c(-1,4))+
	geom_segment(data=rugsGRT,aes(x=rugsGRT,xend=rugsGRT,
		y=-0.0001,yend=0),stat="identity",lwd=0.1,col="gray50")

f4<-ggplot(data=probDepth.V,aes(x=DepthSeq,y=Pred))+
	geom_segment(data=rugsDepth,aes(x=rugsDepth,xend=rugsDepth,
		y=-0.5,yend=0),stat="identity",lwd=0.1,col="gray50")+
	geom_line(colour="black",lwd=0.8)+
	scale_y_continuous("Standardized Abundance",limits=c(-0.5,22))+
	scale_x_continuous("ln-Depth (m)",limits=c(0,6))

f5<-ggplot(data=probDepth.M,aes(x=DepthSeq,y=Pred))+
	geom_line(colour="black",lwd=0.8,linetype="twodash")+
	scale_y_continuous("",limits=c(-0.0003,0.012))+
	scale_x_continuous("ln-Depth (m)",limits=c(0,6))+
	geom_segment(data=rugsDepth,aes(x=rugsDepth,xend=rugsDepth,
		y=-0.0003,yend=0),stat="identity",lwd=0.1,col="gray50")

f6<-ggplot(data=probNight.V,aes(x=NightSeq,y=Pred))+
	geom_line(colour="black",lwd=0.8)+
	scale_y_continuous("Standardized Abundance",limits=c(-0.008,0.25))+
	scale_x_continuous("Night soak (%)",limits=c(0,1))+
	geom_segment(data=rugsNight,aes(x=rugsNight,xend=rugsNight,
		y=-0.008,yend=0),stat="identity",lwd=0.1,col="gray50")

f7<-ggplot(data=probNight.M,aes(x=NightSeq,y=Pred))+
	geom_line(colour="black",lwd=0.8,linetype="twodash")+
	scale_y_continuous("",limits=c(-0.0003,0.01))+
	scale_x_continuous("Night soak (%)",limits=c(0,1))+
	geom_segment(data=rugsNight,aes(x=rugsNight,xend=rugsNight,
		y=-0.0003,yend=0),stat="identity",lwd=0.1,col="gray50")

#10.2.2# Categorical variables (only Seafloor. Year is plotted later)

newSeafloor<-data.frame(fgear=rep("VETAS",times=3),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=3),
					fyear=rep("2006",times=3),
					Julian=rep(183,times=3),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=3),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=3),
					fZoneO=rep("1",times=3),Seafloor=c("hard","mixed","soft"),
					offs1=rep(200,times=3))

seafloorP<-predict(zinb1,newdata=newSeafloor,type="count")
probSeafloor<-data.frame(seafloorP=seafloorP,SeafType=factor(c("Hard","Mixed","Soft"),levels=c("Hard","Mixed","Soft")))

f8<-ggplot(data=probSeafloor,aes(x=SeafType,y=seafloorP))+
	geom_bar(stat="identity",colour="black",fill="black")+
	scale_y_continuous("",limits=c(0,0.25))+
	scale_x_discrete("Seafloor type")

#10.3# All plots together

pdf(file="/Users/jaimeoterovillar/Desktop/ZIpredictions.pdf",width=7,height=11)

multiplot(f4,f2,f6,f1,f5,f3,f7,f8,cols=2)

dev.off()


# ----------------------------- #
#11# Plotting ZINB model predicted means for YEAR with Bootstrapping
# of predictions (95% CI) by means of shuffling residuals (Thierry Onkelinx code)

#11.1# Set up the data

Fit<-predict(zinb1,type="response")

Pearson<-residuals(zinb1,type="pearson") # Pearson residuals
VarComp<-residuals(zinb1,type="response")/Pearson # Raw residuals/Pearson residuals

fgear<-pollack1$fgear
lGRT<-pollack1$lGRT
fyear<-pollack1$fyear
Julian<-pollack1$Julian
lDepth<-pollack1$lDepth
caladoNight<-pollack1$caladoNight
fZoneO<-pollack1$fZoneO
Seafloor<-pollack1$Seafloor
offs1<-pollack1$offs1

newTrend.RB<-data.frame(fgear=rep("VETAS",times=15),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=15),
					fyear=c("1999","2000","2001","2002","2003","2004","2005",
						"2006","2007","2008","2009","2010","2011","2012","2013"),
					Julian=rep(183,times=15),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=15),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=15),
					fZoneO=rep("1",times=15),Seafloor=rep("hard",times=15),
					offs1=rep(200,times=15))

abundInd.RB<-predict(zinb1,newdata=newTrend.RB,type="count")
years<-seq(1999,2013,1)
pollackAbund.RB<-data.frame(cbind(abundInd.RB,years))
colnames(pollackAbund.RB)<-c("Index","Year")

newTrend.AA<-data.frame(fgear=rep("VETAS",times=15),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=15),
					fyear=c("1999","2000","2001","2002","2003","2004","2005",
						"2006","2007","2008","2009","2010","2011","2012","2013"),
					Julian=rep(183,times=15),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=15),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=15),
					fZoneO=rep("2",times=15),Seafloor=rep("hard",times=15),
					offs1=rep(200,times=15))

abundInd.AA<-predict(zinb1,newdata=newTrend.AA,type="count")
pollackAbund.AA<-data.frame(cbind(abundInd.AA,years))
colnames(pollackAbund.AA)<-c("Index","Year")

newTrend.CA<-data.frame(fgear=rep("VETAS",times=15),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=15),
					fyear=c("1999","2000","2001","2002","2003","2004","2005",
						"2006","2007","2008","2009","2010","2011","2012","2013"),
					Julian=rep(183,times=15),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=15),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=15),
					fZoneO=rep("3",times=15),Seafloor=rep("hard",times=15),
					offs1=rep(200,times=15))

abundInd.CA<-predict(zinb1,newdata=newTrend.CA,type="count")
pollackAbund.CA<-data.frame(cbind(abundInd.CA,years))
colnames(pollackAbund.CA)<-c("Index","Year")

pollackAbund<-rbind(pollackAbund.RB,pollackAbund.AA,pollackAbund.CA)
pollackAbund$Zone<-rep(c("Rías Baixas","Arco Ártabro","Cantábrico"),each=15)
pollackAbund$Zone<-factor(pollackAbund$Zone,levels=c("Rías Baixas","Arco Ártabro","Cantábrico"))

#11.2# Bootstrap

RR<-1000 # Number of resamples (reduce this number for testing the code)

#11.2.1# Bootstrap Rias Baixas trend

bootstrap.RB<-replicate(n=RR,{ 
	
    yStar<-pmax(round(Fit+sample(Pearson)*VarComp,0),0)
    try(mod<-zeroinfl(yStar~offset(log(offs1))+fgear+lGRT+fyear+poly(Julian,2)+
    						lDepth+caladoNight+fZoneO+Seafloor+
    						fgear:lDepth+fgear:lGRT+fgear:caladoNight+
    						fZoneO:fyear | lDepth,
    						dist="negbin",link="logit"))
    
    if (exists("mod")) {
    	
    	predict(mod,newdata=newTrend.RB,type="response")
    	
    } else {rep(NA,times=15)} # If the above model crashes this fills in the gaps
    						  # with NA and the algorithm continues

})

CI.RB<-t(apply(X=bootstrap.RB,MARGIN=1,FUN=quantile,c(0.025,0.975),na.rm=T))
colnames(CI.RB)<-c("ciLow","ciUpp")
newTrend.RB$fit<-predict(zinb1,newdata=newTrend.RB,type="response")
newTrend.RB<-cbind(newTrend.RB,CI.RB)
newTrend.RB$Year<-seq(1999,2013,1)

write.table(x=as.data.frame(bootstrap.RB),file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootMatrix.RB.txt",row.names=F)
write.table(x=as.data.frame(newTrend.RB),file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootResults.RB.txt",row.names=F)

#11.2.2# Bootstrap Arco Artabro trend

bootstrap.AA<-replicate(n=RR,{ 
	
    yStar<-pmax(round(Fit+sample(Pearson)*VarComp,0),0)
    try(mod<-zeroinfl(yStar~offset(log(offs1))+fgear+lGRT+fyear+poly(Julian,2)+
    						lDepth+caladoNight+fZoneO+Seafloor+
    						fgear:lDepth+fgear:lGRT+fgear:caladoNight+
    						fZoneO:fyear | lDepth,
    						dist="negbin",link="logit"))
    
    if (exists("mod")) {
    	
    	predict(mod,newdata=newTrend.AA,type="response")
    	
    } else {rep(NA,times=15)} 

})

CI.AA<-t(apply(X=bootstrap.AA,MARGIN=1,FUN=quantile,c(0.025,0.975),na.rm=T))
colnames(CI.AA)<-c("ciLow","ciUpp")
newTrend.AA$fit<-predict(zinb1,newdata=newTrend.AA,type="response")
newTrend.AA<-cbind(newTrend.AA,CI.AA)
newTrend.AA$Year<-seq(1999,2013,1)

write.table(x=as.data.frame(bootstrap.AA),file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootMatrix.AA.txt",row.names=F)
write.table(x=as.data.frame(newTrend.AA),file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootResults.AA.txt",row.names=F)

#11.2.3# Bootstrap Cantabrico trend

bootstrap.CA<-replicate(n=RR,{ 
	
    yStar<-pmax(round(Fit+sample(Pearson)*VarComp,0),0)
    try(mod<-zeroinfl(yStar~offset(log(offs1))+fgear+lGRT+fyear+poly(Julian,2)+
    						lDepth+caladoNight+fZoneO+Seafloor+
    						fgear:lDepth+fgear:lGRT+fgear:caladoNight+
    						fZoneO:fyear | lDepth,
    						dist="negbin",link="logit"))
    
    if (exists("mod")) {
    	
    	predict(mod,newdata=newTrend.CA,type="response")
    	
    } else {rep(NA,times=15)}

})

CI.CA<-t(apply(X=bootstrap.CA,MARGIN=1,FUN=quantile,c(0.025,0.975),na.rm=T))
colnames(CI.CA)<-c("ciLow","ciUpp")
newTrend.CA$fit<-predict(zinb1,newdata=newTrend.CA,type="response")
newTrend.CA<-cbind(newTrend.CA,CI.CA)
newTrend.CA$Year<-seq(1999,2013,1)

write.table(x=as.data.frame(bootstrap.CA),file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootMatrix.CA.txt",row.names=F)
write.table(x=as.data.frame(newTrend.CA),file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootResults.CA.txt",row.names=F)

#11.3# Put all standardized abundances together

newTrend.RB<-read.table(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootResults.RB.txt",header=T,sep=" ",dec=".")
newTrend.AA<-read.table(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootResults.AA.txt",header=T,sep=" ",dec=".")
newTrend.CA<-read.table(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/Bootstrap_pollack/bootResults.CA.txt",header=T,sep=" ",dec=".")

abundance.GA<-rbind(newTrend.RB[,10:12],newTrend.AA[,10:12],newTrend.CA[,10:12])
abundance.GA$Year<-rep(1999:2013,times=3)
abundance.GA$Zone<-rep(c("Rías Baixas","Arco Ártabro","Cantábrico"),each=15)
abundance.GA$Zone<-factor(abundance.GA$Zone,levels=c("Rías Baixas","Arco Ártabro","Cantábrico"))
colnames(abundance.GA)<-c("Index","ciLow","ciUpp","Year","Zone")

#11.4# Year-averaged nominal cpue for each oceanographic zone

head(pollack1)

pollack1.RB<-pollack1[pollack1$ZoneO==1,]
cpues.RB.M<-tapply(pollack1.RB$cpue,pollack1.RB$Year,mean,na.rm=T)
#cpues.RB.L<-tapply(pollack1.RB$cpue,pollack1.RB$Year,ci95Low)
#cpues.RB.U<-tapply(pollack1.RB$cpue,pollack1.RB$Year,ci95Upp)
cpues.RB<-data.frame(cbind(cpues.RB.M,rep(NA,15),rep(NA,15),seq(1999,2013,1)))
colnames(cpues.RB)<-c("Index","ciLow","ciUpp","Year")

pollack1.AA<-pollack1[pollack1$ZoneO==2,]
cpues.AA.M<-tapply(pollack1.AA$cpue,pollack1.AA$Year,mean,na.rm=T)
#cpues.AA.L<-tapply(pollack1.AA$cpue,pollack1.AA$Year,ci95Low)
#cpues.AA.U<-tapply(pollack1.AA$cpue,pollack1.AA$Year,ci95Upp)
cpues.AA<-data.frame(cbind(cpues.AA.M,rep(NA,15),rep(NA,15),seq(1999,2013,1)))
colnames(cpues.AA)<-c("Index","ciLow","ciUpp","Year")

pollack1.CA<-pollack1[pollack1$ZoneO==3,]
cpues.CA.M<-tapply(pollack1.CA$cpue,pollack1.CA$Year,mean,na.rm=T)
#cpues.CA.L<-tapply(pollack1.CA$cpue,pollack1.CA$Year,ci95Low)
#cpues.CA.U<-tapply(pollack1.CA$cpue,pollack1.CA$Year,ci95Upp)
cpues.CA<-data.frame(cbind(cpues.CA.M,rep(NA,15),rep(NA,15),seq(1999,2013,1)))
colnames(cpues.CA)<-c("Index","ciLow","ciUpp","Year")

cpues.GA<-rbind(cpues.RB,cpues.AA,cpues.CA)
cpues.GA$Zone<-rep(c("Rías Baixas","Arco Ártabro","Cantábrico"),each=15)
cpues.GA$Zone<-factor(cpues.GA$Zone,levels=c("Rías Baixas","Arco Ártabro","Cantábrico"))
cpues.GA$Year<-rep(1999:2013,times=3)

#11.5# Official landings

catches<-read.csv2(file="/Users/jaimeoterovillar/Documents/Proyectos/ICES/ICES Science Fund/R/2_abadejo landings.csv",header=T,dec=",",sep=",")

catches

catches2<-stack(catches[,c(10:12)])
catches2$Year<-rep(seq(1999,2013,1),times=3)
catches2$ciLow<-rep(NA,times=dim(catches2)[1])
catches2$ciUpp<-rep(NA,times=dim(catches2)[1])

catches.GA<-catches2[,c(1,4,5,3)]
colnames(catches.GA)<-c("Index","ciLow","ciUpp","Year")
catches.GA$Zone<-rep(c("Rías Baixas","Arco Ártabro","Cantábrico"),each=15)
catches.GA$Zone<-factor(catches.GA$Zone,levels=c("Rías Baixas","Arco Ártabro","Cantábrico"))
catches.GA$Year<-rep(1999:2013,times=3)

#11.6# Plot all panels together

# pollackTrends<-rbind(abundance.GA,cpues.GA,catches.GA)
# pollackTrends$Indices<-rep(c("Standardized Abundance","Nominal CPUE","Official landings"),each=45)
# pollackTrends$Indices<-factor(pollackTrends$Indices,levels=c("Standardized Abundance","Nominal CPUE","Official landings"))

fish1<-ggplot(data=abundance.GA,aes(x=Year,y=Index))+
	geom_segment(aes(x=Year,y=Index,xend=Year,yend=ciUpp),col="black",lwd=0.5)+
	geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=Index),col="black",lwd=0.5)+
	geom_line(lwd=0.5,linetype="dotted",col="gray20")+
	geom_point(size=5,col="gray90")+
	geom_point(size=3,col="black")+
	facet_wrap(~Zone,scales="free_y")+
	scale_y_continuous("Standardized Abundance")+
	scale_x_continuous("",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

fish2<-ggplot(data=cpues.GA,aes(x=Year,y=Index))+
	geom_line(lwd=0.5,linetype="dotted",col="gray20")+
	geom_point(size=5,col="gray90")+
	geom_point(size=3,col="black")+
	facet_wrap(~Zone,scales="free_y")+
	scale_y_continuous("Nominal CPUE")+
	scale_x_continuous("",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

fish3<-ggplot(data=catches.GA,aes(x=Year,y=Index))+
	geom_line(lwd=0.5,linetype="dotted",col="gray20")+
	geom_point(size=5,col="gray90")+
	geom_point(size=3,col="black")+
	facet_wrap(~Zone,scales="free_y")+
	scale_y_continuous("Official landings")+
	scale_x_continuous("",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

pdf(file="/Users/jaimeoterovillar/Desktop/PollackTrends.pdf",width=10,height=10)

multiplot(fish1,fish2,fish3,cols=1)

dev.off()

# ggplot(data=pollackTrends,aes(x=Year,y=Index))+
	# geom_segment(aes(x=Year,y=Index,xend=Year,yend=ciUpp),col="orange",lwd=0.5)+
	# geom_segment(aes(x=Year,y=ciLow,xend=Year,yend=Index),col="orange",lwd=0.5)+
	# geom_line(lwd=0.5,linetype="dotted",col="blue")+
	# geom_point(size=5,col="gray90")+
	# geom_point(size=3,col="orange")+
	# facet_wrap(Indices~Zone,scales="free_y")+
	# scale_x_continuous("Year",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

#11.7# Correlation between abundance, cpue and catches

cor(x=cpues.GA$Index[1:15],y=abundance.GA$Index[1:15],use="pairwise.complete.obs",method="spearman")
cor(x=catches.GA$Index[1:15],y=abundance.GA$Index[1:15],use="pairwise.complete.obs",method="spearman")

cor(x=cpues.GA$Index[16:30],y=abundance.GA$Index[16:30],use="pairwise.complete.obs",method="spearman")
cor(x=catches.GA$Index[16:30],y=abundance.GA$Index[16:30],use="pairwise.complete.obs",method="spearman")

cor(x=cpues.GA$Index[31:45],y=abundance.GA$Index[31:45],use="pairwise.complete.obs",method="spearman")
cor(x=catches.GA$Index[31:45],y=abundance.GA$Index[31:45],use="pairwise.complete.obs",method="spearman")

cor(x=abundance.GA$Index[1:15],y=abundance.GA$Index[16:30],use="pairwise.complete.obs",method="spearman")
cor(x=abundance.GA$Index[1:15],y=abundance.GA$Index[31:45],use="pairwise.complete.obs",method="spearman")
cor(x=abundance.GA$Index[16:30],y=abundance.GA$Index[31:45],use="pairwise.complete.obs",method="spearman")

#11.8# Plot of biomass (NOT TO BE FINALLY USED!!!!)

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

#11.9# Plotting maps

#11.9.1# Load contour Galicia

galicia.coast<-read.csv2(file="/Users/jaimeoterovillar/Documents/Imagenes_mapas/Mapas/Otros/Coastline/Galicia_coast.csv",header=F,dec=",")
names(galicia.coast)<-c('lon','lat')

head(galicia.coast)

galicia.coast$lon<-ifelse(galicia.coast$lon > -8.85 & galicia.coast$lat < 42.1,NA,galicia.coast$lon)
galicia.coast$lon<-ifelse(galicia.coast$lon > -8.4 & galicia.coast$lat < 42.5,NA,galicia.coast$lon)

#11.9.2# Load bathymetry

galicia.bathy<-read.table(file="/Users/jaimeoterovillar/Documents/Imagenes_mapas/Mapas/Otros/Bathymetry/Galicia_bathy.txt",sep="",header=F,dec=".",fill=T)
names(galicia.bathy)<-c('lon','lat','depth')

galicia.bathy$depth[galicia.bathy$depth>0]<-NA

head(galicia.bathy)

galicia.bathy.mat<-matrix(galicia.bathy$depth,
	nrow=length(unique(galicia.bathy$lon)),
	ncol=length(unique(galicia.bathy$lat)))[,order(unique(galicia.bathy$lat))]

#11.9.3# Load Europe

europa.coast<-read.csv2(file="/Users/jaimeoterovillar/Documents/Imagenes_mapas/Mapas/Otros/Coastline/NE Atlantic_enlarge.csv",header=F,dec=",")
names(europa.coast)<-c('area','lon','lat')

head(europa.coast)

europa.coast1<-europa.coast[europa.coast$area==1,]
europa.coast2<-subset(europa.coast,europa.coast$area!=1)
areas<-unique(europa.coast2$area)

#11.9.4# Galician map

pdf(file="/Users/jaimeoterovillar/Desktop/galiciaMap.pdf",width=8,height=9)

par(mar=c(5,5,3,3))
plot(galicia.coast,type="n",xlim=c(-9.5,-7.1),ylim=c(41.9,43.9),xlab="Longitude (ºW)",ylab="Latitude (ºN)",cex.lab=1.2,cex.axis=1.1,lwd=1.3)

contour(unique(galicia.bathy$lon),sort(unique(galicia.bathy$lat)),
	galicia.bathy.mat,levels=-seq(0,500,by=50),
	labcex=0.4,col='gray80',add=T)

points(Lat~Lon,data=pollack1[pollack1$ZoneO==1 & pollack1$Gear=="VETAS",],
	pch=1,col="lightblue",cex=0.5) # Hauls in Rias Baixas
points(Lat~Lon,data=pollack1[pollack1$ZoneO==1 & pollack1$Gear=="MINOS",],
	pch=3,col="lightblue",cex=0.5)
points(Lat~Lon,data=pollack1[pollack1$ZoneO==2 & pollack1$Gear=="VETAS",],
	pch=1,col="darkgreen",cex=0.5) # Hauls in Arco Artabro
points(Lat~Lon,data=pollack1[pollack1$ZoneO==2 & pollack1$Gear=="MINOS",],
	pch=3,col="darkgreen",cex=0.5)
points(Lat~Lon,data=pollack1[pollack1$ZoneO==3 & pollack1$Gear=="VETAS",],
	pch=1,col="orange",cex=0.5) # Hauls in Cantabrico
points(Lat~Lon,data=pollack1[pollack1$ZoneO==3 & pollack1$Gear=="MINOS",],
	pch=3,col="orange",cex=0.5) # Hauls in Cantabrico

points(galicia.coast,type="l") # Plot Map on top

points(-8.64639,42.435,cex=2,pch=16,col="red") # Photoperiod location (Pontevedra)

#points(-8.875,42.125,cex=2,pch=4,col="black") # centre of SST cells
#points(-9.125,42.375,cex=2,pch=4,col="black")
#points(-8.875,42.375,cex=2,pch=4,col="black")
#points(-9.125,42.625,cex=2,pch=4,col="black")
#points(-9.125,42.875,cex=2,pch=4,col="black")
#points(-9.125,43.375,cex=2,pch=4,col="black")
#points(-8.375,43.375,cex=2,pch=4,col="black")
#points(-8.125,43.625,cex=2,pch=4,col="black")
#points(-7.375,43.875,cex=2,pch=4,col="black")

legend(-9.55,43.95,legend=c("Cantábrico","Arco Ártabro","Rías Baixas"),bty="n",cex=0.8,pch=16,col=c("orange","darkgreen","lightblue"),title="Oceanographic Zone")

par(new=T)
par(mar=c(6,25,20,4))
plot(europa.coast1$lon,europa.coast1$lat,type="l",xlim=c(-15,5),ylim=c(35,60),xlab="",ylab="",xaxt="n",yaxt="n",bg='transparent',lwd=1.3)

for (i in 1:length(areas)) {
	par(mar=c(6,25,20,4),new=T)
	plot(europa.coast2$lon[europa.coast2$area==areas[i]],
	europa.coast2$lat[europa.coast2$area==areas[i]],
	type="l",xlab="",ylab="",xaxt="n",yaxt="n",
	axes=F,xlim=c(-15,5),ylim=c(35,60),lwd=1.3)
	}

# points(-11,43,cex=1.5,pch=16,col="black") # Upwelling cell
rect(-7,42,-9.5,44,border="red",lwd=2)

dev.off()


# ----------------------------- #
#12# Plotting NB model results

summary(negbin1)

#12.1# Continuous variables

z1<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS",
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=seq(from=1,to=366,by=1), # Julian
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDOY<-data.frame(z1$fit,z1$se.fit,seq(from=1,to=366,by=1))
colnames(probDOY)<-c("z1","z1SE","DOYSeq")

z2<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS", # VETAS
					   lGRT=seq(min(pollack1$lGRT,na.rm=T),
					   			max(pollack1$lGRT,na.rm=T),length=100), # lGRT
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probGRT.V<-data.frame(z2$fit,z2$se.fit,seq(min(pollack1$lGRT),max(pollack1$lGRT),length=100))
colnames(probGRT.V)<-c("Pred","SE","GRTSeq")

z3<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="MINOS", # MIÑOS
					   lGRT=seq(min(pollack1$lGRT,na.rm=T),
					   			max(pollack1$lGRT,na.rm=T),length=100), # lGRT
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probGRT.M<-data.frame(z3$fit,z3$se.fit,seq(min(pollack1$lGRT),max(pollack1$lGRT),length=100))
colnames(probGRT.M)<-c("Pred","SE","GRTSeq")

probGRT<-rbind(probGRT.V,probGRT.M)
probGRT$Gear<-rep(c("Veta","Miño"),each=100)

z4<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS", # VETAS
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(pollack1$lDepth,na.rm=T),
					   			  max(pollack1$lDepth,na.rm=T),length=100), # lDepth
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDepth.V<-data.frame(z4$fit,z4$se.fit,seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100))
colnames(probDepth.V)<-c("Pred","SE","DepthSeq")

z5<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="MINOS", # MIÑOS
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=seq(min(pollack1$lDepth,na.rm=T),
					   			  max(pollack1$lDepth,na.rm=T),length=100), # lDepth
					   caladoNight=mean(pollack1$caladoNight,na.rm=T),
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probDepth.M<-data.frame(z5$fit,z5$se.fit,seq(min(pollack1$lDepth),max(pollack1$lDepth),length=100))
colnames(probDepth.M)<-c("Pred","SE","DepthSeq")

probDepth<-rbind(probDepth.V,probDepth.M)
probDepth$Gear<-rep(c("Veta","Miño"),each=100)

z6<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="VETAS", # VETAS
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=seq(min(pollack1$caladoNight,na.rm=T),
					   			   max(pollack1$caladoNight,na.rm=T),length=100), # caladoNight
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probNight.V<-data.frame(z6$fit,z6$se.fit,seq(min(pollack1$caladoNight),max(pollack1$caladoNight),length=100))
colnames(probNight.V)<-c("Pred","SE","NightSeq")

z7<-predict(negbin1,type="link",se=T,
	newdata=data.frame(fgear="MINOS", # MIÑOS
					   lGRT=mean(pollack1$lGRT,na.rm=T),
					   fyear="2006",
					   Julian=183,
					   lDepth=mean(pollack1$lDepth,na.rm=T),
					   caladoNight=seq(min(pollack1$caladoNight,na.rm=T),
					   			   max(pollack1$caladoNight,na.rm=T),length=100), # caladoNight
					   fZoneO="1",
					   Seafloor="hard",
					   offs1=200))

probNight.M<-data.frame(z7$fit,z7$se.fit,seq(min(pollack1$caladoNight),max(pollack1$caladoNight),length=100))
colnames(probNight.M)<-c("Pred","SE","NightSeq")

probNight<-rbind(probNight.V,probNight.M)
probNight$Gear<-rep(c("Veta","Miño"),each=100)

newSeafloor<-data.frame(fgear=rep("VETAS",times=3),
					lGRT=rep(mean(pollack1$lGRT,na.rm=T),times=3),
					fyear=rep("2006",times=3),
					Julian=rep(183,times=3),
					lDepth=rep(mean(pollack1$lDepth,na.rm=T),times=3),
					caladoNight=rep(mean(pollack1$caladoNight,na.rm=T),times=3),
					fZoneO=rep("1",times=3),Seafloor=c("hard","mixed","soft"),
					offs1=rep(200,times=3))

seafloorP<-predict(negbin1,type="link",se=T,newdata=newSeafloor)
probSeafloor<-data.frame(seafloorP=seafloorP$fit,seafloorP$se.fit,SeafType=factor(c("Hard","Mixed","Soft"),levels=c("Hard","Mixed","Soft")))
colnames(probSeafloor)<-c("Pred","SE","SeafType")

f1<-ggplot(data=probDOY,aes(x=DOYSeq,y=exp(z1)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(z1-1.96*z1SE),ymax=exp(z1+1.96*z1SE)),
		fill="#66CCFF",alpha=0.15)+
	scale_y_continuous("Standardized Abundance",limits=c(0,0.3))+
	scale_x_continuous("Day of the Year",limits=c(1,366),breaks=c(1,100,200,300))

f2<-ggplot(data=probGRT,aes(x=GRTSeq,y=exp(Pred)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(Pred-1.96*SE),ymax=exp(Pred+1.96*SE)),
		fill="#66CCFF",alpha=0.15)+
	facet_grid(Gear ~ .,scales="free_y")+
	scale_y_continuous("")+
	scale_x_continuous("ln-GRT",limits=c(-1,4))

f3<-ggplot(data=probDepth,aes(x=DepthSeq,y=exp(Pred)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(Pred-1.96*SE),ymax=exp(Pred+1.96*SE)),
		fill="#66CCFF",alpha=0.15)+
	facet_grid(Gear ~ .,scales="free_y")+
	scale_y_continuous("Standardized Abundance")+
	scale_x_continuous("ln-Depth (m)",limits=c(0,6))

f4<-ggplot(data=probNight,aes(x=NightSeq,y=exp(Pred)))+
	geom_line(colour="blue",lwd=1)+
	geom_ribbon(aes(ymin=exp(Pred-1.96*SE),ymax=exp(Pred+1.96*SE)),
		fill="#66CCFF",alpha=0.15)+
	facet_grid(Gear ~ .,scales="free_y")+
	scale_y_continuous("")+
	scale_x_continuous("Night soak (%)",limits=c(0,1))

pdf(file="/Users/jaimeoterovillar/Desktop/NBprobCont.pdf",width=10,height=10)

multiplot(f1,f3,f2,f4,cols=2)

dev.off()

#12.2# Categorical variables

#12.2.1# Seafloor

pdf(file="/Users/jaimeoterovillar/Desktop/NBprobFloor.pdf")

ggplot(data=probSeafloor,aes(x=SeafType,y=exp(Pred)))+
	geom_bar(stat="identity",colour="blue",fill="blue",fil="black")+
	scale_y_continuous("Standardized Abundance",limits=c(0,0.3))+
	scale_x_discrete("Seafloor type")+
	geom_errorbar(aes(ymax=exp(Pred+1.96*SE),ymin=exp(Pred-1.96*SE)),width=0,lwd=1)

dev.off()

#12.2.1# Year trend for each oceanographic zone

zYears.RB<-predict(negbin1,type="link",se=T,newdata=newTrend.RB)
zYears.AA<-predict(negbin1,type="link",se=T,newdata=newTrend.AA)
zYears.CA<-predict(negbin1,type="link",se=T,newdata=newTrend.CA)

abundInd.RB.fit<-zYears.RB$fit
abundInd.RB.SE<-zYears.RB$se.fit
abundInd.RB<-data.frame(cbind(abundInd.RB.fit,abundInd.RB.SE))
colnames(abundInd.RB)<-c("Fit","SE")

abundInd.AA.fit<-zYears.AA$fit
abundInd.AA.SE<-zYears.AA$se.fit
abundInd.AA<-data.frame(cbind(abundInd.AA.fit,abundInd.AA.SE))
colnames(abundInd.AA)<-c("Fit","SE")

abundInd.CA.fit<-zYears.CA$fit
abundInd.CA.SE<-zYears.CA$se.fit
abundInd.CA<-data.frame(cbind(abundInd.CA.fit,abundInd.CA.SE))
colnames(abundInd.CA)<-c("Fit","SE")

abundInd.GA<-rbind(abundInd.RB,abundInd.AA,abundInd.CA)

abundInd.GA$Year<-rep(1999:2013,times=3)
abundInd.GA$Zone<-rep(c("Rías Baixas","Arco Ártabro","Cantábrico"),each=15)
abundInd.GA$Zone<-factor(abundInd.GA$Zone,levels=c("Rías Baixas","Arco Ártabro","Cantábrico"))

pdf(file="/Users/jaimeoterovillar/Desktop/NBprobTrend.pdf",width=12,height=4)

ggplot(data=abundInd.GA,aes(x=Year,y=exp(Fit)))+
	geom_segment(aes(x=Year,y=exp(Fit),
				xend=Year,yend=exp(Fit+1.96*SE)),
				col="orange",lwd=0.5)+
	geom_segment(aes(x=Year,y=exp(Fit-1.96*SE),
				xend=Year,yend=exp(Fit)),
				col="orange",lwd=0.5)+
	geom_line(lwd=0.5,linetype="dotted",col="blue")+
	geom_point(size=5,col="gray90")+
	geom_point(size=3,col="orange")+
	facet_wrap(~Zone,scales="free_y")+
	scale_y_continuous("Standardized Abundance")+
	scale_x_continuous("Year",breaks=c(1999,2001,2003,2005,2007,2009,2011,2013))

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

library(glmmADMB)

negbinADMB<-glmmadmb(Ntot~fgear+
						  lGRT+
						  fyear+poly(Julian,2)+
						  lDepth+
						  caladoNight+fZoneO+Seafloor+
						  offset(log(offs1))+
						  (1|fvessel),
			 	 	      family="nbinom",data=pollack1)

summary(negbinADMB)

negbinADMBzi<-glmmadmb(Ntot~fgear+
						    lGRT+
						    fyear+poly(Julian,2)+
						    lDepth+
						    caladoNight+fZoneO+Seafloor+
						    offset(log(offs1))+
						    (1|fvessel),zeroInflation=T,
			 	 	        family="nbinom",data=pollack1)

summary(negbinADMBzi)

AIC(negbinADMB,negbinADMBzi)

