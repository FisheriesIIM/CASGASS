#################################################################
# Data base UTPB: CASGASS project (ICES Science Fund 2014-2015) #
#################################################################
# Date start: September 2014
# Last update: 20-01-2015
#updated data: 04-09-2014
#General objective: preliminary exploration, remove errors and prepare data for analysis. 
#
# Project: CASGASS
#
# Authors: Jaime Otero (jotero@iim.csic.es)
# Further modifications by: Alex Alonso (alex@iim.csic.es)
# Institute of Marine Research (Vigo)

#working directory
setwd("C:\\Users\\alex\\Dropbox\\casgass")

library(ggplot2)

#1# Cargar datos guardados en csv (desde el archivo ACCESS enviado por UTPB y conveniente modificado por Alex)
# separado por ";" para mantener formato de fechas
#datos originales en la correspondiente carpeta compartida de dropbox "CASGASS"

utpb<-read.csv(file="casgass_utpb_99-13.csv",header=T,dec=".",sep=";")

head(utpb)
dim(utpb)
summary(utpb)
str(utpb)


#2# Comprobar errores en variables clave

sort(unique(utpb$Zone)) # Ok
sort(unique(utpb$Harbour)) # Ok
sort(unique(utpb$Gear)) # RASTRO-VIEIRA,-VOLANDEIRA,-ZAMBURINA,-OSTRA??
sort(unique(utpb$Seafloor)) # algas,-roca-y-arena??
sort(unique(utpb$Species)) # Ok

naSoak<-which(is.na(utpb$Soak))
utpb<-utpb[-naSoak,] # Eliminar NAs en el tiempo de calado calculado por la UTPB
# Esos lances tienen errores o bien en la fecha de largada o de virada que
# no podemos solucionar.... Hemos solucionado muchos otros errores usando
# su "soak" y comparandolo con uno calculado por nosotros usando las fechas y
# horas de largada y virada (ver calculo mas abajo).
# una vez que vemos que no hay muchas desviaciones entre ambos
#utilizamos el calculado por ellos 


#3# Adaptar fechas (ejemplos en:
# http://www.noamross.net/blog/2014/2/10/using-times-and-dates-in-r---presentation-code.html) y calcular tiempo de calado del arte

utpb$FLARG<-as.POSIXct(utpb$FLARG,format="%Y/%m/%d %H:%M:%S")
utpb$FVIR<-as.POSIXct(utpb$FVIR,format="%Y/%m/%d %H:%M:%S")

utpb$Soaktime<-utpb$FVIR-utpb$FLARG
utpb$Soaktime<-as.numeric(difftime(utpb$FVIR,utpb$FLARG,units="mins"))

summary(utpb$Soak) # tiempo de calado UTPB
summary(utpb$Soaktime) # tiempo de calado IIM

plot(Soak~Soaktime,utpb,xlab="Soaktime (IIM)",ylab="Soak (UTPB)",xlim=c(0,60000),ylim=c(0,60000)) # Una vez que fuimos corrigiendo 
# manualmente algunos errores en las fechas de lances particulares (Alex los tiene anotados)
#se ve que el tiempo de calado proporcionado por la UTPB
# y el que calculamos nosotros se corresponden bien, pero usaremos el suyo

#4# Usamos como fecha de referencia la de la virada porque es el dia de muestreo

utpb$FVIR<-as.POSIXlt(utpb$FVIR,format="%Y/%m/%d %H:%M:%S")
utpb$Year<-utpb$FVIR$year+1900
summary(utpb$Year)
utpb$Month<-utpb$FVIR$mon+1
summary(utpb$Month)
utpb$Julian<-utpb$FVIR$yday+1
summary(utpb$Julian)

#5# Reducir los niveles de los tipos de fondo

utpb$Seafloor2<-ifelse(utpb$Seafloor=="piedra" |
                         utpb$Seafloor=="roca-y-algas" |
                         utpb$Seafloor=="piedra-dura","hard",
                       ifelse(utpb$Seafloor=="cascos-de-barcos-y/o-bateas" |
                                utpb$Seafloor=="piedra-y-arena" |
                                utpb$Seafloor=="algas,-roca-y-arena" |
                                utpb$Seafloor=="piedra-y-fango" |
                                utpb$Seafloor=="fango-pedra-cascallo" |
                                utpb$Seafloor=="piedra-y-coral" |
                                utpb$Seafloor=="pedra-cascallo" |
                                utpb$Seafloor=="arena-coral-piedra","mixed","soft"))
utpb$Seafloor2<-as.factor(utpb$Seafloor2)

#6# Comprobar los NA de las profundidades

plot(DepthMax~DepthMin,data=utpb)
abline(1,1) # Parece que hay datos cambiados porque DepthMin nunca puede
# ser mayor que DepthMax

# Primero recolocar columnas
utpb$DepthMinNew<-ifelse(utpb$DepthMin>utpb$DepthMax,utpb$DepthMax,utpb$DepthMin)
utpb$DepthMaxNew<-ifelse(utpb$DepthMax<utpb$DepthMin,utpb$DepthMin,utpb$DepthMax)

plot(DepthMaxNew~DepthMinNew,data=utpb)
abline(1,1)

# Segundo combinar columnas 
utpb$DepthMinNew2<-ifelse(is.na(utpb$DepthMinNew),utpb$DepthMin,utpb$DepthMinNew)
utpb$DepthMaxNew2<-ifelse(is.na(utpb$DepthMaxNew),utpb$DepthMax,utpb$DepthMaxNew)

plot(DepthMaxNew2~DepthMinNew2,data=utpb)
abline(1,1)

# Tercero obtener media de profundidad y cubrir NA con el dato que haya
utpb$depth<-(utpb$DepthMaxNew2+utpb$DepthMinNew2)/2
utpb$depthNew1<-ifelse(is.na(utpb$depth),utpb$DepthMinNew2,utpb$depth)
utpb$depthNew2<-ifelse(is.na(utpb$depthNew1),utpb$DepthMaxNew2,utpb$depthNew1)

par(mfcol=c(1,2))
plot(DepthMaxNew2~depthNew2,data=utpb)
abline(1,1)
plot(DepthMinNew2~depthNew2,data=utpb)
abline(1,1)


#7# Añadir otras variables de interes

utpb$Ntot<-utpb$NUMc+utpb$NUMd # Numero total de individuos (retained+discarded)

utpb$Wtot<-(utpb$PESOc+utpb$PESOd)/1000 # Peso total (kg) de los individuos (retained+discarded)
utpb$Wtot<-ifelse(utpb$Wtot==0,NA,utpb$Wtot) # kg con valor 0 son NA


#8# Generar columnas nuevas para combinar LATin y LONin con LATfin y LONfin
# ya que hay NAs en coordenadas. Cambiar unidades de las coordenadas

utpb$LATdef<-ifelse(is.na(utpb$LATin),utpb$LATfin,utpb$LATin)
utpb$LONdef<-ifelse(is.na(utpb$LONin),utpb$LONfin,utpb$LONin)

plot(LATdef~LONdef,data=utpb)
utpb$Lat<-trunc(utpb$LATdef/100000)+((utpb$LATdef-(100000*trunc(utpb$LATdef/100000)))/1000)/60
utpb$Lon<--(trunc(utpb$LONdef/100000)+((utpb$LONdef-(100000*trunc(utpb$LONdef/100000)))/1000)/60)
plot(Lat~Lon,data=utpb,pch=16,cex=0.4,col=Zone)
# Todos los datos que a priori parecian raros fueron
# corregidos usando la profundidad de calado, el tipo de arte, puerto etc 
# (anotados por Alex). Los que ahora parece que quedan alejados creemos que 
# estan bien. No obstante hay errores en algunas de las clasificaciones por zonas
# con lo que las arreglamos manualmente en funcion de la Lat y Lon


#9# Nueva asignacion de Zonas Administrativas y Zonas Oceanograficas

#9.1# Coordenadas decimales de los puntos de separación

BaresLat<-43.791336
BaresLon<--7.688769

PriorinhoLat<-43.463275
PriorinhoLon<--8.348032

LangosteiraLat<-43.361084
LangosteiraLon<--8.486756

TourinhanLat<-43.059324
TourinhanLon<--9.293350

InsuaLat<-42.770942
InsuaLon<--9.127085

SieiraLat<-42.653324
SieiraLon<--9.042177

FaxildaLat<-42.414997
FaxildaLon<--8.881116

SoavelaLat<-42.277784
SoavelaLon<--8.864851

#9.2# Zonas Administrativas

utpb$ZoneA<-ifelse(utpb$Lon >= BaresLon,9,8) # Zona 9 (A Mariña)
utpb$ZoneA<-ifelse(utpb$Lat < BaresLon & utpb$Lat >= 43.51,8,utpb$ZoneA) # Zona 8 (Cedeira)
utpb$ZoneA<-ifelse(utpb$Lon >= LangosteiraLon & utpb$Lat <= 43.51,7,utpb$ZoneA) # Zona 7 (Coruña-Ferrol)
utpb$ZoneA<-ifelse(utpb$Lon < LangosteiraLon & utpb$Lat >= TourinhanLat,6,utpb$ZoneA) # Zona 6 (Costa da Morte)
utpb$ZoneA<-ifelse(utpb$Lat < TourinhanLat & utpb$Lat >= InsuaLat & utpb$Lon <= -9.08,5,utpb$ZoneA) # Zona 5 (Fisterra)
utpb$ZoneA<-ifelse(utpb$Lat < InsuaLat & utpb$Lat >= SieiraLat,4,utpb$ZoneA) # Zona 4 (Muros)
utpb$ZoneA<-ifelse(utpb$Lon > -9.08 & utpb$Lon <= -8.8 & utpb$Lat > InsuaLat & utpb$Lat < 42.85,4,utpb$ZoneA) # Zona 4 (Interior Ria Muros)
utpb$ZoneA<-ifelse(utpb$Lat < SieiraLat & utpb$Lat >= FaxildaLat,3,utpb$ZoneA) # Zona 3 (Arousa)
utpb$ZoneA<-ifelse(utpb$Lon > -8.9 & utpb$Lon <= -8.7 & utpb$Lat > SieiraLat & utpb$Lat < 42.7,3,utpb$ZoneA) # Zona 3 (Interior Ria Arousa)
utpb$ZoneA<-ifelse(utpb$Lat < FaxildaLat & utpb$Lat >= SoavelaLat,2,utpb$ZoneA) # Zona 2 (Pontevedra)
utpb$ZoneA<-ifelse(utpb$Lon > -8.75 & utpb$Lon <= -8.65 & utpb$Lat > FaxildaLat & utpb$Lat < 42.45,2,utpb$ZoneA) # Zona 2 (Interior Ria Pontevedra 1)
utpb$ZoneA<-ifelse(utpb$Lon > -8.8229 & utpb$Lon <= -8.822 & utpb$Lat > 42.277 & utpb$Lat < 42.278,2,utpb$ZoneA) # Zona 2 (Interior Ria Pontevedra 2)
utpb$ZoneA<-ifelse(utpb$Lat < SoavelaLat,1,utpb$ZoneA) # Zona 1 (Vigo)
utpb$ZoneA<-ifelse(utpb$Lon > -8.68 & utpb$Lon <= -8.5 & utpb$Lat > SoavelaLat & utpb$Lat < 42.36,1,utpb$ZoneA) # Zona 1 (Interior Ria Vigo 1)
utpb$ZoneA<-ifelse(utpb$Lon > -8.74 & utpb$Lon <= -8.67 & utpb$Lat > SoavelaLat & utpb$Lat < 42.3,1,utpb$ZoneA) # Zona 1 (Interior Ria Vigo 2)

#9.3# Zonas Oceanograficas

utpb$ZoneO<-ifelse(utpb$ZoneA >= 8,3,ifelse(utpb$ZoneA==7 | utpb$ZoneA==6,2,1))

#9.4# Checking Plot

p1<-ggplot(data=utpb,aes(x=Lon,y=Lat))
p1+geom_point(aes(colour=factor(ZoneA)))

p2<-ggplot(data=utpb,aes(x=Lon,y=Lat))
p2+geom_point(aes(colour=factor(ZoneO)))


#10# Reordenar columnas y renombrar

head(utpb)
names(utpb)

utpb<-utpb[,c("Idjornada","Idlance","Year","Month","Julian","Observer",
              "Harbour","Idflota","TRB","Trips","Gear","Pieces","FLARG","FVIR",
              "Soak","Lat","Lon","ZoneA","ZoneO","Seafloor2","depthNew2","Species",
              "Ntot","Wtot")] # Reordenar

colnames(utpb)<-c("Idjornada","Idlance","Year","Month","Julian","Observer",
                  "Harbour","Idflota","GRT","Crew","Gear","Pieces","Deployment",
                  "Retrieval","Soak","Lat","Lon","ZoneA","ZoneO","Seafloor","Depth",
                  "Species","Ntot","Wtot") # Renombrar

head(utpb)
dim(utpb)
summary(utpb)
str(utpb)


#11# Base de datos independiente para cada especie
#abrimos la base de datos completa para separar por especies

#estadistica descriptiva
library(pastecs)
stat.desc(utpb, norm=F)

congrio<-utpb[utpb$Species=="Conger-conger",]
lubina<-utpb[utpb$Species=="Dicentrarchus-labrax",]
sargo<-utpb[utpb$Species=="Diplodus-sargus",]
pinto<-utpb[utpb$Species=="Labrus-bergylta",]
centolla<-utpb[utpb$Species=="Maja-brachydactyla",]
salmonete<-utpb[utpb$Species=="Mullus-surmuletus",]
pulpo<-utpb[utpb$Species=="Octopus-vulgaris",]
platija<-utpb[utpb$Species=="Platichthys-flesus",]
abadejo<-utpb[utpb$Species=="Pollachius-pollachius",]
raya<-utpb[utpb$Species=="Raja-undulata",]
pintaroja<-utpb[utpb$Species=="Scyliorhinus-canicula",]
lenguado<-utpb[utpb$Species=="Solea-vulgaris",]
faneca<-utpb[utpb$Species=="Trisopterus-luscus",]


#12# Ejemplos de distribucion de capturas por artes y especies

congrios<-as.numeric(tapply(congrio$Ntot,congrio$Gear,sum,na.rm=T))
lubinas<-as.numeric(tapply(lubina$Ntot,lubina$Gear,sum,na.rm=T))
sargos<-as.numeric(tapply(sargo$Ntot,sargo$Gear,sum,na.rm=T))
pintos<-as.numeric(tapply(pinto$Ntot,pinto$Gear,sum,na.rm=T))
centollas<-as.numeric(tapply(centolla$Ntot,centolla$Gear,sum,na.rm=T))
salmonetes<-as.numeric(tapply(salmonete$Ntot,salmonete$Gear,sum,na.rm=T))
pulpos<-as.numeric(tapply(pulpo$Ntot,pulpo$Gear,sum,na.rm=T))
platijas<-as.numeric(tapply(platija$Ntot,platija$Gear,sum,na.rm=T))
abadejos<-as.numeric(tapply(abadejo$Ntot,abadejo$Gear,sum,na.rm=T))
rayas<-as.numeric(tapply(raya$Ntot,raya$Gear,sum,na.rm=T))
pintarojas<-as.numeric(tapply(pintaroja$Ntot,pintaroja$Gear,sum,na.rm=T))
lenguados<-as.numeric(tapply(lenguado$Ntot,lenguado$Gear,sum,na.rm=T))
fanecas<-as.numeric(tapply(faneca$Ntot,faneca$Gear,sum,na.rm=T))

par(mar=c(15,5,2,2))
barplot((congrios*100)/sum(congrio$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Congrio")
#quartz()
par(mar=c(15,5,2,2))
barplot((lubinas*100)/sum(lubina$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Lubina")
#quartz()
par(mar=c(15,5,2,2))
barplot((sargos*100)/sum(sargo$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Sargo")
#quartz()
par(mar=c(15,5,2,2))
barplot((pintos*100)/sum(pinto$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Pinto")
#quartz()
par(mar=c(15,5,2,2))
barplot((centollas*100)/sum(centolla$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Centolla")
#quartz()
par(mar=c(15,5,2,2))
barplot((salmonetes*100)/sum(salmonete$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Salmonete")
#quartz()
par(mar=c(15,5,2,2))
barplot((pulpos*100)/sum(pulpo$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Pulpo")
#quartz()
par(mar=c(15,5,2,2))
barplot((platijas*100)/sum(platija$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Platija")
#quartz()
par(mar=c(15,5,2,2))
barplot((abadejos*100)/sum(abadejo$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Abadejo")
#quartz()
par(mar=c(15,5,2,2))
barplot((rayas*100)/sum(raya$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Raya")
#quartz()
par(mar=c(15,5,2,2))
barplot((pintarojas*100)/sum(pintaroja$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Pintaroja")
#quartz()
par(mar=c(15,5,2,2))
barplot((lenguados*100)/sum(lenguado$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Lenguado")
#quartz()
par(mar=c(15,5,2,2))
barplot((fanecas*100)/sum(faneca$Ntot,na.rm=T),ylim=c(0,100),names.arg=sort(unique(utpb$Gear)),ylab="Frecuencia",xlab="",las=2,main="Faneca")

#graph para el ICES progress report
#Alex
ggplot(utpb, aes(Gear))+ geom_bar()+facet_wrap(~ Species, ncol=5)+xlab("Fishing gear type")+ylab("Numaber of hauls sampled")+
  theme(axis.text.x = element_blank(),axis.title.x = element_text(face="bold", size=15))

#13# Exportar datos para cada especies

write.table(congrio,file="congrio.txt",sep=",",dec=".",row.names=F)

write.table(lubina,file="lubina.txt",sep=",",dec=".",row.names=F)

write.table(sargo,file="sargo.txt",sep=",",dec=".",row.names=F)

write.table(pinto,file="pinto.txt",sep=",",dec=".",row.names=F)

write.table(centolla,file="centolla.txt",sep=",",dec=".",row.names=F)

write.table(salmonete,file="salmonete.txt",sep=",",dec=".",row.names=F)

write.table(pulpo,file="pulpo.txt",sep=",",dec=".",row.names=F)

write.table(platija,file="platija.txt",sep=",",dec=".",row.names=F)

write.table(abadejo,file="abadejo.txt",sep=",",dec=".",row.names=F)

write.table(raya,file="raya.txt",sep=",",dec=".",row.names=F)

write.table(pintaroja,file="pintaroja.txt",sep=",",dec=".",row.names=F)

write.table(lenguado,file="lenguado.txt",sep=",",dec=".",row.names=F)

write.table(faneca,file="faneca.txt",sep=",",dec=".",row.names=F)
