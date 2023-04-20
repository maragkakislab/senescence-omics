
library(readxl)
library(tidyverse)

Etop <- read_excel("2022-2-8 Etop.xlsx")
Etop_ISRIB<- read_excel("2022-2-8 Etop ISRIB.xlsx")
Cyc <- read_excel("2022-2-8 Cycling.xlsx")
Cyc_ISRIB <- read_excel("2022-2-8 Cycling ISRIB.xlsx")
Cyc_TG<- read_excel("2022-2-8 Cycling TG.xlsx")
Cyc_TG_ISRIB<- read_excel("2022-2-8 Cycling TG ISRIB.xlsx")


Etop$Time[which.max(Etop$Absorbance)]
Etop_ISRIB$Time[which.max(Etop_ISRIB$Absorbance)]
Cyc$Time[which.max(Cyc$Absorbance)]
Cyc_ISRIB$Time[which.max(Cyc_ISRIB$Absorbance)]
Cyc_TG$Time[which.max(Cyc_TG_ISRIB$Absorbance)]
Cyc_TG_ISRIB$Time[which.max(Cyc_TG_ISRIB$Absorbance)]



ggplot()+
  theme_classic()+
  theme(panel.border = element_rect(color="black",size=2,fill=NA))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(size=16,color="black")) +
  theme(axis.ticks.y = element_line(size=1.5))+
  theme(axis.title = element_text(size=16,color="black"))+
  theme(plot.title = element_text(size=14,color="black", hjust = 0.5))+
  geom_line(data=Etop,aes(x=Time+10,y=Absorbance),color="#B22222",size=1)+
  geom_line(data=Etop_ISRIB,aes(x=Time,y=Absorbance),color="steelblue",size=1)+
  geom_line(data=Cyc_TG,aes(x=Time,y=Absorbance),color="black",size=1)+
  geom_line(data=Cyc_TG_ISRIB,aes(x=Time,y=Absorbance),color="#004D40",size=1)+
  ylab("Normalized UV Absorbance")+
  ylab("UV Absorbance")+coord_cartesian(xlim=c(220,600),ylim=c(0,0.30))+
  xlab("")+geom_rect(aes(xmin=420,xmax=590,ymin=0,ymax=0.10),alpha=0.2,fill="blue", color="black")+
  geom_rect(aes(xmin=360,xmax=400,ymin=0,ymax=0.25),alpha=0.2,fill="grey", color="black")


Etop_Mono<-Etop%>%filter(between(Time,350,390))
Etop_Poly<-Etop%>%filter(between(Time,410,580))
sum(Etop_Poly$Absorbance)
sum(Etop_Mono$Absorbance)
sum(Etop_Poly$Absorbance)/sum(Etop_Mono$Absorbance)

2.463328/2.459551

EtopISRIB_Poly<-Etop_ISRIB%>%filter(between(Time,420,590))
EtopISRIB_Mono<-Etop_ISRIB%>%filter(between(Time,360,400))
sum(EtopISRIB_Poly$Absorbance)
sum(EtopISRIB_Mono$Absorbance)
sum(EtopISRIB_Poly$Absorbance)/sum(EtopISRIB_Mono$Absorbance)

CycTg_Poly<-Cyc_TG%>%filter(between(Time,420,590))
CycTg_Mono<-Cyc_TG%>%filter(between(Time,360,400))
sum(CycTg_Poly$Absorbance)
sum(CycTg_Mono$Absorbance)
sum(CycTg_Poly$Absorbance)/sum(CycTg_Mono$Absorbance)

CycTgISRIB_Poly<-Cyc_TG_ISRIB%>%filter(between(Time,420,590))
CycTgISRIB_Mono<-Cyc_TG_ISRIB%>%filter(between(Time,360,400))
sum(CycTgISRIB_Poly$Absorbance)
sum(CycTgISRIB_Mono$Absorbance)
sum(CycTgISRIB_Poly$Absorbance)/sum(CycTgISRIB_Mono$Absorbance)


2.330602/1.640682


sum(A$Absorbance)/sum(E$Absorbance)
sum(B$Absorbance)/sum(E$Absorbance)
sum(C$Absorbance)/sum(E$Absorbance)
sum(D$Absorbance)/sum(E$Absorbance)
sum(F$Absorbance)/sum(E$Absorbance)




A<-A%>%mutate(Normal=Absorbance*1.146822^-1)
B<-B%>%mutate(Normal=Absorbance*1.225099^-1)
C<-C%>%mutate(Normal=Absorbance*1.155575^-1)
D<-D%>%mutate(Normal=Absorbance*1.100739^-1)
E<-E%>%mutate(Normal=Absorbance*1)
F<-F%>%mutate(Normal=Absorbance*1.081467^-1)




ggplot()+
  theme_classic()+
  theme(panel.border = element_rect(color="black",size=2,fill=NA))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(size=26,color="black")) +
  theme(axis.ticks.y = element_line(size=1.5))+
  theme(axis.title = element_text(size=26,color="black"))+
  theme(plot.title = element_text(size=16,color="black", hjust = 0.5))+
  geom_line(data=A,aes(x=Time,y=Normal),color="#B22222",size=1)+
  geom_line(data=B,aes(x=Time,y=Normal),color="steelblue",size=1)+
  geom_line(data=C,aes(x=Time,y=Normal),color="black",size=1)+
  geom_line(data=D,aes(x=Time,y=Normal),color="#004D40",size=1)+
  ylab("UV Absorbance")+coord_cartesian(xlim=c(200,590),ylim=c(0,0.20))+
  xlab("")+geom_rect(aes(xmin=390,xmax=590,ymin=0,ymax=0.10),alpha=0.2,fill="blue", color="black")+
  geom_rect(aes(xmin=335,xmax=370,ymin=0,ymax=0.25),alpha=0.2,fill="grey", color="black")















################### Old Way
###############Fraction 5 is Fraction 1 in RTPCR
7+5

sum(E$Normal)

A<-Etop%>%filter(between(Fraction,8,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
B<-Etop_ISRIB%>%filter(between(Fraction,8,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
C<-Cyc%>%filter(between(Fraction,8,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
D<-Cyc_ISRIB%>%filter(between(Fraction,8,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
E<-Cyc_TG%>%filter(between(Fraction,8,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
F<-Cyc_TG_ISRIB%>%filter(between(Fraction,8,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)


ggplot()+geom_line(data=A,aes(x=Time,y=Absorbance))+
  geom_line(data=B,aes(x=Time,y=Absorbance),color="red")+scale_x_continuous(breaks=seq(200,600,25))


ggplot()+geom_line(data=C,aes(x=Time,y=Absorbance))+
  geom_line(data=D,aes(x=Time,y=Absorbance),color="red")+scale_x_continuous(breaks=seq(100,600,25))





ggplot()+
  theme_classic()+
  theme(panel.border = element_rect(color="black",size=2,fill=NA))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(face="bold",size=26,color="black")) +
  theme(axis.ticks.y = element_line(size=1.5))+
  theme(axis.title = element_text(face="bold",size=26,color="black"))+
  theme(plot.title = element_text(face="bold",size=16,color="black", hjust = 0.5))+
  geom_line(data=Etop,aes(x=Time,y=Absorbance),color="#B22222",size=2)+
  geom_line(data=Cyc,aes(x=Time,y=Absorbance),color="steelblue",size=2)+
  geom_line(data=Etop_ISRIB,aes(x=Time,y=Absorbance),color="black",size=2)+
  geom_line(data=Cyc_ISRIB,aes(x=Time,y=Absorbance),color="#004D40",size=2)+
  scale_x_continuous(breaks=seq(200,600,10))+ylab("Normalized UV Absorbance")+
  scale_y_continuous(breaks=seq(0,0.2,0.05))+ coord_cartesian(ylim=c(0.00,0.3),xlim = c(200,600))+
  xlab("")+ geom_rect(aes(xmin=410,xmax=550,ymin=0,ymax=0.18),alpha=0.2,fill="blue", color="black")+
  geom_rect(aes(xmin=350,xmax=390,ymin=0,ymax=0.18),alpha=0.2,fill="blue", color="black")


ggplot()+
  theme_classic()+
  theme(panel.border = element_rect(color="black",size=2,fill=NA))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(size=16,color="black"))+
  theme(axis.ticks.y = element_line(size=1))+
  theme(axis.title = element_text(size=16,color="black"))+
  theme(plot.title = element_text(size=14,color="black", hjust = 0.5))+
  geom_line(data=A,aes(x=Time,y=Normal-0.03),color="#B22222",size=2)+
  geom_line(data=B,aes(x=Time,y=Normal-0.03),color="steelblue",size=2)+
  geom_line(data=C,aes(x=Time,y=Normal-0.03),color="black",size=2)+
  geom_line(data=D,aes(x=Time,y=Normal-0.03),color="#004D40",size=2)+
  ylab("UV Absorbance")+coord_cartesian(xlim=c(200,590),ylim=c(0,0.20))+
  xlab("")+geom_rect(aes(xmin=390,xmax=590,ymin=0,ymax=0.10),alpha=0.2,fill="blue", color="black")+
  geom_rect(aes(xmin=335,xmax=370,ymin=0,ymax=0.25),alpha=0.2,fill="grey", color="black")


ggplot()+
  theme_classic()+
  theme(panel.border = element_rect(color="black",size=2,fill=NA))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(face="bold",size=26,color="black")) +
  theme(axis.ticks.y = element_line(size=1.5))+
  theme(axis.title = element_text(face="bold",size=26,color="black"))+
  theme(plot.title = element_text(face="bold",size=16,color="black", hjust = 0.5))+
  geom_line(data=E,aes(x=Time,y=Normal-0.03),color="#FFC107",size=2)+
  geom_line(data=F,aes(x=Time,y=Normal-0.03),color="steelblue",size=2)+
  geom_line(data=C,aes(x=Time,y=Normal-0.03),color="black",size=2)+
  geom_line(data=D,aes(x=Time,y=Normal-0.03),color="#004D40",size=2)+
  scale_x_continuous(breaks=seq(200,600,10))+ylab("Normalized UV Absorbance")+
  scale_y_continuous(breaks=seq(0,0.15,0.05))+ coord_cartesian(ylim=c(0.00,0.15))+
  xlab("")+geom_rect(aes(xmin=410,xmax=550,ymin=0,ymax=0.18),alpha=0.2,fill="blue", color="black")+
  geom_rect(aes(xmin=350,xmax=400,ymin=0,ymax=0.18),alpha=0.2,fill="blue", color="black")



sum(C$Normal[between(C$Time,410,550)])/sum(C$Normal[between(C$Time,350,400)])
sum(D$Normal[between(D$Time,410,550)])/sum(D$Normal[between(D$Time,350,400)])

sum(A$Normal[between(A$Time,410,550)])/sum(A$Normal[between(A$Time,350,400)])
sum(B$Normal[between(B$Time,410,550)])/sum(B$Normal[between(B$Time,350,400)])


sum(E$Normal[between(E$Time,410,550)])/sum(E$Normal[between(E$Time,350,400)])
sum(F$Normal[between(F$Time,410,550)])/sum(F$Normal[between(F$Time,350,400)])


