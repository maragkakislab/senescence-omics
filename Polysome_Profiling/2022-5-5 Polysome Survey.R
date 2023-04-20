
library(readxl)
library(tidyverse)
getwd()

SIQ <- read_excel("2022-5-2 SIQ.xlsx")
Qui<- read_excel("2022-5-2 Qui.xlsx")
IR <- read_excel("2022-5-2 IR.xlsx")
Etop <- read_excel("2022-5-2 Etop.xlsx")
Tg<- read_excel("2022-5-2 Tg.xlsx")
Cyc<- read_excel("2022-5-2 Cyc.xlsx")


colnames(SIQ)<-c("Distance","Absorbance","Fraction")
colnames(Qui)<-c("Distance","Absorbance","Fraction")
colnames(IR)<-c("Distance","Absorbance","Fraction")
colnames(Etop)<-c("Distance","Absorbance","Fraction")
colnames(Tg)<-c("Distance","Absorbance","Fraction")
colnames(Cyc)<-c("Distance","Absorbance","Fraction")




###############Fraction XXX is Fraction 1 in RTPCR


A<-SIQ%>%filter(between(Fraction,8,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
B<-Qui%>%filter(between(Fraction,7,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
C<-IR%>%filter(between(Fraction,7,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
D<-Etop%>%filter(between(Fraction,7,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
E<-Tg%>%filter(between(Fraction,8,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
F<-Cyc%>%filter(between(Fraction,7,15))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)


ggplot()+geom_line(data=E,aes(x=Distance,y=Normal))+
  geom_line(data=F,aes(x=Distance,y=Normal),color="red")+scale_x_continuous(breaks=seq(200,600,50))



A<-SIQ%>%filter(between(Distance,225,600))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
B<-Qui%>%filter(between(Distance,232,600))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
C<-IR%>%filter(between(Distance,176,550))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
D<-Etop%>%filter(between(Distance,219,600))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
E<-Tg%>%filter(between(Distance,224,600))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)
F<-Cyc%>%filter(between(Distance,225,600))%>%mutate(Normal=Absorbance/sum(Absorbance)*100)



##########################

ggplot()+
  theme_classic()+
  theme(panel.border = element_rect(color="black",size=2,fill=NA))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(face="bold",size=16,color="black")) +
  theme(axis.ticks.y = element_line(size=1.5))+
  theme(axis.title = element_text(face="bold",size=16,color="black"))+
  theme(plot.title = element_text(face="bold",size=14,color="black", hjust = 0.5))+
  geom_line(data=D,aes(x=Distance,y=Normal),color="Firebrick",size=1)+
  geom_line(data=C,aes(x=Distance+40,y=Normal),color="Steelblue",size=1)+
  geom_line(data=F,aes(x=Distance,y=Normal),color="black",size=1)+
  scale_x_continuous(breaks=seq(200,600,10))+ylab("Normalized UV Absorbance")+
  xlab("")+
  geom_rect(aes(xmin=390,xmax=550,ymin=0,ymax=0.18),alpha=0.2,fill="blue", color="black")+
  geom_rect(aes(xmin=290,xmax=335,ymin=0,ymax=0.18),alpha=0.2,fill="grey", color="black")

############## Distances
#polysome: +165 - +325
#monosome: +65 - +110


F_Poly<-F%>%filter(between(Distance,390,600))
F_Mono<-F%>%filter(between(Distance,290,335))
sum(F_Poly$Normal)
sum(F_Mono$Normal)
sum(F_Poly$Absorbance)/sum(F_Mono$Absorbance)


E_Poly<-E%>%filter(between(Distance,400,610))
E_Mono<-E%>%filter(between(Distance,300,345))
sum(E_Poly$Absorbance)/sum(E_Mono$Absorbance)

D_Poly<-D%>%filter(between(Distance,400,610))
D_Mono<-D%>%filter(between(Distance,300,345))
sum(D_Poly$Normal)
sum(D_Mono$Normal)
sum(D_Poly$Absorbance)/sum(D_Mono$Absorbance)

C_Poly<-C%>%filter(between(Distance,350,560))
C_Mono<-C%>%filter(between(Distance,250,295))
sum(C_Poly$Normal)
sum(C_Mono$Normal)
sum(C_Poly$Absorbance)/sum(C_Mono$Absorbance)

B_Poly<-B%>%filter(between(Distance,400,610))
B_Mono<-B%>%filter(between(Distance,300,345))
sum(B_Poly$Normal)
sum(B_Mono$Normal)
sum(B_Poly$Absorbance)/sum(B_Mono$Absorbance)



sum(B_Poly$Normal)
sum(B_Mono$Normal)

ggplot()+
  theme_classic()+
  theme(panel.border = element_rect(color="black",size=1,fill=NA))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(face="bold",size=26,color="black")) +
  theme(axis.ticks.y = element_line(size=1.5))+
  theme(axis.title = element_text(face="bold",size=26,color="black"))+
  theme(plot.title = element_text(face="bold",size=16,color="black", hjust = 0.5))+
  geom_line(data=B,aes(x=Distance-10,y=Normal),color="purple",size=2)+
  geom_line(data=D,aes(x=Distance,y=Normal),color="red",size=2)+
  geom_line(data=C,aes(x=Distance+40,y=Normal),color="green4",size=2)+
  scale_x_continuous(breaks=seq(200,600,10))+ylab("Normalized UV Absorbance")+
  xlab("")

