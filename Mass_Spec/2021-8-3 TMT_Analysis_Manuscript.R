library(tidyverse)
library(readxl)
library(writexl)

list.files()
TMT_16<-read_xlsx("Appendix2-P759-VT23-TMT16PLEX-Data-Statistic-analysis.xlsx")

summary(TMT_16)

sdt<-function(x,y,z){
  x<-c(x,y,z)
  sd(x)}

tTest<-function(a,b,c,d,e,f){
  x<-c(a,b,c)
  y<-c(d,e,f)
  t.test(x,y)$p.value
}


TMT_High<-TMT_16%>%filter(`Protein FDR Confidence: Combined`!= "Low")%>%mutate(sd_Pro=mapply(sdt,P10,P11,P12),
                                                                               sd_Qui=mapply(sdt,Q4,Q5,Q6),
                                                                               sd_EM=mapply(sdt,EM1,EM2,EM3),
                                                                               p_Pro_Qui=mapply(tTest,Q4,Q5,Q6,P10,P11,P12),
                                                                               p_EM_Qui=mapply(tTest,EM1,EM2,EM3,Q4,Q5,Q6),
                                                                               p_EM_Pro=mapply(tTest,EM1,EM2,EM3,P10,P11,P12))

data.frame(colnames(TMT_High))
TMT_averages<-TMT_High%>%dplyr::select(2,9,26,27,31,32,33,60:65)
names(TMT_averages)<-c("Accession","Description","Pro","Qui","Etop","EtopM","EtopT","sd_Pro","sd_Qui","sd_EM",
                       "p_Pro_Qui","p_EM_Qui","p_EM_Pro")



Volcano<-TMT_averages%>%mutate(FC2_QuiPro=log2(Qui/Pro),
                               FC2_EtopMPro=log2(EtopM/Pro),
                               FC2_EtopMQui=log2(EtopM/Qui),
                               FC2_EtopMEtopT=log2(EtopM/EtopT))%>%
  filter(sd_EM/EtopM<=0.10, sd_Pro/Pro <=0.1, sd_Qui/Qui <=0.1)

######################################HGNC names
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

filters<-listFilters(ensembl)
attributes<-listAttributes(ensembl)

Accession<-Volcano$Accession

ID<-getBM(attributes = c('hgnc_symbol','uniprot_gn_id'),
          filters = 'uniprot_gn_id',
          values = Accession, 
          mart = ensembl)

Volcano<-Volcano%>%left_join(ID,by=c("Accession"="uniprot_gn_id"))



################################################################################################ Volcano Plots
######Saving 8.73 x 5.73 in image
sd(Volcano$FC2_EtopMPro)

sd(Volcano$FC2_EtopMQui)



low<-log2(0.75)
high<-log2(1.25)
ulow<-log2(0.9)
uhigh<-log2(1.1)

cycle<-c("CDKN1A","CDKN2A","CDKN1C","PCNA","LMNB1")
Cycle_Select<-Volcano%>%filter(hgnc_symbol%in%cycle)


Cyc_common<-Volcano%>%filter(between(FC2_EtopMPro,ulow,uhigh))
Cyc_modest<-Volcano%>%filter(!between(Volcano$FC2_EtopMPro,ulow,uhigh)&between(Volcano$FC2_EtopMPro,low,high))
Cyc_substantial<-Volcano%>%filter(!between(Volcano$FC2_EtopMPro,low,high))
Qui_common<-Volcano%>%filter(between(FC2_EtopMQui,ulow,uhigh))
Qui_modest<-Volcano%>%filter(!between(Volcano$FC2_EtopMQui,ulow,uhigh)&between(Volcano$FC2_EtopMQui,low,high))
qui_substantial<-Volcano%>%filter(!between(Volcano$FC2_EtopMQui,low,high))






######## Etop vs. Cyc
Fig3A<-ggplot()+
geom_point(data=Volcano,
             aes(x=FC2_EtopMPro,y=-log10(p_EM_Pro)),
             color=case_when(!between(Volcano$FC2_EtopMPro,ulow,uhigh)&between(Volcano$FC2_EtopMPro,low,high)~"tan3",
                             !between(Volcano$FC2_EtopMPro,low,high)~"black",
                             TRUE~"grey70"),
             alpha=0.6)+
  geom_point(data=Cycle_Select,
             aes(x=FC2_EtopMPro,y=-log10(p_EM_Pro)),
             color="black",fill="red",size=3,shape=22)+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  geom_hline(yintercept = -log10(0.01),color="dark red",size=2,linetype=2)+
  xlab("Log2FC (Senescent/Cycling)")+ ylab("-log10 p-value")+
  theme_classic()+
  coord_cartesian(xlim = c(-3,3),ylim=c(0,8))+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+scale_x_continuous(breaks=c(seq(-3,3,1)))+
  geom_label(data=Cycle_Select,aes(x=FC2_EtopMPro,y=-log10(p_EM_Pro),label=hgnc_symbol),nudge_y = 0.4,
             nudge_x=-0.1,fontface="bold",fill="white",size=6)
ggsave("Fig. 3A.tiff", Fig3A,dpi=600)
Fig3A
?ggsave()
options("device")

######## Etop vs. Qui
Fig3B<-ggplot()+
  geom_point(data=Volcano,
             aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui)),
             color=case_when(!between(Volcano$FC2_EtopMQui,ulow,uhigh)&between(Volcano$FC2_EtopMQui,low,high)~"tan3",
                             !between(Volcano$FC2_EtopMQui,low,high)~"black",
                             TRUE~"grey"),
             alpha=0.6)+
geom_point(data=Cycle_Select,
           aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui)),
           color="black",fill="red",size=3,shape=22)+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  geom_hline(yintercept = -log10(0.01),color="dark red",size=2,linetype=2)+
  xlab("Log2FC (Senescent/Quiescent)")+ ylab("-log10 p-value")+
  theme_classic()+
  coord_cartesian(xlim = c(-3,3),ylim=c(0,8))+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+scale_x_continuous(breaks=c(seq(-3,3,1)))+
  geom_label(data=Cycle_Select,aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui),label=hgnc_symbol),nudge_y = 0.4,
             nudge_x=-0.1,fontface="bold",fill="white",size=6)
Fig3B
ggsave("Fig. 3B.tiff", Fig3B,dpi=600)


################################################ Ribosomal Proteins
#################### 

RProteins<-Volcano%>%filter(str_detect(Volcano$Description,"ribosomal"))%>%
  filter(str_detect(hgnc_symbol,"^RPS|RPL"))%>%filter(!str_detect(Description,"mitochondrial"))%>%
  filter(!str_detect(hgnc_symbol,"M"))%>%filter(!str_detect(Description,"like"))


Volcano_Ribo<-Volcano%>%filter(hgnc_symbol%in%RProteins$hgnc_symbol)


Fig3D<-ggplot()+
  geom_point(data=Volcano,
             aes(x=FC2_EtopMPro,y=-log10(p_EM_Pro)),
             color=case_when(!between(Volcano$FC2_EtopMPro,ulow,uhigh)&between(Volcano$FC2_EtopMPro,low,high)~"tan3",
                             !between(Volcano$FC2_EtopMPro,low,high)~"black",
                             TRUE~"grey"),
             alpha=0.6)+
  geom_point(data=Volcano_Ribo,
             aes(x=FC2_EtopMPro,y=-log10(p_EM_Pro)),
             color="black",fill="red",size=3,shape=22)+
  coord_cartesian(xlim = c(-3,3),ylim=c(0,8))+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  geom_hline(yintercept = -log10(0.01),color="dark red",size=2,linetype=2)+
  xlab("Log2FC (Senescent/Cycling)")+ ylab("-log10 p-value")+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))


Fig3D
ggsave("Fig. 3D.tiff", Fig3D,dpi=600)


Fig3E<-ggplot()+
  geom_point(data=Volcano,
             aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui)),
             color=case_when(!between(Volcano$FC2_EtopMQui,ulow,uhigh)&between(Volcano$FC2_EtopMQui,low,high)~"tan3",
                             !between(Volcano$FC2_EtopMQui,low,high)~"black",
                             TRUE~"grey"),
             alpha=0.6)+
  geom_point(data=Volcano_Ribo,
             aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui)),
             color="black",fill="red",size=3,shape=22)+
  coord_cartesian(xlim = c(-3,3),ylim=c(0,8))+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  geom_hline(yintercept = -log10(0.01),color="dark red",size=2,linetype=2)+
  xlab("Log2FC (Senescent/Quiescent)")+ ylab("-log10 p-value")+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))
Fig3E
ggsave("Fig3E.tiff", Fig3E,dpi=600)


################################################ rRNA Processing Proteins
rRNA<-read_tsv("rRNA_Processing_Reactome.tsv")
rRNA_Prot<-Volcano%>%filter(Accession%in%rRNA$Identifier)%>%filter(!hgnc_symbol%in%RProteins$hgnc_symbol)

Fig4A<-ggplot()+
  geom_point(data=Volcano,
             aes(x=FC2_EtopMPro,y=-log10(p_EM_Pro)),
             color=case_when(!between(Volcano$FC2_EtopMPro,ulow,uhigh)&between(Volcano$FC2_EtopMPro,low,high)~"tan3",
                             !between(Volcano$FC2_EtopMPro,low,high)~"black",
                             TRUE~"grey"),
             alpha=0.6)+
  geom_point(data=rRNA_Prot,
             aes(x=FC2_EtopMPro,y=-log10(p_EM_Pro)),
             color="black",fill="red",size=3,shape=22)+
  coord_cartesian(xlim = c(-3,3),ylim=c(0,8))+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  geom_hline(yintercept = -log10(0.01),color="dark red",size=2,linetype=2)+
  xlab("Log2FC (Senescent/Cycling)")+ ylab("-log10 p-value")+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))


Fig4A
ggsave("Fig4A.tiff", Fig4A,dpi=600)


Fig4B<-ggplot()+
  geom_point(data=Volcano,
             aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui)),
             color=case_when(!between(Volcano$FC2_EtopMQui,ulow,uhigh)&between(Volcano$FC2_EtopMQui,low,high)~"tan3",
                             !between(Volcano$FC2_EtopMQui,low,high)~"black",
                             TRUE~"grey"),
             alpha=0.6)+
  geom_point(data=rRNA_Prot,
             aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui)),
             color="black",fill="red",size=3,shape=22)+
  coord_cartesian(xlim = c(-3,3),ylim=c(0,8))+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  geom_hline(yintercept = -log10(0.01),color="dark red",size=2,linetype=2)+
  xlab("Log2FC (Senescent/Quiescent)")+ ylab("-log10 p-value")+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))
Fig4B
ggsave("Fig4B.tiff", Fig4B,dpi=600)


################################################ Translation
select<-c("EIF2AK4","EEF2K","EEF1A1","EEF1A2","EIF5","EIF4G2","EIF2D")
ISR_Select<-Volcano%>%filter(hgnc_symbol%in%select)

######## Etop vs. Qui
Fig4C<-ggplot()+
  geom_point(data=Volcano,
             aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui)),
             color=case_when(!between(Volcano$FC2_EtopMQui,ulow,uhigh)&between(Volcano$FC2_EtopMQui,low,high)~"tan3",
                             !between(Volcano$FC2_EtopMQui,low,high)~"black",
                             TRUE~"grey"),
             alpha=0.6)+
  geom_point(data=ISR_Select,
             aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui)),
             color="black",fill="red",size=3,shape=22)+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  geom_hline(yintercept = -log10(0.01),color="dark red",size=2,linetype=2)+
  xlab("Log2FC (Senescent/Quiescent)")+ ylab("-log10 p-value")+
  theme_classic()+
  coord_cartesian(xlim = c(-3,3),ylim=c(0,8))+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+scale_x_continuous(breaks=c(seq(-3,3,1)))+
  geom_label(data=ISR_Select,aes(x=FC2_EtopMQui,y=-log10(p_EM_Qui),label=hgnc_symbol),
             nudge_y = c(0.1,0.1,0.3,0.1,0.1,0.1,0.1),
             nudge_x=-0.3,fontface="bold",fill="white",size=4)
Fig4C
ggsave("Fig4C.tiff", Fig4C,dpi=600)

