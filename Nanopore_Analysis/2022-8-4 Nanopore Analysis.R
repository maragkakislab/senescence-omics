library(tidyverse)
library(readxl)
library(biomaRt)
library(writexl)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(reshape)


RNA<-read_tsv("CYC_SEN_DESEQ2.tab")%>%filter(baseMean>1)%>%na.omit()


RNA%>%filter(log2FoldChange>=1|log2FoldChange<=-1)%>%filter(padj<=0.05)%>%nrow()
plot(RNA$log2FoldChange,log(RNA$baseMean,2))
plot(RNA$log2FoldChange*-log10(RNA$pvalue),RNA$baseMean)


RNA%>%filter(log2FoldChange>=1)%>%filter(padj<=0.05)%>%nrow()
RNA%>%filter(log2FoldChange<=-1)%>%filter(padj<=0.05)%>%nrow()


writexl::write_xlsx(RNA,"E_C_dRNA.xlsx")


############################################################ Volcano Plot
RNA$Norm<-RNA$log2FoldChange*log(RNA$baseMean,2)
head(RNA[order(RNA$Norm,decreasing = T),],10)


topten<-c("FTL","FTH1","MT2A","QSOX1","CD44","CDKN1A","SERPINE1","BEST1","CCND1","RND3",
          "TGFBI","PTMA","S100A10","H2AFZ","HMGN2","HMGB1","PCOLCE","NREP","PTGS1","EIF4EBP1")


topnormal<-c("MT2A","CCND2","CCND1","IGFBP2","KRTAP2-3","SERPINE1","CLDN1","HMGA2","KRTAP1-5","CDKN1A",
             "MYOCD","TK1","HMGB2","LXN","FOS","TMSB15A","PIMREG","PCOLCE","TOP2A","UBE2C")

topten_RNA<-RNA%>%filter(genename%in%topten)
topnormal_RNA<-RNA%>%filter(genename%in%topnormal)

ggplot()+
  geom_point(data=RNA,
             aes(x=log2FoldChange,y=-log10(padj)),alpha=0.6,
             color=case_when(RNA$log2FoldChange< -1 & RNA$pvalue<0.05~"#B22222",
                             RNA$log2FoldChange>1&RNA$pvalue<0.05~"#004D40",
                             T~"black"))+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  xlab("Log2FC RNA (Etop/Cyc)")+ ylab("-log10 p-value")+
  theme_classic()+
  coord_cartesian(xlim=c(-6,6),ylim=c(0,100))+
  theme(axis.text.x = element_text(size=24,color="black"),
        axis.text.y = element_text(size=24,color="black"),
        axis.title.x = element_text(size=20,color="black"),
        axis.title.y = element_text(size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  geom_text_repel(data=topnormal_RNA,aes(x=log2FoldChange,y=-log10(padj),label=genename),
                   nudge_x=0,size=3,min.segment.length = 0)


######################## SEN Validation
SEN<-c("LMNB1","CDKN1A","CDKN2A","CXCL8","GDF15","CCL2","CCND1")

SEN_RNA<-RNA%>%filter(genename%in%SEN)

##################################################### RNA vs. Protiesn

Protein<-read_xlsx("Volcano_2.xlsx")

Protein_RNA<-merge(Protein,RNA,by.x="hgnc_symbol",by.y="genename")

SenUnique<-read_xlsx("SenUnique_Proteins.xlsx")
Protein_RNA_Unique<-Protein_RNA%>%filter(hgnc_symbol%in%SenUnique$hgnc_symbol)



cor(Protein_RNA$log2FoldChange,Protein_RNA$FC2_EtopMPro,method = "pearson")
cor(Protein_RNA$log2FoldChange,Protein_RNA$FC2_EtopMPro,method = "spearman")

summary(lm(Protein_RNA$log2FoldChange~Protein_RNA$FC2_EtopMPro))

Protein_RNA_Outs<-Protein_RNA%>%filter(hgnc_symbol!="SERPINB2")%>%filter(hgnc_symbol!="ACTA1")

Protein_RNA%>%
  ggplot(aes(x=log2FoldChange,y=FC2_EtopMPro))+
  geom_point(shape=21)+stat_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(-5,5),ylim=c(-5,5))+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
                         axis.text.y = element_text(face="bold",size=24,color="black"),
                         axis.title.x = element_text(face="bold",size=20,color="black"),
                         axis.title.y = element_text(face="bold",size=20,color="black"),
                         axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
                         axis.ticks.length = unit(0.2,"cm"))+
  xlab("RNA Log2FC")+ylab("Protein Log2FC")


summary(lm(Protein_RNA$FC2_EtopMPro~Protein_RNA$log2FoldChange))




Protein_RNA_Unique%>%
  ggplot()+geom_boxplot(aes(group=sign(FC2_EtopMPro),x=FC2_EtopMPro,y=log2FoldChange),
                        fill=c("#B22222","#004D40"),alpha=0.6)+
  theme_classic()+
  scale_x_discrete()+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  coord_cartesian(ylim=c(-3,3))+
  scale_y_continuous(breaks=seq(-4,4,1))+ylab("RNA Expression Log2FC (Etop/Cyc)")+xlab("")

Protein_RNA_Unique%>%
  ggplot()+geom_violin(aes(group=sign(FC2_EtopMPro),x=FC2_EtopMPro,y=log2FoldChange),alpha=0.6)+
  theme_classic()+
  scale_x_discrete()+
  theme(axis.text.x = element_text(face="bold",size=24,color="black"),
        axis.text.y = element_text(face="bold",size=24,color="black"),
        axis.title.x = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  coord_cartesian(ylim=c(-3,3))+
  scale_y_continuous(breaks=seq(-4,4,1))+ylab("RNA Expression Log2FC (Etop/Cyc)")+xlab("")




cor(Protein_RNA_Unique$log2FoldChange,Protein_RNA_Unique$FC2_EtopMPro)^2
cor(Protein_RNA_Unique$log2FoldChange,Protein_RNA_Unique$FC2_EtopMPro,method="spearman")
plot(Protein_RNA$log2FoldChange,Protein_RNA$FC2_EtopMPro)
plot(Protein_RNA_Unique$log2FoldChange,Protein_RNA_Unique$FC2_EtopMPro)

RProteins_RNAProt<-Protein_RNA%>%filter(str_detect(Protein_RNA$Description,"ribosomal"))%>%
  filter(str_detect(hgnc_symbol,"^RPS|^RPL"))%>%filter(!str_detect(Description,"like"))

Factor_RNAProt<-Protein_RNA%>%filter(hgnc_symbol%in%c("EIF2S1","EEF1A1","EEF2",
                                                      "EIF5","EIF1"))


RProtein_SigSpli<-read_xlsx("EC_Ribo.xlsx")
RProtein_SigSpli_Prot<-RProteins_RNAProt%>%filter(hgnc_symbol%in%c("RPS19","RPL31","RPL28","RPS10",
                                                                   "RPS15A","RPS6","RPS14","RPL27A",
                                                                   "RPL38","RPL2","RPS17"))


ggplot()+geom_point(data=RProteins_RNAProt,aes(x=log2FoldChange,y=FC2_EtopMPro),color="red",shape=18,size=4,alpha=0.6)+
  geom_point(data=Factor_RNAProt,aes(x=log2FoldChange,y=FC2_EtopMPro),color="orange",size=4,alpha=0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),
        axis.text.y = element_text(size=24,color="black"),
        axis.title.x = element_text(size=20,color="black"),
        axis.title.y = element_text(size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  geom_hline(yintercept =0,color="blue", size=2,linetype=3)+
  geom_hline(yintercept =-0.2,color="black", size=2,linetype=3)+
  xlab("Log2FC RNA (Etop/Cyc)")+ylab("Log2FC Protein (Etop/Cyc)")+
  geom_text_repel(data=Factor_RNAProt,aes(x=log2FoldChange,y=FC2_EtopMPro,label=hgnc_symbol),
                   nudge_x=0,fontface="bold",size=3,min.segment.length = 0,max.overlaps=20)



RProtein_SigSpli%>%filter(Isoform%in%unlist(ExonSkip$Iso))

x<-ExonSkip[1,15]


############################################################################## Tg Analysis

ETg_RNA<-read_tsv("ETG_Etop.tsv")%>%filter(baseMean>1)
CTg_RNA<-read_tsv("CTG_Cyc.tsv")%>%filter(baseMean>1)

ETg_RNA$Norm<-ETg_RNA$log2FoldChange*log(ETg_RNA$baseMean,10)
CTg_RNA$Norm<-CTg_RNA$log2FoldChange*log(CTg_RNA$baseMean,10)

ETg_RNA%>%filter(log2FoldChange>=1|log2FoldChange<=-1)%>%filter(pvalue<=0.05)%>%nrow()
CTg_RNA%>%filter(log2FoldChange>=1|log2FoldChange<=-1)%>%filter(pvalue<=0.05)%>%nrow()

ETg_Short<-ETg_RNA%>%dplyr::select(1,8,3,6)%>%na.omit
names(ETg_Short)<-c("ENSG","hgnc_symbol","log2FC_ETg","pvalue")

CTg_Short<-CTg_RNA%>%dplyr::select(1,8,3,6)%>%na.omit
names(CTg_Short)<-c("ENSG","hgnc_symbol","log2FC_CTg","pvalue")



Tg_Compare<-merge(CTg_Short,ETg_Short,by="ENSG")
Tg_Compare_UPR%>%filter(log2FC_ETg>=1)%>%filter(pvalue.y<=0.05)%>%nrow()


head(ETg_RNA[order(ETg_RNA$Norm,decreasing = T),],10)
head(CTg_RNA[order(CTg_RNA$Norm,decreasing = T),],10)

ATF4<-read_xlsx("Han_ATF4_HFold.xlsx")
Xbp1_ATF6<-read_xlsx("Shoulders_ATF6_Xbp1.xlsx")
Xbp1<-Xbp1_ATF6%>%filter(`Fold-Change (XBP1s)`>2)
ATF6<-Xbp1_ATF6%>%filter(`Fold-Change ATF6`>2)



ETg_RNA_Up<-ETg_RNA%>%filter(log2FoldChange>=1)%>%filter(pvalue<=0.05)
CTg_RNA_Up<-CTg_RNA%>%filter(log2FoldChange>=1)%>%filter(pvalue<=0.05)

Tg_RNAs<-toupper(c(Xbp1$Symbol,ATF6$Symbol,ATF4$Symbol))%>%unique()

Reactome_UPR<-read_tsv("UPR_Reactome.tsv")

Reactome_UPR_2<-separate(Reactome_UPR,col = "MoleculeName",into = c("Junk","hgnc"),sep = " ")



Tg_Compare_UPR<-Tg_Compare%>%filter(hgnc_symbol.x%in%Tg_RNAs)


UIR_ETg<-ETg_RNA%>%filter(log2FoldChange>=1&genename%in%Tg_RNAs)%>%filter(pvalue<=0.05)
UIR_CTg<-CTg_RNA%>%filter(log2FoldChange>=1&genename%in%Tg_RNAs)%>%filter(pvalue<=0.05)


Xbp1_RNA<-Tg_Compare_UPR%>%filter(hgnc_symbol.x%in%Xbp1$Symbol)
ATF6_RNA<-Tg_Compare_UPR%>%filter(hgnc_symbol.x%in%ATF6$Symbol)
ATF4_RNA<-Tg_Compare_UPR%>%filter(hgnc_symbol.x%in%toupper(ATF4$Symbol))








CTg_UPR<-Tg_Compare_UPR%>%filter(log2FC_CTg>=1|log2FC_ETg>=1)


Tg_Compare_UPR%>%ggplot()+
  geom_point(aes(x=log2FC_CTg,y=log2FC_ETg),alpha=0.6,size=3,
             color=case_when(Tg_Compare_UPR$log2FC_CTg>=1 & Tg_Compare_UPR$log2FC_ETg <=1~"Firebrick",
                             Tg_Compare_UPR$log2FC_CTg>=1 & Tg_Compare_UPR$log2FC_ETg>=1~"steel blue",
                             Tg_Compare_UPR$log2FC_CTg<=1 & Tg_Compare_UPR$log2FC_ETg>=1~"#FFC107",
                             TRUE~"grey"))+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=20,color="black"),
        axis.text.y = element_text(face="bold",size=20,color="black"),
        axis.title.x = element_text(face="bold",size=18,color="black"),
        axis.title.y = element_text(face="bold",size=18,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  geom_hline(yintercept = c(1,0),color="black",size=1,linetype=2)+
  geom_vline(xintercept = c(1,0),color="black",size=1,linetype=2)+
  xlab("log2FC CTg/Cyc")+ylab("log2FC ETg/Etop")+
  coord_cartesian(xlim=c(-2,5),ylim=c(-2,5))+
  geom_text_repel(data=CTg_UPR,aes(x=log2FC_CTg,y=log2FC_ETg,label=hgnc_symbol.x,segment.size=0.7),
                  nudge_x=0,fontface="bold",size=3,min.segment.length = 0,max.overlaps = 20)




Tg_Compare_Red<-Tg_Compare_UPR%>%filter(Tg_Compare_UPR$log2FC_CTg>=1 & Tg_Compare_UPR$log2FC_ETg <=1)
Tg_Compare_Yellow<-Tg_Compare%>%filter(Tg_Compare$log2FC_CTg <1 & 
                                             Tg_Compare$log2FC_ETg >=1)





Tg_Compare_Blue<-Tg_Compare_UPR%>%filter(Tg_Compare_UPR$log2FC_CTg>=1 & Tg_Compare_UPR$log2FC_ETg>=1)
UpsetTg<-Tg_Compare_Red%>%
  mutate(
    Xbp1=ifelse(Tg_Compare_Red$hgnc_symbol.x%in%Xbp1_RNA$hgnc_symbol.x,1,0),
    Atf6=ifelse(Tg_Compare_Red$hgnc_symbol.x%in%ATF6_RNA$hgnc_symbol.x,1,0),
    Atf4=ifelse(Tg_Compare_Red$hgnc_symbol.x%in%ATF4_RNA$hgnc_symbol.x,1,0))

UpsetTg<-Tg_Compare_Blue%>%
  mutate(
    Xbp1=ifelse(Tg_Compare_Blue$hgnc_symbol.x%in%Xbp1_RNA$hgnc_symbol.x,1,0),
    Atf6=ifelse(Tg_Compare_Blue$hgnc_symbol.x%in%ATF6_RNA$hgnc_symbol.x,1,0),
    Atf4=ifelse(Tg_Compare_Blue$hgnc_symbol.x%in%ATF4_RNA$hgnc_symbol.x,1,0))

upset(UpsetTg[,8:10],keep.order = T,sets = c("Xbp1","Atf6","Atf4"))



  geom_label_repel(data=Compare_Tg_NormTen,aes(x=log2FoldChange.x,y=log2FoldChange.y,label=genename),
                   nudge_x=0,fontface="bold",fill="goldenrod1",size=3,min.segment.length = 0)


#################################################### Starvation




  E_SV<-read_tsv("SEN23_ETSV12_DESeq2.tab")%>%filter(baseMean>1)
  ESv_Short<-E_SV%>%dplyr::select(1,8,3,6)%>%na.omit
  names(ESv_Short)<-c("ENSG","hgnc_symbol","log2FC_ESv","pvalue")
  
  C_SV<-read_tsv("CSVA12_Cyc123_DESEQ2.tab")%>%filter(baseMean>1)
  CSv_Short<-C_SV%>%dplyr::select(1,8,3,6)%>%na.omit
  names(CSv_Short)<-c("ENSG","hgnc_symbol","log2FC_CSv","pvalue")

  
  
  ggplot()+
    geom_point(data=E_SV,
               aes(x=log2FoldChange,y=-log10(padj)),alpha=0.6,
               color=case_when(E_SV$log2FoldChange< -1 & E_SV$padj<0.05~"#B22222",
                               E_SV$log2FoldChange>1&E_SV$padj<0.05~"#004D40",
                               T~"black"))+
    geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
    xlab("Log2FC RNA (ESv/Etop)")+ ylab("-log10 p-value")+
    theme_classic()+
    theme(axis.text.x = element_text(size=24,color="black"),
          axis.text.y = element_text(size=24,color="black"),
          axis.title.x = element_text(size=20,color="black"),
          axis.title.y = element_text(size=20,color="black"),
          axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
          axis.ticks.length = unit(0.2,"cm"))
  
  
  ggplot()+
    geom_point(data=C_SV,
               aes(x=log2FoldChange,y=-log10(padj)),alpha=0.6,
               color=case_when(C_SV$log2FoldChange< -1 & C_SV$padj<0.05~"#B22222",
                               C_SV$log2FoldChange>1&C_SV$padj<0.05~"#004D40",
                               T~"black"))+
    geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
    xlab("Log2FC RNA (CSv/Cyc)")+ ylab("-log10 p-value")+
    theme_classic()+
    theme(axis.text.x = element_text(size=24,color="black"),
          axis.text.y = element_text(size=24,color="black"),
          axis.title.x = element_text(size=20,color="black"),
          axis.title.y = element_text(size=20,color="black"),
          axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
          axis.ticks.length = unit(0.2,"cm"))
  
  
  E_SV%>%filter(log2FoldChange>=1|log2FoldChange<=-1)%>%filter(pvalue<=0.05)%>%nrow()
  
  C_SV%>%filter(log2FoldChange>=1|log2FoldChange<=-1)%>%filter(pvalue<=0.05)%>%nrow()



  
  
  Starve_Genes_Murine<-c("ATF3","ATF5","EGR1","EGR2","EGR3","EGR4","FOS","FOSB","FOS1L",
                         "FOXQ1","GBX1","CUN","KLF6","MAFF","NR4A3","RAX","SRF","TCEAL7","VGLL3","ZMYND15",
                         "CCL20","CD244","CXCL11","CXCL10","EDN1","EDN2","IL11","LAG3","TFF2","TNFRS11B",
                         "TNFRSF25","THBS1", "CYR61","F3","FLRT3","HBEGF","HGF","ITGA1","JAM3","PLAU","PLAUR",
                         "CTGF","GDF9","SLITRK6","ATF4","WWTR1","DCT","AGR3","APCDD1","GRAMD1A","SP5",
                         "CDK6","ESRRG","TNFRSF19","ANGPT2","SORBS2","CACHD1","VAV3","CLCA1","NOTCH1")
  
  Tang<-read_xlsx("Tang_MCF7_AAR.xlsx")
  Tang_Genes<-Tang$`Gene Symbol`[Tang$`-Met`>1]
  
  AllStarveGenes<-c(Tang_Genes,Starve_Genes_Murine)


  
 Sv_Compare<-merge(CSv_Short,ESv_Short,by="ENSG")


 Sv_Compare_Starve<-Sv_Compare%>%filter(hgnc_symbol.x%in%AllStarveGenes)


 
 Sv_Compare_Starve%>%filter(log2FC_ESv>=1|log2FC_ESv<=-1 &pvalue.x<=0.05)%>%nrow()  
  
Sv_Compare_Starve_Genes<-Sv_Compare_Starve%>%filter(log2FC_CSv>=1)
 
 
  
 Sv_Compare%>%ggplot()+
   geom_point(aes(x=log2FC_CSv,y=log2FC_ESv),alpha=0.6,size=3,
              color=case_when(Sv_Compare$log2FC_CSv<=-1 & Sv_Compare$log2FC_ESv <=-1~"red",
                              Sv_Compare$log2FC_CSv>=1 & Sv_Compare$log2FC_ESv>=1~"blue",
                              TRUE~"black"))+
   theme_classic()+
   theme(axis.text.x = element_text(face="bold",size=20,color="black"),
         axis.text.y = element_text(face="bold",size=20,color="black"),
         axis.title.x = element_text(face="bold",size=18,color="black"),
         axis.title.y = element_text(face="bold",size=18,color="black"),
         axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
         axis.ticks.length = unit(0.2,"cm"))+
   coord_cartesian(xlim=c(-4,5))+
   geom_hline(yintercept = c(-1,1),color="black",size=0.5,linetype=2)+
   geom_vline(xintercept = c(-1,1),color="black",size=0.5,linetype=2)+
   xlab("log2FC CSv/Cyc")+ylab("log2FC ESv/Etop")
 
 
 
   geom_text_repel(data=Sv_Compare_Starve_Genes,aes(x=log2FC_CSv,y=log2FC_ESv,label=hgnc_symbol.x),
                   nudge_x=0,fontface="bold",size=3,min.segment.length = 0)

 
 ######################################################################################## Splicing
 
Col<-function(a,b,c,d){
  x<-sum(a,b)/2
  y<-sum(c,d)/2
  log(x/y,2)
}

Col(2,2,7,7)



Alt5<-read_tsv("Alt5_AllLibs.tsv")%>%filter(str_detect(feature_id,"inclusion"))
names(Alt5)<-c("ID","Coordinate","Sen_1","Cyc_1","Sen_2","Cyc_2","ETg_1","ETg_2","CSv_1","CSv_2",
               "ESv_1","ESv_2","CTg_1","CTg_2","Iso")


Alt5_FC<-Alt5%>%mutate(E_C=mapply(Col,Sen_1,Sen_2,Cyc_1,Cyc_2),
                          ETg_E=mapply(Col,ETg_1,ETg_2,Sen_1,Sen_2),
                          ESv_E=mapply(Col,ESv_1,ESv_2,Sen_1,Sen_2),
                          CTg_C=mapply(Col,CTg_1,CTg_2,Cyc_1,Cyc_2),
                          CSv_C=mapply(Col,CSv_1,CSv_2,Cyc_1,Cyc_2))%>%
  select(2,16,17,18,19,20)

Alt5_FC_Melt<-melt(as.data.frame(Alt5_FC))
Alt5_FC_Melt$variable<-factor(Alt5_FC_Melt$variable,levels=c("E_C","CTg_C","ETg_E","CSv_C","ESv_E"))

Alt5_FC_Melt%>%
  ggplot(aes(x=value,y=variable))+geom_boxplot()+coord_flip()+
  theme_classic()+
  theme(
    axis.text.y = element_text(face="bold",size=20,color="black"),
    axis.title.y = element_text(face="bold",size=18,color="black"),
    axis.line = element_line(size=1.5),axis.ticks.y = element_line(size=1.5),
    axis.ticks.length = unit(0.2,"cm"),
    axis.ticks.x=element_blank(),
    axis.text.x = element_text(face="bold",size=20,color="black"),
    axis.title.x = element_blank())+
  xlab("log2FC Alt. 5' SS Events")

ExonSkip<-read_tsv("ExonSkip_All_Libs.tsv")%>%filter(str_detect(feature_id,"inclusion"))
names(ExonSkip)<-c("ID","Coordinate","Sen_1","Cyc_1","Sen_2","Cyc_2","ETg_1","ETg_2","CSv_1","CSv_2",
               "ESv_1","ESv_2","CTg_1","CTg_2","Iso")
ExonSkip_FC<-ExonSkip%>%mutate(E_C=mapply(Col,Sen_1,Sen_2,Cyc_1,Cyc_2),
                       ETg_E=mapply(Col,ETg_1,ETg_2,Sen_1,Sen_2),
                       ESv_E=mapply(Col,ESv_1,ESv_2,Sen_1,Sen_2),
                       CTg_C=mapply(Col,CTg_1,CTg_2,Cyc_1,Cyc_2),
                       CSv_C=mapply(Col,CSv_1,CSv_2,Cyc_1,Cyc_2))%>%
  select(2,16,17,18,19,20)

ExonSkip_FC_Melt<-melt(as.data.frame(ExonSkip_FC))
ExonSkip_FC_Melt$variable<-factor(ExonSkip_FC_Melt$variable,levels=c("E_C","CTg_C","ETg_E","CSv_C","ESv_E"))

ExonSkip_FC_Melt%>%
  ggplot(aes(x=value,y=variable))+geom_boxplot()+coord_flip()+
  theme_classic()+
  theme(
        axis.text.y = element_text(face="bold",size=20,color="black"),
        axis.title.y = element_text(face="bold",size=18,color="black"),
        axis.line = element_line(size=1.5),axis.ticks.y = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(face="bold",size=20,color="black"),
        axis.title.x = element_blank())+
  xlab("log2FC Exon Skip Events")



Alt3<-read_tsv("Alt3_AllLibs.tsv")%>%filter(str_detect(feature_id,"inclusion"))
names(Alt3)<-c("ID","Coordinate","Sen_1","Cyc_1","Sen_2","Cyc_2","ETg_1","ETg_2","CSv_1","CSv_2",
               "ESv_1","ESv_2","CTg_1","CTg_2","Iso")


Alt3_FC<-Alt3%>%mutate(E_C=mapply(Col,Sen_1,Sen_2,Cyc_1,Cyc_2),
                       ETg_E=mapply(Col,ETg_1,ETg_2,Sen_1,Sen_2),
                       ESv_E=mapply(Col,ESv_1,ESv_2,Sen_1,Sen_2),
                       CTg_C=mapply(Col,CTg_1,CTg_2,Cyc_1,Cyc_2),
                       CSv_C=mapply(Col,CSv_1,CSv_2,Cyc_1,Cyc_2))%>%
  select(2,16,17,18,19,20)

Alt3_FC_Melt<-melt(as.data.frame(Alt3_FC))
Alt3_FC_Melt$variable<-factor(Alt3_FC_Melt$variable,levels=c("E_C","CTg_C","ETg_E","CSv_C","ESv_E"))

Alt3_FC_Melt%>%
  ggplot(aes(x=value,y=variable))+geom_boxplot()+coord_flip()+
  theme_classic()+
  theme(
    axis.text.y = element_text(face="bold",size=20,color="black"),
    axis.title.y = element_text(face="bold",size=18,color="black"),
    axis.line = element_line(size=1.5),axis.ticks.y = element_line(size=1.5),
    axis.ticks.length = unit(0.2,"cm"),
    axis.ticks.x=element_blank(),
    axis.text.x = element_text(face="bold",size=20,color="black"),
    axis.title.x = element_blank())+
  xlab("log2FC Alt. 3' SS Events")












ExonSkip<-read_tsv("ExonSkip_All_Libs.tsv")%>%filter(str_detect(feature_id,"inclusion"))
IntronRet<-read_tsv("IntronRet_AllLibs.tsv")%>%filter(str_detect(feature_id,"inclusion")) 







DiffSplice_All<-read_xlsx("DiffSplice_AllComparisons.xlsx")

DiffSplice_CTg<-CTg_Short%>%filter(ENSG%in%DiffSplice_All$`CTg/Cyc`)

DiffSplice_CSv<-CSv_Short%>%filter(ENSG%in%DiffSplice_All$`CSv/Cyc`)

DiffSplice_EC<-RNA_All%>%filter(ENSG%in%DiffSplice_All$`E/C`)
write_xlsx(DiffSplice_EC,"DiffSplice_EC.xlsx")


All_DiffSplice_RProtein<-as.data.frame(c(DiffSplice_CSv$hgnc_symbol,DiffSplice_CTg$hgnc_symbol,
                           DiffSplice_EC$hgnc_symbol),
                           c(DiffSplice_CSv$ENSG,DiffSplice_CTg$ENSG,
                             DiffSplice_EC$ENSG))%>%unique()

All_DiffSplice_RProtein$ENSG<-rownames(All_DiffSplice_RProtein)

write_xlsx(All_DiffSplice_RProtein,"DiffSplice_Proteins.xlsx")




.RNA_DiffSplice_Rprotein<-RNA%>%filter()

########################################################################################## NanoPlen

NanoPlen<-read_tsv("FourSample_TLength.tsv",col_names = T)
names(NanoPlen)<-c("Row","Sample","ENST","MLength","Log2_Length")

NanoPlen_S<-pivot_wider(NanoPlen,id_cols = "ENST",names_from = "Sample",values_from = "MLength")%>%na.omit()
 
values<-NanoPlen_S$ENST
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
x<-getBM(attributes = c('ensembl_transcript_id', "description","hgnc_symbol","gene_biotype"),
         filters = 'ensembl_transcript_id',
         values = values, 
         mart = ensembl)
NanoPlen_N<-merge(NanoPlen_S,x,by.x="ENST",by.y="ensembl_transcript_id")



NanoPlen_N<-NanoPlen_N%>%mutate(SC=log(Sen23/Cyc23,2),SSV=log(ETSV12/Sen23,2),
                                CSV=log(CSVA12/Cyc23,2), ESVC=log(ETSV12/Cyc23),
                                SSV_Diff=ETSV12-Sen23)

Nanoplen_RProteinsLong<-NanoPlen_N%>%filter(hgnc_symbol%in%RProteins$hgnc_symbol)

Nanoplen_RProteins<-NanoPlen_WideN%>%filter(hgnc_symbol%in%RProteins$hgnc_symbol)

NanoPlen_WideN<-merge(NanoPlen,x,by.x="ENST",by.y="ensembl_transcript_id")


NanoPlen_WideN_Secreted<-NanoPlen_WideN%>%filter(hgnc_symbol%in%Secreted$Gene)

NanoPlen_N_Secreted<-NanoPlen_N%>%filter(hgnc_symbol%in%Secreted$Gene)

Nanoplen_RProteins%>%filter(Sample%in%c("Sen23","Cyc23"))%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  scale_x_continuous(breaks = c(seq(0,2000,200)))+coord_cartesian(xlim = c(200,1400))
  
NanoPlen_WideN_Secreted%>%filter(Sample%in%c("Sen23","Cyc23"))%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)

NanoPlen_WideN_Secreted%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,4000))

NanoPlen_WideN_Secreted%>%filter(Sample%in%c("CSVA12","ETSV12"))%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,4000))

NanoPlen_WideN%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,4000))


EC_Transcripts<-read_tsv("DESeq_Sen_Cyc_ENST.tab")%>%filter(baseMean>=5)

NanoPlen_F<-NanoPlen%>%filter(ENST%in%EC_Transcripts$tid)
NanoPlen_Widen_F<-merge(NanoPlen_F,x,by.x="ENST",by.y="ensembl_transcript_id")
NanoPlen_Widen_F%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,4000))

NanoPlen_WideN%>%filter(hgnc_symbol%in%Stress_Up$hgnc_symbol.x)%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,4000))


NanoPlen_WideN%>%filter(hgnc_symbol%in%Protein_RNA$hgnc_symbol[Protein_RNA$FC2_EtopMPro<=-0.2])%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,4000))

NanoPlen_WideN%>%filter(hgnc_symbol%in%Protein_RNA$hgnc_symbol[Protein_RNA$FC2_EtopMPro<=-0.4])%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,4000))


NanoPlen_WideN%>%filter(hgnc_symbol%in%Protein_RNA$hgnc_symbol)%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,4000))

Protein_RNA$hgnc_symbol[Protein_RNA$FC2_EtopMPro>=0.4]

Protein_RNA$hgnc_symbol[Protein_RNA$FC2_EtopMPro<=-0.4]


NanoPlen_WideN_P<-NanoPlen_WideN%>%filter(Sample%in%c("Sen23","Cyc23"))%>%
  mutate(Prot=case_when(hgnc_symbol%in%Protein_RNA$hgnc_symbol[Protein_RNA$FC2_EtopMPro<=-0.2]~"x < -0.2",
                        hgnc_symbol%in%Protein_RNA$hgnc_symbol[between(Protein_RNA$FC2_EtopMPro,-0.19,0.19)]~"-0.2 < x < 0.2",
                        hgnc_symbol%in%Protein_RNA$hgnc_symbol[Protein_RNA$FC2_EtopMPro>=0.2]~"0.2 < x"))%>%na.omit()

NanoPlen_WideN_P$Prot<-factor(NanoPlen_WideN_P$Prot,levels=c("x < -0.2","-0.2 < x < 0.2","0.2 < x"))


NanoPlen_WideN_P%>%
  filter(hgnc_symbol%in%Protein_RNA$hgnc_symbol)%>%
  ggplot(aes(x=MLength,color=Sample))+stat_ecdf(geom="step")+geom_hline(yintercept = 0.5)+
  coord_cartesian(xlim=c(0,3000))+facet_grid(cols = vars(Prot))


NanoPlen_WideN_P%>%
  filter(hgnc_symbol%in%Protein_RNA$hgnc_symbol)%>%
  ggplot(aes(x=MLength,color=Sample))+geom_boxplot()+facet_grid(rows = vars(Prot))

Nanoplen_RProteins%>%filter(Sample%in%c("Sen23","Cyc23"))%>%
  filter(hgnc_symbol%in%Protein_RNA$hgnc_symbol)%>%
  ggplot(aes(x=MLength,color=Sample))+geom_boxplot()


Nanoplen_RProteins%>%filter(hgnc_symbol=="RPS10")

wilcox.test()
################################################################### RNA All
 
 
 
 RNA<-read_tsv("CYC_SEN_DESEQ2.tab")%>%filter(baseMean>1)
 RNA_Short<-RNA%>%dplyr::select(1,8,3,4,6)%>%na.omit
 names(RNA_Short)<-c("ENSG","hgnc_symbol","log2FC_EC","lfcSE","pvalue")
 
 ESvC<-read_tsv("DESeq2_ESv_C.tsv")%>%filter(baseMean>1)
ESv_C_Short<-  ESvC%>%dplyr::select(1,8,3,4,6)%>%na.omit
 names(ESv_C_Short)<-c("ENSG","hgnc_symbol","log2FC_ESvC","lfcSE","pvalue")

ETg_C<-read_tsv("DESeq2_ETg_C.tsv")%>%filter(baseMean>1)
ETg_C_Short<-ETg_C%>%dplyr::select(1,8,3,4,6)%>%na.omit
names(ETg_C_Short)<-c("ENSG","hgnc_symbol","log2FC_ETgC","lfcSE","pvalue")
 
 
   
RNA_All<-merge(RNA_Short,ESv_C_Short,by="ENSG")%>%merge(.,ETg_C_Short,by="ENSG")
RNA_All_Top<-RNA_All%>%filter(hgnc_symbol.x%in%topnormal)


RNA_All%>%ggplot()+
  geom_point(aes(x=log2FC_EC,y=log2FC_ESvC),alpha=0.6,
             color=case_when(RNA_All$log2FC_EC<=-1 & RNA_All$log2FC_ESvC <=-1~"red",
                             RNA_All$log2FC_EC>=1 & RNA_All$log2FC_ESvC>=1~"blue",
                             TRUE~"grey70"))+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=20,color="black"),
        axis.text.y = element_text(face="bold",size=20,color="black"),
        axis.title.x = element_text(face="bold",size=18,color="black"),
        axis.title.y = element_text(face="bold",size=18,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  geom_hline(yintercept = c(-1,1),color="black",size=1,linetype=2)+
  geom_vline(xintercept = c(-1,1),color="black",size=1,linetype=2)+
  xlab("log2FC E/C")+ylab("log2FC ESv/C")+
  geom_label_repel(data=RNA_All_Top,aes(x=log2FC_EC,y=log2FC_ESvC,label=hgnc_symbol.x),
                  nudge_x=0,fontface="bold",fill="goldenrod1",size=3,min.segment.length = 0)


RNA_All%>%ggplot()+
  geom_point(aes(x=log2FC_EC,y=log2FC_ETgC),alpha=0.6,
             color=case_when(RNA_All$log2FC_EC<=-1 & RNA_All$log2FC_ETgC <=-1~"red",
                             RNA_All$log2FC_EC>=1 & RNA_All$log2FC_ETgC>=1~"blue",
                             TRUE~"grey70"))+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=20,color="black"),
        axis.text.y = element_text(face="bold",size=20,color="black"),
        axis.title.x = element_text(face="bold",size=18,color="black"),
        axis.title.y = element_text(face="bold",size=18,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  geom_hline(yintercept = c(-1,1),color="black",size=1,linetype=2)+
  geom_vline(xintercept = c(-1,1),color="black",size=1,linetype=2)+
  xlab("log2FC E/C")+ylab("log2FC ETg/C")+
  geom_label_repel(data=RNA_All_Top,aes(x=log2FC_EC,y=log2FC_ETgC,label=hgnc_symbol.x),
                   nudge_x=0,fontface="bold",fill="goldenrod1",size=3,min.segment.length = 0)




Stress_Up<-RNA_All%>%filter(log2FC_ESvC>1.5*log2FC_EC)%>%
  filter(log2FC_ETgC>1.5*log2FC_EC)%>%filter(pvalue.y<=0.05&pvalue<=0.05)%>%
  filter(log2FC_ESvC>0&log2FC_ETgC>0)
  

Stress_Up_Secreted<-Stress_Up%>%filter(hgnc_symbol.x%in%Secreted$Gene)%>%
  filter(log2FC_ESvC>1&log2FC_ETgC>1)


Stress_Up_SASP<-Stress_Up%>%filter(hgnc_symbol.x%in%SASP)%>%
  filter(log2FC_ESvC>1&log2FC_ETgC>1)

Stress_Up_Basisty<-Stress_Up%>%filter(hgnc_symbol.x%in%Basisty_SASP$Genes)%>%
  filter(log2FC_ESvC>1&log2FC_ETgC>1)

write_xlsx(Stress_Up,"Stress_Up.xlsx")

write_xlsx(RNA_All,"RNA_All.xlsx")

write_xlsx(Stress_Up_Tg,"Stress_Up_Tg.xlsx")


SASP<-c("IL6","CXCL8","CCL2","CCL8","CCL7","CCL13","IGFBP3","CSF2","CXCL1",
                "CXCL2","CXCL3","ICAM1","GDF15","SERPINE1", "STC1",
        "TIMP1","TIMP2","SERPINEB2","IL1A","IL1B","IL13","IL15",
        "CXCL2","CXCL3","CCL8","CCL13","CCL3","CCL20","CCL16","CCL11",
        "CCL25","MIF","AREG","EREG","EGF","FGF2","HGF","FGF7","VEGFA",
        "ANG","KITLG","CXCL12","PGF","IGFBP2","IGFBP3","IGFBP4","IGFPB6","IGFPB7",
        "MMP1","MMP3","MMP10","MMP12","MMP13","MMP14","ICAM1","ICAM3",
        "TNFRSF11B","TNFRSF15A","TNFRSF10C","FAS","TNFRSF1A","PLAUR","IL6ST",
        "EGF")


ETg_Short%>%filter(hgnc_symbol%in%SASP)%>%filter(log2FC_ETg>=1)

ETg_Short%>%filter(log2FC_ETg>=1)%>%nrow()


RNA_All_SASP<-RNA_All%>%filter(hgnc_symbol.x%in%SASP)%>%filter(log2FC_EC>=1)%>%
  filter(pvalue.x<=0.05)

sum(Stress_Up_SASP$hgnc_symbol.x%in%RNA_All_SASP$hgnc_symbol.x)

SASP_Compare<-RNA_All%>%filter(hgnc_symbol.x%in%c(Stress_Up_SASP$hgnc_symbol.x,RNA_All_SASP$hgnc_symbol.x))
write_xlsx(SASP_Compare,"SASP_Compare.xlsx")




Basisty<-read_xlsx("Basisty_IR.xlsx")

Basisty_SASP<-Basisty%>%filter(`AVG Log2 Ratio`>1)



Secreted<-read_tsv("Secreted_Proteins.tsv")


RNA_All_Basisty<-RNA_All%>%filter(hgnc_symbol.x%in%Basisty_SASP$Genes)%>%
  filter(log2FC_EC>1|log2FC_ESvC>1|log2FC_ETgC>1)

RNA_All_Secreted<-RNA_All%>%filter(hgnc_symbol.x%in%Secreted$Gene)%>%
  filter(log2FC_EC>1|log2FC_ESvC>1|log2FC_ETgC>1)


write_xlsx(RNA_All_Secreted,"RNA_All_Secreted.xlsx")


Stress_RNA_Secreted<-RNA_All_Secreted%>%filter(log2FC_ESvC>=1.5*log2FC_EC&log2FC_ETgC>=1.5*log2FC_EC)%>%
  filter(log2FC_ESvC>0&log2FC_ETgC>=0)
Stress_RNA_Basisty<-RNA_All_Basisty%>%filter(log2FC_ESvC>=1.5*log2FC_EC&log2FC_ETgC>=1.5*log2FC_EC)%>%
  filter(log2FC_ESvC>0&log2FC_ETgC>=0)


RNA_All_Secreted%>%filter(log2FC_EC>1&log2FC_ESvC<1&log2FC_ETgC<1)%>%nrow()
RNA_All_Secreted%>%filter(log2FC_EC<1&log2FC_ESvC<1&log2FC_ETgC>1)%>%nrow()
RNA_All_Secreted%>%filter(log2FC_EC>1&log2FC_ESvC>1&log2FC_ETgC>1)%>%nrow()
RNA_All_Secreted%>%filter(log2FC_EC<1&log2FC_ESvC<1&log2FC_ETgC>1)%>%nrow()

library(UpSetR)
UpsetSecreted<-RNA_All_Secreted%>%
  mutate(
    StressSen=ifelse(log2FC_ESvC>=1.5*log2FC_EC&log2FC_ETgC>=1.5*log2FC_EC,1,0),
    Common=ifelse(log2FC_ESvC>=1&log2FC_EC>=1&log2FC_ETgC>=1,1,0),
    EC=ifelse(log2FC_EC>=1,1,0),
    ESv=ifelse(log2FC_ESvC>=1|log2FC_ESvC>=1.5*log2FC_EC,1,0),
    ETg=ifelse(log2FC_ETgC>=1|log2FC_ETgC>=1.5*log2FC_EC,1,0))






upset(UpsetSecreted[,14:18],keep.order = T,order.by = "freq")




write_xlsx(RNA_All_Secreted,"RNA_All_Secreted.xlsx")

install.packages("reshape")
library(reshape)

RNA_All_Secreted_Short<-RNA_All_Secreted%>%select(2,3,7,11)

meltData <- melt(RNA_All_Secreted_Short)
boxplot(data=meltData, value~variable)

meltData_Stress<-meltData%>%filter(hgnc_symbol.x%in%c("HBEGF","LIF","NAMPT","FGF2","ADAMTS19"))



meltData%>%group_by(variable)%>%ggplot+geom_boxplot(aes(x=variable,y=value))+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",size=14,color="black"),
        axis.text.y = element_text(face="bold",size=14,color="black"),
        axis.title.x = element_text(face="bold",size=18,color="black"),
        axis.title.y = element_text(face="bold",size=18,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  geom_hline(yintercept = c(0,1),color="black",size=1,linetype=2)+
  xlab("Comparison")+ylab("log2FC")+
  geom_label_repel(data=meltData_Stress,aes(x=variable,y=value,label=hgnc_symbol.x),
                   nudge_x=0,fontface="bold",fill="goldenrod1",size=3,min.segment.length = 0)


mean(RNA_All_Secreted$log2FC_EC)
mean(RNA_All_Secreted$log2FC_ESvC)
mean(RNA_All_Secreted$log2FC_ETgC)


t.test(RNA_All_Secreted$log2FC_EC,RNA_All_Secreted$log2FC_ESvC)

t.test(RNA_All_Secreted$log2FC_EC,RNA_All_Secreted$log2FC_ETgC)




5############################################################# Spliceosome RNAs
2.8/2.1
 
 
 
 
 
Spliceosome<-c("SNRNP40","WBP11" ,"CHERP","RBM22","SF3B2","FUS","HNRNPC","CDC5L","HNRNPM","SF3A1","SNRPD3","SNU13",
               "PHF5A","SNW1","SRSF5","ACIN1","CRNKL1","PQBP1","SF3A2","LSM5","DDX5",
               "EFTUD2","PRPF19","SRSF3","TCERG1","NCBP2","SF3B6","SF3B1","SRSF7","BCAS2",
               "PRPF3","RBM25","SRSF6","SNRPC","SNRPD2","SNRPB2","LSM8","LSM4","SNRPA1",
               "CTNNBL1",  "PRPF38B",  "PRPF38A",  "HNRNPA1",  "SRSF1",  "TRA2B",  "PRPF4",  "NCBP1",  "PPIL1",
               "SNRPF",  "EIF4A3",  "SNRNP200",  "DDX46",  "RBMX",  "CWC15",  "HNRNPU",  "SRSF2",  "MAGOH",
               "U2SURP",  "LSM6",  "HNRNPK",  "SNRPD1",  "CDC40","SF3B5",  "HNRNPA3",
               "LSM3",  "PLRG1",  "PPIH",  "SART1",  "SF3A3",  "SRSF10",  "SF3B3",  "LSM2", "DDX39B",
               "ISY1",  "RBM8A",  "PRPF8")

RNA_Spliceosome<-RNA%>%filter(genename%in%Spliceosome)
median(RNA_Spliceosome$log2FoldChange)
median(RNA_Spliceosome$Norm)


DiffSplice_ES<-read_xlsx("DiffSplice_Starve.xlsx",sheet="common-ETSV-SEN")


DiffSPlice_RNA_ESV<-E_SV%>%filter(Gene%in%DiffSplice_ES$`Common-89 @ 0.01`)



############################### ER increase

ER_Increased<-c("CKAP4","CANX", "S100B","DNAJB1","MAN1B1","EDEM1","BAK1","BAX",
                "ERO1","ERLEC1","BCAP31","TICAM2","DERL1","VCP","NPL4","EIF2Ak3","WFS1",
                "UBE2J1","UBE2J2","UBE2G2","STUB1","SEL1L","SYVN1","RNF5")

Protein_RNA_ER<-Protein_RNA%>%filter(hgnc_symbol%in%ER_Increased)


Protein_RNA_ER%>%ggplot(aes(x=log2FoldChange,y=FC2_EtopMPro))+
  geom_point(color=ifelse(Protein_RNA_ER$padj<=0.05,"red","black"))+theme_classic()+
  stat_smooth(method="lm")

Protein_RNA%>%ggplot(aes(x=log2FoldChange,y=FC2_EtopMPro))+
  geom_point(color=ifelse(Protein_RNA$padj<=0.05,"red","black"))+theme_classic()+
  stat_smooth(method="lm")

lm(Protein_RNA_ER$FC2_EtopMPro~Protein_RNA_ER$log2FoldChange)
cor(x=Protein_RNA_ER$log2FoldChange,y=Protein_RNA_ER$FC2_EtopMPro)^2
cor(x=Protein_RNA_ER$log2FoldChange,y=Protein_RNA_ER$FC2_EtopMPro,method="spearman")


cor(x=Protein_RNA$log2FoldChange,y=Protein_RNA$FC2_EtopMPro)^2



