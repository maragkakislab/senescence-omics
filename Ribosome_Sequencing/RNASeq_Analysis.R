


names(IR_rna)[names(IR_rna)=="log2FoldChange"]<-"log2FC_IR"
names(EC_rna)[names(EC_rna)=="log2FC_RNA"]<-"log2FC_E"
names(QC_rna)[names(QC_rna)=="log2FC_RNA"]<-"log2FC_Q"

IR_rna%>%ggplot()+geom_point(aes(x=log2FC_IR,y=-log10(padj)),
                             color=case_when(IR_rna$log2FC_IR>=1 &IR_rna$padj<=0.05~"blue",
                                             IR_rna$log2FC_IR<= -1 & IR_rna$padj <= 0.05~"red",
                                             T~"black"))+
  theme_classic()

EC_rna%>%ggplot()+geom_point(aes(x=log2FC_E,y=-log10(padj)),
                             color=case_when(EC_rna$log2FC_E>=1 &EC_rna$padj<=0.05~"blue",
                                             EC_rna$log2FC_E<= -1 & EC_rna$padj <= 0.05~"red",
                                             T~"black"))+
  theme_classic()



All_RNA<-merge(IR_rna,EC_rna,by="hgnc_symbol")%>%merge(.,QC_rna,by="hgnc_symbol")%>%
  na.omit()

All_RNA_Signal<-All_RNA%>%filter(hgnc_symbol%in%SignalP$Gene)



plot(All_RNA$log2FC_IR,All_RNA$log2FC_E)


plot(All_RNA$log2FC_IR,All_RNA$log2FC_Q)

ER_Increased<-c("CKAP4","CANX", "S100B","DNAJB1","MAN1B1","EDEM1","BAK1","BAX",
                "ERO1","ERLEC1","BCAP31","TICAM2","DERL1","VCP","NPL4","EIF2Ak3","WFS1",
                "UBE2J1","UBE2J2","UBE2G2","STUB1","SEL1L","SYVN1","RNF5")

All_RNA%>%filter(hgnc_symbol%in%ER_Increased)

SASP<-c("IL6","CXCL8","CCL2","CCL8","CCL7","CCL13","IGFBP3","CSF2","CXCL1",
        "CXCL2","CXCL3","ICAM1","GDF15","SERPINE1", "STC1",
        "TIMP1","TIMP2","SERPINEB2","IL1A","IL1B","IL13","IL15",
        "CXCL2","CXCL3","CCL8","CCL13","CCL3","CCL20","CCL16","CCL11",
        "CCL25","MIF","AREG","EREG","EGF","FGF2","HGF","FGF7","VEGFA",
        "ANG","KITLG","CXCL12","PGF","IGFBP2","IGFBP3","IGFBP4","IGFPB6","IGFPB7",
        "MMP1","MMP3","MMP10","MMP12","MMP13","MMP14","ICAM1","ICAM3",
        "TNFRSF11B","TNFRSF15A","TNFRSF10C","FAS","TNFRSF1A","PLAUR","IL6ST",
        "EGF")

All_RNA%>%filter(hgnc_symbol%in%SASP)


############################## Cluster profiler

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(clusterProfiler)
gene_list<-All_RNA$log2FC_Q
names(gene_list)<-All_RNA$hgnc_symbol
gene_list<-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

head(gene_list)

gse_bp <- gseGO(geneList=gene_list, 
                ont ="BP", 
                keyType = "SYMBOL",
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Hs.eg.db)

gse_bp_simple<-clusterProfiler::simplify(gse_bp)
dotplot(gse_bp,showCategory=5,split=".sign")+theme(axis.text.y=element_text(size=12),
                                                          strip.text.x = element_text(size=18,face="bold"),
                                                          axis.text.x = element_text(size = 10,face="bold"),
                                                          axis.title.x =element_text(size=16,face="bold") )+
  scale_y_discrete(labels=label_wrap_gen(35))+
  facet_grid(.~.sign)



