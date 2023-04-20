library(DESeq2)
library(tidyverse)
library(readxl)
library(writexl)
library(ggrepel)
library(biomaRt)


All_Counts<-merge(RNA_Counts,Ribo_Counts)

All_Counts_Matrix<-All_Counts[,2:23]
rownames(All_Counts_Matrix)<-All_Counts$ENSG
coldata<-read_xlsx("ColData_Trim.xlsx")
ddsMat<-DESeqDataSetFromMatrix(All_Counts_Matrix,
                               colData = coldata,
                               design=~Condition+SeqType+Condition:SeqType)

ddsMat$SeqType=relevel(ddsMat$SeqType,"RNA")

ddsMat<-DESeq(ddsMat)


res_EC<-as.data.frame(results(ddsMat,contrast=list("ConditionEtoposide.SeqTypeRIBO")))
res_EC$ENSG<-rownames(res_EC)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = res_EC$ENSG, 
          mart = ensembl)


res_EC<-merge(ID,res_EC,by.x="ensembl_gene_id",by.y="ENSG")
res_EC<-res_EC%>%na.omit()%>%filter(gene_biotype=="protein_coding")
writexl::write_xlsx(res_EC,"Riboseq_EC.xlsx")



res_EC_sig<-res_EC%>%filter(padj<=0.05)%>%filter(log2FoldChange>=1|log2FoldChange<= -1)

res_EC%>%ggplot()+geom_point(aes(x=log2FoldChange,y=-log10(padj)),
                             color=case_when(res_EC$log2FoldChange>=1 &res_EC$padj<=0.05~"blue",
                                             res_EC$log2FoldChange<= -1 & res_EC$padj <= 0.05~"red",
                             T~"black"))+
  theme_classic()+
  geom_text_repel(data=res_EC_sig,aes(x=log2FoldChange,y=-log10(padj),label=hgnc_symbol),size=3,
                  min.segment.length = 0)


res_EC%>%filter(padj<=0.05)%>%nrow()
writexl::write_xlsx(res_EC_sig,"res_EC_sig.xlsx")
3
#########################
res_IC<-as.data.frame(results(ddsMat,contrast=list("ConditionIonizing.SeqTypeRIBO")))
res_IC$ENSG<-rownames(res_IC)

ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = res_IC$ENSG, 
          mart = ensembl)

res_IC<-merge(ID,res_IC,by.x="ensembl_gene_id",by.y="ENSG")
res_IC<-res_IC%>%na.omit()%>%filter(gene_biotype=="protein_coding")
writexl::write_xlsx(res_IC,"Riboseq_IC.xlsx")

res_IC_sig<-res_IC%>%filter(padj<=0.05)%>%filter(log2FoldChange>=1|log2FoldChange<= -1)
write_xlsx(res_IC_sig,"res_IC_sig.xlsx")




res_IC%>%ggplot()+geom_point(aes(x=log2FoldChange,y=-log10(padj)),shape=21,
                             color=case_when(res_IC$log2FoldChange>=1 &res_IC$padj<=0.05~"blue",
                                             res_IC$log2FoldChange<= -1 & res_IC$padj <= 0.05~"red",
                                             T~"black"))+
  theme_classic()+
  geom_text_repel(data=res_IC_sig,aes(x=log2FoldChange,y=-log10(padj),label=hgnc_symbol),size=3,
                  max.overlaps = 25)


res_IC_sig<-res_IC%>%na.omit()%>%filter(gene_biotype=="protein_coding")

res_IC%>%filter(padj<=0.05)%>%nrow()

###############################

res_QC<-as.data.frame(results(ddsMat,contrast=list("ConditionQuiescent.SeqTypeRIBO")))
res_QC$ENSG<-rownames(res_QC)

ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = res_QC$ENSG, 
          mart = ensembl)

res_QC<-merge(ID,res_QC,by.x="ensembl_gene_id",by.y="ENSG")
res_QC<-res_QC%>%na.omit()%>%filter(gene_biotype=="protein_coding")
writexl::write_xlsx(res_QC,"Riboseq_QC.xlsx")

res_QC_sig<-res_QC%>%filter(padj<=0.05)%>%filter(log2FoldChange>=1|log2FoldChange<= -1)
write_xlsx(res_QC_sig,"res_QC_sig.xlsx")



res_QC%>%ggplot()+geom_point(aes(x=log2FoldChange,y=-log10(padj)),shape=21,
                             color=case_when(res_QC$log2FoldChange>=1 &res_QC$padj<=0.05~"blue",
                                             res_QC$log2FoldChange<= -1 & res_QC$padj <= 0.05~"red",
                                             T~"black"))+
  theme_classic()+
  geom_text_repel(data=res_QC_sig,aes(x=log2FoldChange,y=-log10(padj),label=hgnc_symbol),size=3,
                  max.overlaps = 20)




res_QC%>%filter(padj<=0.05)%>%nrow()


names(res_QC)[names(res_QC)=="log2FoldChange"]<-"log2FC_QC"
names(res_IC)[names(res_IC)=="log2FoldChange"]<-"log2FC_IC"
names(res_EC)[names(res_EC)=="log2FoldChange"]<-"log2FC_EC"



All_Res<-merge(res_EC,res_IC,by="ENSG")%>%merge(.,res_QC,by="ENSG")%>%na.omit()
data.frame(colnames(All_Res))


plot(All_Res$log2FC_EC,All_Res$log2FC_QC)


All_Res_Longer<-pivot_longer(All_Res,cols=c(3,9,15),names_to = "Comparison")%>%
  select(17,18)


All_Res_Longer$Comparison<-factor(All_Res_Longer$Comparison,levels = c("log2FC_EC","log2FC_IC","log2FC_QC"))

All_Res_Longer%>%ggplot()+geom_violin(aes(x=Comparison,y=value,fill=Comparison))

All_Res_Sen<-All_Res%>%filter(between(log2FC_QC,-1,1))%>%
  filter(!between(log2FC_EC,-1,1)&!between(log2FC_IC,-1,1))%>%
  filter(padj.x<=0.05&padj.y<=0.05)

All_Res_Sig<-All_Res%>%filter(padj.x<=0.05 | padj.y<=0.05 | padj <=0.05)


write_xlsx(All_Res_Sig,"All_NormalizedRiboseq.xlsx")
#####################################################




plot(All_RNA$log2FC_IR,All_RNA$log2FC_E)



################################################################################################################
ind = which(coldata$SeqType == "RIBO")
coldata_ribo = coldata[ind,]

rownames(Ribo_Counts)<-Ribo_Counts$ENSG
rownames(RNA_Counts)<-RNA_Counts$ENSG

ddsMat_ribo <- DESeqDataSetFromMatrix(countData = Ribo_Counts[,-1],
                                      colData = coldata_ribo, design =~ Condition)

levels(ddsMat_ribo$Condition)
ddsMat_ribo <- DESeq(ddsMat_ribo)

ind = which(coldata$SeqType == "RNA")
coldata_rna = coldata[ind,]
ddsMat_rna <- DESeqDataSetFromMatrix(countData = RNA_Counts[,-1],
                                     colData = coldata_rna, design =~ Condition)
ddsMat_rna<- DESeq(ddsMat_rna)

################################################################ EC
EC_ribo <- as.data.frame(results(ddsMat_ribo, contrast=list("Condition_Etoposide_vs_Cycling")))
EC_rna <- as.data.frame(results(ddsMat_rna, contrast=list("Condition_Etoposide_vs_Cycling")))

EC_ribo$ENSG<-rownames(EC_ribo)
EC_rna$ENSG<-rownames(EC_rna)


ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = EC_ribo$ENSG, 
          mart = ensembl)

EC_ribo<-merge(EC_ribo,ID,by.x="ENSG",by.y="ensembl_gene_id")%>%
  filter(gene_biotype=="protein_coding")%>%filter(hgnc_symbol!='')

ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = EC_rna$ENSG, 
          mart = ensembl)

EC_rna<-merge(EC_rna,ID,by.x="ENSG",by.y="ensembl_gene_id")%>%
  filter(gene_biotype=="protein_coding")%>%filter(hgnc_symbol!='')

writexl::write_xlsx(EC_rna,"RNA_EC.xlsx")
writexl::write_xlsx(EC_ribo,"ribo_EC.xlsx")




EC_rna%>%ggplot()+geom_point(aes(x=log2FC_E,y=-log10(padj)),
                             color=case_when(EC_rna$log2FC_E>=1 &EC_rna$padj<=0.05~"blue",
                                             EC_rna$log2FC_E<= -1 & EC_rna$padj <= 0.05~"red",
                                             T~"black"))+
  theme_classic()


EC_rna%>%filter(log2FC_E>=1 & padj <= 0.05)%>%nrow()



##############################
QC_rna%>%ggplot()+geom_point(aes(x=log2FC_RNA,y=-log10(padj)),
                             color=case_when(QC_rna$log2FC_RNA>=1 &QC_rna$padj<=0.05~"blue",
                                             QC_rna$log2FC_RNA<= -1 & QC_rna$padj <= 0.05~"red",
                                             T~"black"))+
  theme_classic()



names(EC_rna)[names(EC_rna)=="log2FoldChange"]<-"log2FC_RNA"
names(EC_ribo)[names(EC_ribo)=="log2FoldChange"]<-"log2FC_Ribo"


EC_RiboRNA<-merge(EC_ribo,EC_rna,by="ENSG")%>%na.omit()

EC_RiboRNA%>%ggplot(aes(y=log2FC_Ribo,x=log2FC_RNA))+geom_point(shape=21)+
  theme_classic()+stat_smooth(method="lm")+
  xlab("RNA-Seq Log2FC (Etop/Cyc)")+
  ylab("Ribo-Seq Log2FC (Etop/Cyc)")


write_xlsx(EC_RiboRNA,"EC_RiboRNA.xlsx")

fit_EC<-lm(EC_RiboRNA$log2FC_Ribo~EC_RiboRNA$log2FC_RNA)
summary(fit_EC)



Initation_Factors<-EC_RiboRNA$hgnc_symbol.x[str_detect(EC_RiboRNA$hgnc_symbol.x,"EIF[0-9]")]
Elongation_Factors<-EC_RiboRNA$hgnc_symbol.x[str_detect(EC_RiboRNA$hgnc_symbol.x,"EEF[0-9]")]
Ribosomal_Proteins<-EC_RiboRNA$hgnc_symbol.x[str_detect(EC_RiboRNA$hgnc_symbol.x,"^RPS[0-9]|^RPL[0-9]")]

data.frame(colnames(EC_RiboRNA))

Trans_Seq<-EC_RiboRNA%>%select(1,8,3,7,12,16,18)
Trans_Seq<-Trans_Seq%>%mutate(Type=case_when(Trans_Seq$hgnc_symbol.x%in%Initation_Factors~"Initiation",
                                             Trans_Seq$hgnc_symbol.x%in%Elongation_Factors~"Elongation",
                                             Trans_Seq$hgnc_symbol.x%in%Ribosomal_Proteins~"Ribosomal",
                                             T~"Other"))

Trans_E<-Trans_Seq%>%filter(Type=="Ribosomal")



Trans_E%>%ggplot()+geom_histogram(aes(x=log2FC_RNA),color=
                                    "white")

############################# QC






IR_RiboRNA%>%ggplot()+geom_point(aes(x=log2FoldChange.x,y=log2FoldChange.y),shape=21)+
  coord_cartesian(xlim=c(-8,8),ylim=c(-8,8))+
  theme_classic()+
  xlab("RNA-Seq Log2FC (IR/Cyc")+
  ylab("Ribo-Seq Log2FC (IR/Cyc)")


################################################################ QC

QC_ribo <- as.data.frame(results(ddsMat_ribo, contrast=list("Condition_Quiescent_vs_Cycling")))
QC_rna <- as.data.frame(results(ddsMat_rna, contrast=list("Condition_Quiescent_vs_Cycling")))

QC_ribo$ENSG<-rownames(QC_ribo)
QC_rna$ENSG<-rownames(QC_rna)


ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = QC_ribo$ENSG, 
          mart = ensembl)

QC_ribo<-merge(QC_ribo,ID,by.x="ENSG",by.y="ensembl_gene_id")%>%
  filter(gene_biotype=="protein_coding")%>%filter(hgnc_symbol!='')

ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = QC_rna$ENSG, 
          mart = ensembl)

QC_rna<-merge(QC_rna,ID,by.x="ENSG",by.y="ensembl_gene_id")%>%
  filter(gene_biotype=="protein_coding")%>%filter(hgnc_symbol!='')
writexl::write_xlsx(QC_rna,"RNA_QC.xlsx")
writexl::write_xlsx(QC_ribo,"ribo_QC.xlsx")

names(QC_ribo)[names(QC_ribo)=="log2FoldChange"]<-"log2FC_Ribo"
names(QC_rna)[names(QC_rna)=="log2FoldChange"]<-"log2FC_RNA"


QC_RiboRNA<-merge(QC_ribo,QC_rna,by="ENSG")%>%na.omit()


QC_RiboRNA%>%ggplot(aes(y=log2FC_Ribo,x=log2FC_RNA))+geom_point(shape=21)+
  coord_cartesian(xlim=c(-8,8),ylim=c(-8,8))+
  stat_smooth(method="lm")+
  theme_classic()+
  xlab("RNA-Seq Log2FC (Qui/Cyc")+
  ylab("Ribo-Seq Log2FC (Qui/Cyc)")



fit_QC<-lm(QC_RiboRNA$log2FC_Ribo~QC_RiboRNA$log2FC_RNA)
summary(fit_QC)

################################################################ IR


IR_ribo <- as.data.frame(results(ddsMat_ribo, contrast=list("Condition_Ionizing_vs_Cycling")))
IR_rna <- as.data.frame(results(ddsMat_rna, contrast=list("Condition_Ionizing_vs_Cycling")))

IR_ribo$ENSG<-rownames(IR_ribo)
IR_rna$ENSG<-rownames(IR_rna)


ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = IR_rna$ENSG, 
          mart = ensembl)


IR_ribo<-merge(IR_ribo,ID,by.x="ENSG",by.y="ensembl_gene_id")%>%
  filter(gene_biotype=="protein_coding")%>%filter(hgnc_symbol!='')

IR_rna<-merge(IR_rna,ID,by.x="ENSG",by.y="ensembl_gene_id")%>%
  filter(gene_biotype=="protein_coding")%>%filter(hgnc_symbol!='')
writexl::write_xlsx(IR_rna,"RNA_IR.xlsx")
writexl::write_xlsx(IR_ribo,"ribo_IR.xlsx")




IR_RiboRNA_Sig<-merge(IR_ribo,IR_rna,by="ENSG")%>%
  filter(padj.x<=0.05&padj.y<=0.05)

plot(IR_RiboRNA_Sig$log2FoldChange.x,IR_RiboRNA_Sig$log2FoldChange.y)

IR_RiboRNA<-merge(IR_ribo,IR_rna,by="ENSG")%>%na.omit()
plot(IR_RiboRNA$log2FoldChange.x,IR_RiboRNA$log2FoldChange.y)

IR_RiboRNA%>%ggplot()+geom_point(aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  coord_cartesian(xlim=c(-8,8),ylim=c(-8,8))




############################## Correlation dRNA RNA-seq

Seq_Comp<-merge(EC_rna,E_C_dRNA,by.x="hgnc_symbol", by.y="genename")%>%
  na.omit()


plot(Seq_Comp$log2FC_RNA,Seq_Comp$log2FoldChange)
cor(Seq_Comp$log2FC_RNA,Seq_Comp$log2FoldChange)

summary(lm(Seq_Comp$log2FoldChange~Seq_Comp$log2FC_RNA))


Seq_Comp_sub<-Seq_Comp%>%filter(log2FC_RNA<4&log2FoldChange<4)
Seq_Comp%>%ggplot(aes(x=log2FC_RNA,y=log2FoldChange))+geom_point(shape=21)+
  stat_smooth(method="lm")+theme_classic()


cor(y=Seq_Comp_sub$log2FoldChange,x=Seq_Comp_sub$log2FC_RNA,method="spearman")


