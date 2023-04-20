library(reshape2)
library(tidyverse)
library(readxl)

##################### I call it TPM but it's actually RPKM
##################################### Gene lengths
library(tidyverse)
library(biomaRt)
values<-RNA_Counts$ENSG
x<-listAttributes(ensembl)
ID<-getBM(attributes = c('ensembl_gene_id','transcript_length',"gene_biotype"),
          filters = 'ensembl_gene_id',
          values = values, 
          mart = ensembl)

ID$ensembl_gene_id<-as.factor(ID$ensembl_gene_id)

ID<-ID%>%filter(gene_biotype=="protein_coding")

Gene_Lengths<-ID%>%group_by(ensembl_gene_id)%>%summarize(ENSG=ensembl_gene_id,Length=max(transcript_length))%>%
  unique()

Gene_Lengths$Length<-Gene_Lengths$Length/1000


##################################### RPKM files

####Reads are generated from WC unix on fastq files
Seq_Totals<-read_xlsx("Read_Counts_Lane1.xlsx")


####divide by 1e6; row 1 is ribo row 2 is rna 
Seq_PM<-Seq_Totals[,2:13]/1e6



RNA_Counts_TPM<-RNA_Counts

RNA_Counts_TPM$C1_R<-(RNA_Counts_TPM$C1/Seq_Totals$C1[2])
RNA_Counts_TPM$C2_R<-(RNA_Counts_TPM$C2/Seq_Totals$C2[2])
RNA_Counts_TPM$Q1_R<-(RNA_Counts_TPM$Q1/Seq_Totals$Q1[2])
RNA_Counts_TPM$Q2_R<-(RNA_Counts_TPM$Q2/Seq_Totals$Q2[2])
RNA_Counts_TPM$Q3_R<-(RNA_Counts_TPM$Q3/Seq_Totals$Q3[2])
RNA_Counts_TPM$E1_R<-(RNA_Counts_TPM$E1/Seq_Totals$E1[2])
RNA_Counts_TPM$E2_R<-(RNA_Counts_TPM$E2/Seq_Totals$E2[2])
RNA_Counts_TPM$E3_R<-(RNA_Counts_TPM$E3/Seq_Totals$E3[2])
RNA_Counts_TPM$IR1_R<-(RNA_Counts_TPM$IR1/Seq_Totals$IR1[2])
RNA_Counts_TPM$IR2_R<-(RNA_Counts_TPM$IR2/Seq_Totals$IR2[2])
RNA_Counts_TPM$IR3_R<-(RNA_Counts_TPM$IR3/Seq_Totals$IR3[2])


Ribo_Counts_TPM<-Ribo_Counts
Ribo_Counts_TPM$C1_r<-(Ribo_Counts_TPM$C1/Seq_Totals$C1[1])
Ribo_Counts_TPM$C2_r<-(Ribo_Counts_TPM$C2/Seq_Totals$C2[1])
Ribo_Counts_TPM$Q1_r<-(Ribo_Counts_TPM$Q1/Seq_Totals$Q1[1])
Ribo_Counts_TPM$Q2_r<-(Ribo_Counts_TPM$Q2/Seq_Totals$Q2[1])
Ribo_Counts_TPM$Q3_r<-(Ribo_Counts_TPM$Q3/Seq_Totals$Q3[1])
Ribo_Counts_TPM$E1_r<-(Ribo_Counts_TPM$E1/Seq_Totals$E1[1])
Ribo_Counts_TPM$E2_r<-(Ribo_Counts_TPM$E2/Seq_Totals$E2[1])
Ribo_Counts_TPM$E3_r<-(Ribo_Counts_TPM$E3/Seq_Totals$E3[1])
Ribo_Counts_TPM$IR1_r<-(Ribo_Counts_TPM$IR1/Seq_Totals$IR1[1])
Ribo_Counts_TPM$IR2_r<-(Ribo_Counts_TPM$IR2/Seq_Totals$IR2[1])
Ribo_Counts_TPM$IR3_r<-(Ribo_Counts_TPM$IR3/Seq_Totals$IR3[1])


All_Counts_TPM<-merge(Ribo_Counts_TPM,RNA_Counts_TPM,by="ENSG")%>%
  select(1,13:23,35:45)

####Final convert to RPKM by dividing RPM reads over gen length in Kb


All_Counts_TPM<-merge(All_Counts_TPM,Gene_Lengths,by="ENSG")

All_Counts_RPM<-All_Counts_TPM%>%mutate(C_RNA=mapply(col_avg_Cyc,C1_R,C2_R),
                                        Q_RNA=mapply(col_avg,Q1_R,Q2_R,Q3_R),
                                        E_RNA=mapply(col_avg,E1_R,E2_R,E3_R),
                                        IR_RNA=mapply(col_avg,IR1_R,IR2_R,IR3_R),
                                        C_ribo=mapply(col_avg_Cyc,C1_r,C2_r),
                                        Q_ribo=mapply(col_avg,Q1_r,Q2_r,Q3_r),
                                        E_ribo=mapply(col_avg,E1_r,E2_r,E3_r),
                                        IR_ribo=mapply(col_avg,IR1_r,IR2_r,IR3_r))

ENSG00000276168

data.frame(colnames(All_Counts_TPM))

All_Counts_TPM[,2:23]<-All_Counts_TPM[,2:23]/All_Counts_TPM$Length

############################## Averages and Figures

col_avg<-function(A,B,C){
  x<-c(A,B,C)
  mean(x)
}

col_avg_Cyc<-function(A,B){
  x<-c(A,B)
  mean(x)
}

All_Counts_RPKM<-All_Counts_TPM%>%mutate(C_RNA=mapply(col_avg_Cyc,C1_R,C2_R),
                        Q_RNA=mapply(col_avg,Q1_R,Q2_R,Q3_R),
                        E_RNA=mapply(col_avg,E1_R,E2_R,E3_R),
                        IR_RNA=mapply(col_avg,IR1_R,IR2_R,IR3_R),
                        C_ribo=mapply(col_avg_Cyc,C1_r,C2_r),
                        Q_ribo=mapply(col_avg,Q1_r,Q2_r,Q3_r),
                        E_ribo=mapply(col_avg,E1_r,E2_r,E3_r),
                        IR_ribo=mapply(col_avg,IR1_r,IR2_r,IR3_r))
summary(All_Counts_TPM)
data.frame(colnames(All_Counts_RPKM))
Avg_RPKM<-All_Counts_RPKM%>%select(24:33)

Avg_RPKM[Avg_RPKM==0]<-NA
  
Avg_RPKM<-Avg_RPKM%>%na.omit()


NoLo<-Avg_RPKM

NoLo<-NoLo%>%mutate(C_RNA=ifelse(C_RNA>quantile(C_RNA)[2],C_RNA,NA),
                    Q_RNA=ifelse(Q_RNA>quantile(Q_RNA)[2],Q_RNA,NA),
                    E_RNA=ifelse(E_RNA>quantile(E_RNA)[2],E_RNA,NA),
                    IR_RNA=ifelse(IR_RNA>quantile(IR_RNA)[2],IR_RNA,NA),
                    C_ribo=ifelse(C_ribo>quantile(C_ribo)[2],C_ribo,NA),
                    Q_ribo=ifelse(Q_ribo>quantile(Q_ribo)[2],Q_ribo,NA),
                    E_ribo=ifelse(E_ribo>quantile(E_ribo)[2],E_ribo,NA),
                    IR_ribo=ifelse(IR_ribo>quantile(IR_ribo)[2],IR_ribo,NA)
                    )%>%na.omit()



summary(Avg_TPM)

Avg_RPKM<-Avg_RPKM%>%mutate(C_All=(C_ribo/C_RNA),
                          E_All=(E_ribo/E_RNA),
                          Q_All=(Q_ribo/Q_RNA),
                          IR_All=(IR_ribo/IR_RNA))
Avg_RPKM<-Avg_RPKM%>%na.omit()

Box_RPKM<-Avg_RPKM%>%select(11:14)

Box_TPM_Wide<-pivot_longer(Box_TPM,cols = 1:4,names_to = "Sample",values_to = "Ribo_RNA")
Box_TPM_Wide$Sample<-factor(c(Box_TPM_Wide$Sample),levels = c("C_All","Q_All","E_All","IR_All"))

Box_TPM_Wide%>%ggplot(aes(x=Sample,y=log10(Ribo_RNA)))+geom_boxplot(aes(fill=Sample))+theme_classic()


NoLo%>%ggplot()+stat_smooth(aes(y=E_ribo,x=E_RNA),color="red",method = "lm")+
  stat_smooth(aes(y=C_ribo,x=C_RNA),color="black",method = "lm")+
  stat_smooth(aes(y=IR_ribo,x=IR_RNA),color="dark green",method = "lm")+
  stat_smooth(aes(y=Q_ribo,x=Q_RNA),color="orange",method = "lm")+
  xlab("RNA-Seq (RPKM)")+ylab("Ribo-Seq (RPKM)")+
  theme_classic()

Avg_RPKM%>%ggplot()+stat_smooth(aes(y=E_ribo,x=E_RNA),color="red",method = "lm")+
  stat_smooth(aes(y=C_ribo,x=C_RNA),color="black",method = "lm")+
  stat_smooth(aes(y=IR_ribo,x=IR_RNA),color="dark green",method = "lm")+
  stat_smooth(aes(y=Q_ribo,x=Q_RNA),color="orange",method = "lm")+
  xlab("RNA-Seq (RPKM)")+ylab("Ribo-Seq (RPKM)")+
  theme_classic()


NoLo%>%ggplot()+geom_point(aes(y=E_ribo,x=E_RNA))

  scale_x_log10()+scale_y_log10()

Avg_RPKM%>%ggplot()+geom_point(aes(y=C_ribo,x=C_RNA))+
  scale_x_log10()+scale_y_log10()

Avg_RPKM%>%ggplot()+geom_point(aes(y=Q_ribo,x=Q_RNA))+
  scale_x_log10()+scale_y_log10()


plot(Avg_TPM$C_RNA,Avg_TPM$C_ribo)


fit_c<-summary(lm(Avg_RPKM$C_ribo~Avg_RPKM$C_RNA))
fit_c

fit_q<-summary(lm(Avg_RPKM$Q_ribo~Avg_RPKM$Q_RNA))
fit_q


fit_e<-summary(lm(Avg_RPKM$E_ribo~Avg_RPKM$E_RNA))
fit_e


fit_I<-summary(lm(Avg_RPKM$IR_ribo~Avg_RPKM$IR_RNA))
fit_I

################################


Avg_TPM<-Avg_TPM%>%mutate(EC_RNA=log2(E_RNA/C_RNA),
                          EC_Ribo=log2(E_ribo/C_ribo),
                          IC_RNA=log2(IR_RNA/C_RNA),
                          IC_Ribo=log2(IR_ribo/C_ribo),
                          QC_RNA=log2(Q_RNA/C_RNA),
                          QC_Ribo=log2(Q_ribo/C_ribo))

plot(Avg_TPM$EC_RNA,Avg_TPM$EC_Ribo)
plot(Avg_TPM$QC_RNA,Avg_TPM$QC_Ribo)
plot(Avg_TPM$IC_RNA,Avg_TPM$IC_Ribo)



Avg_TPM%>%ggplot()+geom_point()

#########################################Correlelogram
install.packages("ggcorrplot")
library(ggcorrplot)
data.frame(colnames(RNA_Counts_TPM))




gene_length<-read_tsv("mart_export.txt")%>%na.omit()
gene_length<-gene_length%>%filter(gene_length$`Gene stable ID`%in%RNA_Counts_TPM$ENSG)%>%
  select(1,3,6)


names(gene_length)<-c("ENSG","ENST","Length")

class(gene_length$ENSG)
gene_length$ENSG<-as.factor(gene_length$ENSG)


summary<-gene_length%>%group_by(ENSG)%>%
  summarise(ENSG=ENSG,
            Length=max(Length))%>%
  unique()%>%as.data.frame()

class(summary)




TPM_RNA<-RNA_Counts_TPM

TPM_RNA<-TPM_RNA[order(TPM_RNA$ENSG),]
summary<-summary[order(summary$ENSG),]

merge_TPM<-merge(TPM_RNA,summary,by="ENSG")%>%na.omit()

plot(log10(merge_TPM$IR3_R+1),log10(merge_TPM$Length+1))

TPM_RNA<-TPM_RNA[,2:12]
TPM_Ribo<-Ribo_Counts_TPM[,13:23]

RNA_cor<-cor(TPM_RNA)

Ribo_cor<-cor(TPM_Ribo)

ggcorrplot(RNA_cor)+
  scale_fill_gradient(limits=c(min(as.vector(RNA_cor)),max(as.vector(RNA_cor))),
                      low="red",high="steelblue")

ggcorrplot(Ribo_cor)+
  scale_fill_gradient(limits=c(min(as.vector(Ribo_cor)),max(as.vector(Ribo_cor))),
                      low="red",high="steelblue")



################################################## Rank Correlations



Avg_RPKM_filtered<-Avg_RPKM[rowSums(Avg_RPKM==0)==0,]

Avg_RPKM_filtered$res_E= lm(log10(E_ribo)~log10(E_RNA),data=Avg_RPKM_filtered)$residuals
summary(lm(log10(Avg_RPKM_filtered$E_ribo)~log10(Avg_RPKM_filtered$E_RNA)))


res_2E<-2*sd(Avg_RPKM_filtered$res_E)

Trans_Etop<-Avg_RPKM_filtered%>%filter(res_E>=res_2E)
Repress_Etop<-Avg_RPKM_filtered%>%filter(res_E<=-res_2E)

Avg_RPKM_filtered%>%ggplot()+geom_point(aes(x=log10(E_RNA),y=log10(E_ribo)),
                               color=ifelse(Avg_RPKM_filtered$ensembl_gene_id%in%Trans_Etop$ensembl_gene_id,"blue","black"))



values<-Avg_RPKM_filtered$ensembl_gene_id

ID<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
          filters = 'ensembl_gene_id',
          values = values, 
          mart = ensembl)

Trans_Etop<-merge(Trans_Etop,ID,by="ensembl_gene_id")

Repress_Etop<-merge(Repress_Etop,ID,by="ensembl_gene_id")
writexl::write_xlsx(Repress_Etop,"Repress_Etop.xlsx")
writexl::write_xlsx(Trans_Etop,"Trans_Etop.xlsx")
