library(tidyverse)
library(writexl)

CS_B1<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample1_SEN2_batch1-sample2_CYC2_batch1/sample1_SEN2_batch1-sample2_CYC2_batch1-single.lib_diff_iso.txt"
CS_B2<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample3_SEN3_batch2-sample4_CYC3_batch2/sample3_SEN3_batch2-sample4_CYC3_batch2-single.lib_diff_iso.txt"

ETgE_B1<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample1_SEN2_batch1-sample5_ETGA_1_batch1/sample1_SEN2_batch1-sample5_ETGA_1_batch1-single.lib_diff_iso.txt"
ETgE_B2<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample3_SEN3_batch2-sample6_ETGB_1_batch2/sample3_SEN3_batch2-sample6_ETGB_1_batch2-single.lib_diff_iso.txt"

ESvE_B1<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample1_SEN2_batch1-sample9_ETSV_1_batch1/sample1_SEN2_batch1-sample9_ETSV_1_batch1-single.lib_diff_iso.txt"
ESvE_B2<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample3_SEN3_batch2-sample10_ETSV_2_batch2/sample3_SEN3_batch2-sample10_ETSV_2_batch2-single.lib_diff_iso.txt"

CTgC_B1<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample2_CYC2_batch1-sample11_CTGA_2_batch1/sample2_CYC2_batch1-sample11_CTGA_2_batch1-single.lib_diff_iso.txt"
CTgC_B2<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample4_CYC3_batch2-sample12_CTGB_2_batch2/sample4_CYC3_batch2-sample12_CTGB_2_batch2-single.lib_diff_iso.txt"

CSvC_B1<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample2_CYC2_batch1-sample7_CSVA_1_batch1/sample2_CYC2_batch1-sample7_CSVA_1_batch1-single.lib_diff_iso.txt"
CSvC_B2<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample4_CYC3_batch2-sample8_CSVA_2_batch2/sample4_CYC3_batch2-sample8_CSVA_2_batch2-single.lib_diff_iso.txt"
  
  
ETgC_B1<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample2_CYC2_batch1-sample5_ETGA_1_batch1/sample2_CYC2_batch1-sample5_ETGA_1_batch1-single.lib_diff_iso.txt"
ETgC_B2<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample4_CYC3_batch2-sample6_ETGB_1_batch2/sample4_CYC3_batch2-sample6_ETGB_1_batch2-single.lib_diff_iso.txt"

ESvC_B1<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample2_CYC2_batch1-sample9_ETSV_1_batch1/sample2_CYC2_batch1-sample9_ETSV_1_batch1-single.lib_diff_iso.txt"
ESvC_B2<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/isoform_FLAIR/results/concated_all_samples/diff_iso_usage/sample4_CYC3_batch2-sample10_ETSV_2_batch2/sample4_CYC3_batch2-sample10_ETSV_2_batch2-single.lib_diff_iso.txt"


Cyc_Sen_Batch1<-CS_B1%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","Sen_Count","Alt_Cyc_Count","Alt_Sen_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%
  mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(Sen_Count/(Sen_Count+Alt_Sen_Count)),
         Ratio_C=Cyc_Count/(Cyc_Count+Alt_Cyc_Count),Ratio_S=(Sen_Count/(Sen_Count+Alt_Sen_Count)))

Cyc_Sen_Batch2<-CS_B2%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","Sen_Count","Alt_Cyc_Count","Alt_Sen_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%
  mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(Sen_Count/(Sen_Count+Alt_Sen_Count)),
         Ratio_C=Cyc_Count/(Cyc_Count+Alt_Cyc_Count),Ratio_S=(Sen_Count/(Sen_Count+Alt_Sen_Count)))





Cyc_Sen_Sig<-merge(Cyc_Sen_Batch1,Cyc_Sen_Batch2,by="Isoform")%>%filter(pValue.x<=0.05&pValue.y<=0.05)%>%
  filter(sign(Usage_Diff.x)==sign(Usage_Diff.y))%>%
  dplyr::select("Isoform","ENSG.x","pValue.x","pValue.y","Usage_Diff.x","Usage_Diff.y")

FC<-function(A,B){
  x<-c(A,B)
  mean(x)
}

Cyc_Sen_Sig<-Cyc_Sen_Sig%>%mutate(FC=mapply(FC,Usage_Diff.x,Usage_Diff.y))

Cyc_Sen_Sig%>%ggplot()+geom_point(aes(x=Usage_Diff.x,y=-log10(pValue.x)))

writexl::write_xlsx(Cyc_Sen_Sig,"Cyc_Sen_Sig.xlsx")



Cyc_Sen<-merge(Cyc_Sen_Batch1,Cyc_Sen_Batch2,by="Isoform")%>%pull(ENSG.x)%>%unique()

Cyc_Sen_Sig<-merge(Cyc_Sen_Sig,ID,by.x = "ENSG.x",by.y="ensembl_gene_id")
write_xlsx(Cyc_Sen_Sig,"Table_S3_CycSen.xlsx")


Cyc_Sen_Sig$ENSG.x%>%unique()%>%length()
Cyc_Sen_Sig$Isoform%>%unique()%>%str_detect("ENST")%>%sum()
EC_Ribo<-Cyc_Sen_Sig%>%filter(ENSG.x%in%Ribo_RNA$ENSG)
EC_Ribo$ENSG.x%>%unique()%>%length()



###################

ETg_Sen_Batch1<-ETgE_B1%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Sen_Count","ETg_Count","Alt_Sen_Count","Alt_ETg_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Sen_Count/(Sen_Count+Alt_Sen_Count))-(ETg_Count/(ETg_Count+Alt_ETg_Count)))

ETg_Sen_Batch2<-ETgE_B2%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Sen_Count","ETg_Count","Alt_Sen_Count","Alt_ETg_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Sen_Count/(Sen_Count+Alt_Sen_Count))-(ETg_Count/(ETg_Count+Alt_ETg_Count)))


ETg_Sen_Sig<-merge(ETg_Sen_Batch1,ETg_Sen_Batch2,by="Isoform")%>%filter(pValue.x<=0.05&pValue.y<=0.05)%>%
  filter(sign(Usage_Diff.x)==sign(Usage_Diff.y))%>%dplyr::select("Isoform","ENSG.x","pValue.x","pValue.y","Usage_Diff.x","Usage_Diff.y")

ETg_Sen_Sig<-merge(ETg_Sen_Sig,ID,by.x = "ENSG.x",by.y="ensembl_gene_id")
write_xlsx(ETg_Sen_Sig,"Table_S3_ETg.xlsx")

ETg_Sen_Sig$ENSG.x%>%unique()%>%length()
ETg_Sen_Sig$Isoform%>%unique()%>%str_detect("ENST")%>%sum()

#########################

ESv_Sen_Batch1<-ESvE_B1%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Sen_Count","ESv_Count","Alt_Sen_Count","Alt_ESv_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Sen_Count/(Sen_Count+Alt_Sen_Count))-(ESv_Count/(ESv_Count+Alt_ESv_Count)))

ESv_Sen_Batch2<-ESvE_B2%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Sen_Count","ESv_Count","Alt_Sen_Count","Alt_ESv_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Sen_Count/(Sen_Count+Alt_Sen_Count))-(ESv_Count/(ESv_Count+Alt_ESv_Count)))


ESv_Sen_Sig<-merge(ESv_Sen_Batch1,ESv_Sen_Batch2,by="Isoform")%>%filter(pValue.x<=0.05&pValue.y<=0.05)%>%
  filter(sign(Usage_Diff.x)==sign(Usage_Diff.y))%>%dplyr::select("Isoform","ENSG.x","pValue.x","pValue.y","Usage_Diff.x","Usage_Diff.y")

ESv_Sen_Sig<-merge(ESv_Sen_Sig,ID,by.x = "ENSG.x",by.y="ensembl_gene_id")

write_xlsx(ESv_Sen_Sig,"Table_S3_ESv.xlsx")


ESv_Sen_Sig$ENSG.x%>%unique()%>%length()
ESv_Sen_Sig$Isoform%>%unique()%>%str_detect("ENST")%>%sum()

###############################

CTg_Cyc_Batch1<-CTgC_B1%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","CTg_Count","Alt_Cyc_Count","Alt_CTg_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(CTg_Count/(CTg_Count+Alt_CTg_Count)))

CTg_Cyc_Batch2<-CTgC_B2%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","CTg_Count","Alt_Cyc_Count","Alt_CTg_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(CTg_Count/(CTg_Count+Alt_CTg_Count)))

CTg_Cyc_Sig<-merge(CTg_Cyc_Batch1,CTg_Cyc_Batch2,by="Isoform")%>%filter(pValue.x<=0.05&pValue.y<=0.05)%>%
  dplyr::select("Isoform","ENSG.x","pValue.x","Usage_Diff.x","Usage_Diff.y")%>%filter(sign(Usage_Diff.x)==sign(Usage_Diff.y))

CTg_Cyc_Sig<-merge(CTg_Cyc_Sig,ID,by.x = "ENSG.x",by.y="ensembl_gene_id")
write_xlsx(CTg_Cyc_Sig,"Table_S3_CTg.xlsx")



CTg_Cyc_Sig$ENSG.x%>%unique()%>%length()
CTg_Cyc_Sig$Isoform%>%unique()%>%str_detect("ENST")%>%sum()


#########################

CSv_Cyc_Batch1<-CSvC_B1%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","CSv_Count","Alt_Cyc_Count","Alt_CSv_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(CSv_Count/(CSv_Count+Alt_CSv_Count)))

CSv_Cyc_Batch2<-CSvC_B2%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","CSv_Count","Alt_Cyc_Count","Alt_CSv_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(CSv_Count/(CSv_Count+Alt_CSv_Count)))

CSv_Cyc_Sig<-merge(CSv_Cyc_Batch1,CSv_Cyc_Batch2,by="Isoform")%>%filter(pValue.x<=0.05&pValue.y<=0.05)%>%
filter(sign(Usage_Diff.x)==sign(Usage_Diff.y))%>%
  dplyr::select("Isoform","ENSG.x","pValue.x","pValue.y","Usage_Diff.x","Usage_Diff.y")


CSv_Cyc_Sig<-merge(CSv_Cyc_Sig,ID,by.x = "ENSG.x",by.y="ensembl_gene_id")
write_xlsx(CSv_Cyc_Sig,"Table_S3_CSv.xlsx")

CSv_Cyc_Sig$ENSG.x%>%unique()%>%length()
CSv_Cyc_Sig$Isoform%>%unique()%>%str_detect("ENST")%>%sum()


##########################

ESv_Cyc_Batch1<-ESvC_B1%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","ESv_Count","Alt_Cyc_Count","Alt_ESv_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(ESv_Count/(ESv_Count+Alt_ESv_Count)))

ESv_Cyc_Batch2<-ESvC_B2%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","ESv_Count","Alt_Cyc_Count","Alt_ESv_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(ESv_Count/(ESv_Count+Alt_ESv_Count)))

ESv_Cyc_Sig<-merge(ESv_Cyc_Batch1,ESv_Cyc_Batch2,by="Isoform")%>%filter(pValue.x<=0.05&pValue.y<=0.05)%>%
dplyr::select("Isoform","ENSG.x","pValue.x","Usage_Diff.x","Usage_Diff.y")
ESv_Cyc_Sig<-merge(ESv_Cyc_Sig,ID,by.x = "ENSG.x",by.y="ensembl_gene_id")
ESv_Cyc_Sig$ENSG.x%>%unique()%>%length()
ESv_Cyc_Sig$Isoform%>%unique()%>%str_detect("ENST")%>%sum()



ETg_Cyc_Batch1<-ETgC_B1%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","ETg_Count","Alt_Cyc_Count","Alt_ETg_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(ETg_Count/(ETg_Count+Alt_ETg_Count)))

ETg_Cyc_Batch2<-ETgC_B2%>%
  map_dfr(read_tsv,col_names=c("ENSG","Isoform","pValue","Cyc_Count","ETg_Count","Alt_Cyc_Count","Alt_ETg_Count"))%>%
  filter(str_detect(ENSG,"ENSG"))%>%mutate(Usage_Diff=(Cyc_Count/(Cyc_Count+Alt_Cyc_Count))-(ETg_Count/(ETg_Count+Alt_ETg_Count)))

ETg_Cyc_Sig<-merge(ETg_Cyc_Batch1,ETg_Cyc_Batch2,by="Isoform")%>%filter(pValue.x<=0.05&pValue.y<=0.05)%>%
  dplyr::select("Isoform","ENSG.x","pValue.x","Usage_Diff.x","Usage_Diff.y")%>%filter(sign(Usage_Diff.x)==sign(Usage_Diff.y))


ETg_Cyc_Sig$ENSG.x%>%unique()%>%length()
ETg_Cyc_Sig$Isoform%>%unique()%>%str_detect("ENST")%>%sum()



#######################################################################

library(readxl)
library(writexl)
write_xlsx(ETg_Cyc_Sig,"ETg_Cyc_Sig.xlsx")


write_xlsx(ESv_Cyc_Sig,"ESv_Cyc_Sig.xlsx")

write_xlsx(Cyc_Sen_Sig,"Cyc_Sen_Sig.xlsx")

write_xlsx(EC_Ribo,"EC_Ribo.xlsx")


RNA_All<-read_xlsx("RNA_All.xlsx")

EC<-RNA_All%>%filter(ENSG%in%Cyc_Sen_Sig$ENSG.x)
EC_Ribo<-EC%>%filter(ENSG%in%Ribo_RNA$ENSG)
write_xlsx(EC,"EC.xlsx")

CSvC<-RNA_All%>%filter(ENSG%in%CSv_Cyc_Sig$ENSG.x)
write_xlsx(CSvC,"CSvC.xlsx")

ESv<-RNA_All%>%filter(ENSG%in%ESv_Sen_Sig$ENSG.x)
write_xlsx(ESv,"ESv.xlsx")
ETg<-RNA_All%>%filter(ENSG%in%ETg_Sen_Sig$ENSG.x)
write_xlsx(ETg,"ETg.xlsx")

CTgC<-RNA_All%>%filter(ENSG%in%CTg_Cyc_Sig$ENSG.x)
write_xlsx(CTgC,"CTgC.xlsx")

CTgC_Ribo<-CTgC%>%filter(ENSG%in%Ribo_RNA$ENSG)
ETgE_Ribo<-ETg%>%filter(ENSG%in%Ribo_RNA$ENSG)






ESvC<-RNA_All%>%filter(ENSG%in%ESv_Cyc_Sig$ENSG.x)
write_xlsx(ESvC,"ESvC.xlsx")
ETgC<-RNA_All%>%filter(ENSG%in%ETg_Cyc_Sig$ENSG.x)
write_xlsx(ETgC,"ETgC.xlsx")

Common<-EC%>%filter(ENSG%in%ESv$ENSG)%>%filter(ENSG%in%ETg$ENSG)
Common_Tg<-merge(CTg_Cyc_Sig,ETg_Sen_Sig,by="Isoform")%>%filter(sign(Usage_Diff.x.x)==sign(Usage_Diff.x.y))
Exclusive_Tg_Etop<-ETg_Sen_Sig%>%filter(!(ETg_Sen_Sig$Isoform%in%Common_Tg$Isoform))


Common_Tg$ENSG.x.x%>%unique()%>%length()
Exclusive_Tg_Cyc$ENSG.x%>%unique()%>%length()
Exclusive_Tg_Etop$ENSG.x%>%unique()%>%length()

Common_ESv<-merge(CSv_Cyc_Sig,ESv_Sen_Sig,by="Isoform")%>%filter(sign(Usage_Diff.x.x)==sign(Usage_Diff.x.y))


write_xlsx(Common_Tg,"Common_Tg.xlsx")

write_xlsx(Common_ESv,"Common_Sv.xlsx")

write_xlsx(Exclusive_Tg_Etop,"Excusive_Tg_Etop.xlsx")

Exclusive_Tg_Cyc<-CTg_Cyc_Sig%>%filter(!(CTg_Cyc_Sig$Isoform%in%Common_Tg$Isoform))

Exclusive_Tg_Cyc%>%filter()


write_xlsx(Exclusive_Tg_Cyc,"Excusive_Tg_Cyc.xlsx")

Exclusive_Sv_Cyc<-CSv_Cyc_Sig%>%filter(!(CSv_Cyc_Sig$Isoform%in%Common_ESv$Isoform))
write_xlsx(Exclusive_Sv_Cyc,"Excusive_Sv_Cyc.xlsx")


Exclusive_ESv_Cyc<-ESv_Sen_Sig%>%filter(!(ESv_Sen_Sig$Isoform%in%Common_ESv$Isoform))

#####################################################################################
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

filters<-listFilters(ensembl)
attributes<-listAttributes(ensembl)

ENSG<-RNA_All$ENSG

ID<-getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = ENSG, 
          mart = ensembl)

Ribo_RNA<-RNA_All%>%filter(str_detect(hgnc_symbol.x,"^RPS|^RPL"))%>%filter(!str_detect(hgnc_symbol.x,"6K"))

Ex_ETg<-Exclusive_Tg_Etop$ENSG.x%>%unique()%>%.[.%in%Ribo_RNA$ENSG]
Ex_ETg<-Exclusive_Tg_Cyc%>%unique()%>%filter(.$ENSG.x%in%Ribo_RNA$ENSG)

Exclusive_Tg_Cyc$ENSG.x%>%unique()%>%.[.%in%Ribo_RNA$ENSG]%>%length()
Exclusive_Tg_Etop$ENSG.x%>%unique()%>%.[.%in%Ribo_RNA$ENSG]%>%length()

Exclusive_Sv_Cyc$ENSG.x%>%unique()%>%.[.%in%Ribo_RNA$ENSG]%>%length()
Exclusive_ESv_Cyc$ENSG.x%>%unique()%>%.[.%in%Ribo_RNA$ENSG]%>%length()

Common_Tg$ENSG.x.x%>%unique()%>%.[.%in%Ribo_RNA$ENSG]%>%length()

Common_ESv$ENSG.x.x%>%unique()%>%.[.%in%Ribo_RNA$ENSG]%>%length()


ETg_Cyc_Sig$ENSG.x%>%unique()%>%.[.%in%Ribo_RNA$ENSG]%>%length()


Common_Ribo<-ETg_Cyc_Sig%>%filter(ENSG.x%in%Cyc_Sen_Sig$ENSG.x & ENSG.x%in%ESv_Cyc_Sig$ENSG.x)%>%
  filter(ENSG.x%in%Ribo_RNA$ENSG)%>%merge(.,ID,by.x="ENSG.x",by.y="ensembl_gene_id")

Common_Ribo$ENSG.x%>%unique()%>%length()






############################################ Xbp1

Cyc_Sen_XBP1<-merge(Cyc_Sen_Batch1,Cyc_Sen_Batch2,by="Isoform")%>%filter(ENSG.x=="ENSG00000100219")
CTg_Cyc_XBP1<-merge(CTg_Cyc_Batch1,CTg_Cyc_Batch2,by="Isoform")%>%filter(ENSG.x=="ENSG00000100219")
ETg_Sen_XBP1<-merge(ETg_Sen_Batch1,ETg_Sen_Batch2,by="Isoform")%>%filter(ENSG.x=="ENSG00000100219")

write_xlsx(Cyc_Sen_XBP1,"EC_Xbp1.xlsx")

     