library(tidyverse)
cds_genic<-read_tsv("cds_genic.bed",col_names = F)
mRNA_genic<-read_tsv("mRNA.bed",col_names = F)

values<-cds_genic$X1

library(biomaRt)
ID<-getBM(attributes = c('hgnc_symbol','description','gene_biotype','ensembl_transcript_id'),
          filters = 'ensembl_transcript_id',
          values = values, 
          mart = ensembl)



cds_genic_labeled<-merge(cds_genic,ID,by.x="X1",by.y="ensembl_transcript_id")
mRNA_genic_labeled<-merge(mRNA_genic,ID,by.x="X1",by.y="ensembl_transcript_id")


RProtein_genc<-cds_genic_labeled%>%filter(hgnc_symbol%in%Ribosomal_Proteins)%>%select(1:6)

write_tsv(RProtein_genc,"RProtein_genic.bed")


CEBPB_genic<-cds_genic_labeled%>%filter(hgnc_symbol=="CEBPB")%>%select(1:6)
write_tsv(CEBPB_genic,"CEBPB_genic.bed")

SASP<-c("IL6","CXCL8","CCL2","CCL8","CCL7","CCL13","IGFBP3","CSF2","CXCL1",
        "CXCL2","CXCL3","ICAM1","GDF15","SERPINE1", "STC1",
        "TIMP1","TIMP2","SERPINEB2","IL1A","IL1B","IL13","IL15",
        "CXCL2","CXCL3","CCL8","CCL13","CCL3","CCL20","CCL16","CCL11",
        "CCL25","MIF","AREG","EREG","EGF","FGF2","HGF","FGF7","VEGFA",
        "ANG","KITLG","CXCL12","PGF","IGFBP2","IGFBP3","IGFBP4","IGFPB6","IGFPB7",
        "MMP1","MMP3","MMP10","MMP12","MMP13","MMP14","ICAM1","ICAM3",
        "TNFRSF11B","TNFRSF15A","TNFRSF10C","FAS","TNFRSF1A","PLAUR","IL6ST",
        "EGF")

SASP_genic<-cds_genic_labeled%>%filter(hgnc_symbol%in%SASP)%>%select(1:6)
write_tsv(SASP_genic,"SASP_genic.bed")

SASP_genic_mRNA<-mRNA_genic_labeled%>%filter(hgnc_symbol%in%SASP)%>%select(1:6)
write_tsv(SASP_genic_mRNA,"SASP_genic_mRNA.bed")


SignalP<-read_tsv("protein_class_SignalP.tsv")

Signal_genic<-cds_genic_labeled%>%filter(hgnc_symbol%in%SignalP$Gene)%>%select(1:6)
write_tsv(Signal_genic,"Signal_Genic.bed")

########################################################## Trans Bed

Trans_E_Bed<-cds_genic_labeled%>%filter(hgnc_symbol%in%Trans_Etop$hgnc_symbol)%>%select(1:6)
write_tsv(Trans_E_Bed,"Trans_E.bed")




