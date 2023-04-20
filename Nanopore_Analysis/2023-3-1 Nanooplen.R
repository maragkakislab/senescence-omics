library(tidyverse)
library(writexl)

EC_Nanoplen_URL<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/nanoplen/results/CYC23-Sen23/transcripts/diff_len_cyc23-sen23_LMM-test_log2.tab"

EC_Nanoplen<-EC_Nanoplen_URL%>%
  map_dfr(read_tsv,col_names=T)%>%
  na.omit()

EC_Nanoplen<-EC_Nanoplen%>%filter(name%in%EC_ENST_Gene$tid)%>%
  filter(n.control>=5 & n.alt>=5)

ggplot()+
  geom_point(data=EC_Nanoplen,
             aes(x=log2FC,y=-log10(qvalue)),alpha=0.6)+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  xlab("Log2FC Transcript Length (Etop/Cyc)")+ ylab("-log10 q-value")+
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),
        axis.text.y = element_text(size=24,color="black"),
        axis.title.x = element_text(size=20,color="black"),
        axis.title.y = element_text(size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))


EC_SigShort<-EC_Nanoplen%>%filter(log2FC<= -0.25& pvalue<=0.05)

values<-EC_SigShort$name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
x<-getBM(attributes = c('ensembl_transcript_id', "description","hgnc_symbol","gene_biotype"),
         filters = 'ensembl_transcript_id',
         values = values, 
         mart = ensembl)


EC_SigShort<-merge(EC_SigShort,x,by.x="name",by.y="ensembl_transcript_id")

EC_ENST_ShortFold<-EC_ENST_Gene%>%filter(tid%in%EC_SigShort$name)

plot(EC_SigShort$log2FC,EC_ENST_ShortFold$log2FoldChange)

write_xlsx(EC_SigShort,"EC_Sigshort.xlsx")



############################################################################## E/C Poly A

EC_PolyA_url<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/nanoplen_polyA/results/cyc23-sen23/polyA-diff_len_cyc23-sen23-LMM.tab"

EC_PolyA<-EC_PolyA_url%>%
  map_dfr(read_tsv,col_names=T)%>%
  na.omit()

EC_PolyA<-EC_PolyA%>%filter(name%in%EC_ENST_Gene$tid)%>%
  filter(n.control>=5 & n.alt>=5)


ggplot()+
  geom_point(data=EC_PolyA,
             aes(x=log2FC,y=-log10(qvalue)),alpha=0.6)+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  xlab("Log2FC Poly-A Tail Length (Etop/Cyc)")+ ylab("-log10 q-value")+
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),
        axis.text.y = element_text(size=24,color="black"),
        axis.title.x = element_text(size=20,color="black"),
        axis.title.y = element_text(size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))




############################################################################## C/CTg

CTgC_Nanoplen_URL<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/nanoplen/results/ctgab_cyc23/transcripts/diff_len_ctgab_cyc23_LMM-test_log2.tab"

CTgC_Nanoplen<-CTgC_Nanoplen_URL%>%
  map_dfr(read_tsv,col_names=T)%>%
  na.omit()


CTgc_Nanoplen<-CTgC_Nanoplen%>%filter(name%in%CTgC_ENST_Gene$tid)%>%
  filter(n.control>=5 & n.alt>=5)

ggplot()+
  geom_point(data=CTgc_Nanoplen,
             aes(x=log2FC,y=-log10(qvalue)),alpha=0.6)+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  xlab("Log2FC Length (CTg/Cyc)")+ ylab("-log10 q-value")+
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),
        axis.text.y = element_text(size=24,color="black"),
        axis.title.x = element_text(size=20,color="black"),
        axis.title.y = element_text(size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))



############################################################################ ETg/CTg



CTgETg_Nanoplen_URL<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/nanoplen/results/ctgab_etgab/transcripts/diff_len_ctgab_etgab_LMM-test_log2.tab"
CTgETg_Nanoplen<-CTgETg_Nanoplen_URL%>%
  map_dfr(read_tsv,col_names=T)%>%
  na.omit()

CTgETg_Nanoplen<-CTgETg_Nanoplen%>%
  filter(n.control>=5 & n.alt>=5)


ggplot()+
  geom_point(data=CTgETg_Nanoplen,
             aes(x=log2FC,y=-log10(qvalue)),alpha=0.6)+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  xlab("Log2FC Length (CTg/Cyc)")+ ylab("-log10 q-value")+
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),
        axis.text.y = element_text(size=24,color="black"),
        axis.title.x = element_text(size=20,color="black"),
        axis.title.y = element_text(size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))



############################################################################ ETg/E

ETgE_URL<-"https://niairplggcgu.irp.nia.nih.gov/shared/pIiIoFaWFANXoDZN3JoBH24AmvzxOD0c/cyc_sen/analysis/nanoplen/results/etgab_sen23/transcripts/diff_len_etgab_sen23_LMM-test_log2.tab"

ETgE_Nanoplen<-ETgE_URL%>%
  map_dfr(read_tsv,col_names=T)%>%
  na.omit()

ETgE_Nanoplen<-ETgE_Nanoplen%>%
  filter(n.control>=5 & n.alt>=5)


ggplot()+
  geom_point(data=ETgE_Nanoplen,
             aes(x=log2FC,y=-log10(qvalue)),alpha=0.6)+
  geom_hline(yintercept = -log10(0.05),color="blue",size=2,linetype=2)+
  xlab("Log2FC Length (ETg/Etop)")+ ylab("-log10 q-value")+
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),
        axis.text.y = element_text(size=24,color="black"),
        axis.title.x = element_text(size=20,color="black"),
        axis.title.y = element_text(size=20,color="black"),
        axis.line = element_line(size=1.5),axis.ticks = element_line(size=1.5),
        axis.ticks.length = unit(0.2,"cm"))

