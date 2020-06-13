library(ggplot2)
library(twoBit)
library(stabilityHMM)
library(rqhmm)
setwd("/HMM")
library(cowplot)

fontSize=40
theme_set(theme_cowplot(font_size=fontSize))
#########
###   design each set of gene as group 
### here example for Intron-containing and intron-less mRNA and lincRNA
###
### each of "mRNA_spliced" ; "mRNA_unspliced" ; "linc_spliced"  ; "linc_unspliced" is set of genes in bed format. 

group1<-mRNA_spliced
group2<-mRNA_unspliced
group3<-linc_spliced
group4<-linc_unspliced

group1$type<-"mRNA"
group2$type<-"mRNA"
group3$type<-"lincRNA"
group4$type<-"lincRNA"

group1$status<-"Intron-containing"
group2$status<-"Intron-less"
group3$status<-"Intron-containing"
group4$status<-"Intron-less"


mydata<-rbind(group1,group2,group3,group4)

seqLen=1000 
hg19 <-"/sonas-hs/siepel/nlsas/data/data/home/ablumber/genomes/hg19.2bit" 

seqs <- collectSequences(hg19, mydata, seq.length = seqLen)
mData <- prepareData(seqs)
stable_score_spliced <- unstableScore(mData)
mydata$HMM<-stable_score_spliced

instability_Score<-mydata

instability_Score<-mydata
instability_Score$HMM=1-instability_Score$HMM

pd <- position_dodge(0.9) # move them .05 to the left and right

HMM3<-ggplot(instability_Score, aes(x=type, y=HMM, fill=status)) + geom_violin(alpha=0.5, draw_quantiles=0.5) + ylab("SSI") + xlab("") 
HMM3<-HMM3 + theme(legend.title=element_blank())
HMM3<-HMM3 + guides(fill=guide_legend( keywidth=0.5, keyheight=0.8, default.unit="inch"))

save_plot("HMM_K562_mRNS_spliced_unspliced.pdf",HMM3,base_width=16,base_height=8)
