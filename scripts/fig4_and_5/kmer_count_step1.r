
library(twoBit)
#library(plyr)
#library(ggplot2)
#library(cowplot)


collect_sequences <- function(twobit.filename, bed, seq.length = 1000,start_pos=0) {
  twobit = twoBit::twobit.load(path.expand(twobit.filename))

  N = dim(bed)[1]
  result = vector(mode="list", length = N)

  is.minus = bed[,6] == '-'
  starts = bed[,2] +start_pos
  ends = bed[,2] + seq.length +start_pos
  starts[is.minus] = bed[is.minus, 3] - seq.length - start_pos
  ends[is.minus] = bed[is.minus, 3] - start_pos


  chroms = as.character(bed[,1])

  for (i in 1:N) {
    chrom = chroms[i]

    seq = twoBit::twobit.sequence(twobit, chrom, starts[i], ends[i])
    if (is.minus[i]){
      seq = twoBit::twobit.reverse.complement(seq)
    }
    result[i] <- seq
  }

  return(result)
}

final_data<-data.frame()
range_matrix<-data.frame()

bases<-c("A","G","C","T")

hexamers_3mer=vector()

#########################################################
####
#### Generate K-mer index
#### Change for the proper K-mer length
#### number of loops
#### and 'paste' line
####
#########################################################

for (a in bases) {
  for (b in bases) {
    for (c in bases) {
#      for (d in bases) {
#        for (e in bases) {
#          for (f in bases) {
#            hexamer <- paste(a,b,c,d,e,f,sep="")
            hexamer <- paste(a,b,c,sep="")
            hexamers_3mer<-rbind(hexamers_3mer,hexamer)

#          }
 #       }
  #    }
    }
  }
}

hexamers_2mer=vector()

for (a in bases) {
  for (b in bases) {
            hexamer <- paste(a,b,sep="")
            hexamers_2mer<-rbind(hexamers_2mer,hexamer)
  }
}

prepare_kmer_data<-function(stable_genes,unstable_genes,pos_start=0,pos_end=1500,steps=500,seqLen=1000,hexamers_in){

seqLen=seqLen
hexamers=hexamers_in
pos_start=pos_start
pos_end=pos_end
steps=steps
FILE_stable=stable_genes
FILE_UNstable=unstable_genes

hg19 <-"/sonas-hs/siepel/nlsas/data/home/ablumber/genomes/hg19.2bit"
lim_range=vector()


for (i in seq(pos_start,pos_end,steps)){
genes<-read.table(FILE_stable)
nam<-paste("plot_",i,sep="")
N=nrow(genes)
seqs<- collect_sequences(hg19, genes, seq.length = seqLen ,start_pos=i)

results_sense<-data.frame()
results_position=1
N=length(seqs)

for (seq_test in seqs){
count_position<-0
hexamer_count=vector()
hexamer_count_as=vector()

#print(results_position)


for (hexamer in hexamers) {
            count=0
            x <- seq_test
            m <- gregexpr(paste("(?=(", hexamer, "))", sep=""), x, perl=TRUE)
            m <- lapply(m, function(i) {
    attr(i,"match.length") <- attr(i,"capture.length")
    i
})

            count_sense <- length(regmatches(x,m)[[1]])
       #     hexamer_count[count_position]<-count
            count_position=count_position+1


            hexamer_count[count_position]<-count_sense

#if (sum(hexamer_count)<=(nchar(seq_test)-(nchar(hexamer))+1)) {
}
results_sense<-rbind(results_sense,hexamer_count)
#}
}


#write.table(results_sense,"hexamer_lincRNA_K562_sense_raw_data.txt",quote=F,row.names=F,col.names=F,sep="\t")

col.sum.stable<-apply(results_sense, 2 ,sum)
N=nrow(results_sense)
col.sum.stable<-col.sum.stable/N
print(N)
print(sum(col.sum.stable))
#col.sum.stable<-as.data.frame(col.sum.stable)
names(col.sum.stable)<-hexamers

genes<-read.table(FILE_UNstable)
n=nrow(genes)
seqs<- collect_sequences(hg19, genes,seq.length = seqLen ,start_pos=i )

results_antisense<-data.frame()
results_position=1
N=length(seqs)

for (seq_test in seqs){
count_position<-0
hexamer_count=vector()
hexamer_count_as=vector()

#print(results_position)


for (hexamer in hexamers) {
            count=0
            x <- seq_test
            m <- gregexpr(paste("(?=(", hexamer, "))", sep=""), x, perl=TRUE)
            m <- lapply(m, function(i) {
    attr(i,"match.length") <- attr(i,"capture.length")
    i
})

            count_sense <- length(regmatches(x,m)[[1]])
       #     hexamer_count[count_position]<-count
            count_position=count_position+1


            hexamer_count[count_position]<-count_sense

#if (sum(hexamer_count)<=(nchar(seq_test)-(nchar(hexamer))+1)) {
}
results_antisense<-rbind(results_antisense,hexamer_count)
#}
}


#write.table(results_antisense,"hexamer_lincRNA_K562_sense_raw_data.txt",quote=F,row.names=F,col.names=F,sep="\t")

col.sum.unstable<-apply(results_antisense, 2 ,sum)
print(sum(col.sum.unstable))
N=nrow(results_antisense)
col.sum.unstable<-col.sum.unstable/N
print(N)
print(sum(col.sum.unstable))
#col.sum.unstable<-as.data.frame(col.sum.unstable)
names(col.sum.unstable)<-hexamers
col.sum.stable<-as.data.frame(col.sum.stable)
col.sum.unstable<-as.data.frame(col.sum.unstable)

mydata<-cbind(col.sum.stable,col.sum.unstable)

mydata$ratio<-mydata$col.sum.stable/mydata$col.sum.unstable
mydata$V2<-hexamers
mydata$log_ratio<-log(mydata$ratio,2)
mydata<-mydata[order(mydata$log_ratio),]
mydata$V2 <- factor(mydata$V2, levels = mydata$V2[order(mydata$log_ratio)])
mydata$position<-i
lim_range<-c(lim_range, range(mydata$log_ratio))
temp<-range(mydata$log_ratio)
temp[3]<-i
range_matrix<-rbind(range_matrix,temp)

final_data<-rbind(final_data,mydata)
#p1<-ggplot(mydata, aes(x=V2,y=log_ratio)) + geom_point(size=6,color="blue")
#+ geom_bar(stat="identity")
#p1<-p1+xlab("") + ylab("Stable/Unstable (log scale) " ) + ylim(range(lim_range)) + theme(axis.text.x  = element_text(angle=45, vjust=0.5))
#assign(nam,p1)
print(i)
print("Hello")
}
return(final_data)

}

###
##
#
#
#      mRNA
#

FILE_stable<-"/local1/home/ablumber/K562/updated_data/K562_updated_matched_spliced_protein_class_5.bed"

FILE_UNstable<-"/local1/home/ablumber/K562/updated_data/K562_updated_matched_spliced_protein_class_1.bed"


setwd("/local1/home/ablumber/K562/updated_data/kmer_count/temp_files/")
final_data_3mer<-prepare_kmer_data(stable_genes=FILE_stable,unstable_genes=FILE_UNstable,hexamers_in=hexamers_3mer)
write.table(final_data_3mer,"final_data_from_mRNA_3mer_step1.txt",quote=F,row.names=F,col.names=T,sep="\t")
write.table(range_matrix,"range_matrix_from_step1.txt",quote=F,row.names=F,col.names=T,sep="\t")
final_data_2mer<-prepare_kmer_data(stable_genes=FILE_stable,unstable_genes=FILE_UNstable,hexamers_in=hexamers_2mer)
write.table(final_data_2mer,"final_data_from_mRNA_2mer_step1.txt",quote=F,row.names=F,col.names=T,sep="\t")


FILE_stable<-"/local1/home/ablumber/K562/updated_data/K562_updated_matched_spliced_lincRNA_class_5.bed"

FILE_UNstable<-"/local1/home/ablumber/K562/updated_data/K562_updated_matched_spliced_lincRNA_class_1.bed"

final_data_3mer<-prepare_kmer_data(stable_genes=FILE_stable,unstable_genes=FILE_UNstable,hexamers_in=hexamers_3mer)
write.table(final_data_3mer,"final_data_from_lincs_3mer_step1.txt",quote=F,row.names=F,col.names=T,sep="\t")
final_data_2mer<-prepare_kmer_data(stable_genes=FILE_stable,unstable_genes=FILE_UNstable,hexamers_in=hexamers_2mer)
write.table(final_data_2mer,"final_data_from_lincs_2mer_step1.txt",quote=F,row.names=F,col.names=T,sep="\t")

FILE_stable<-"/local1/home/ablumber/CAGE/TSS_data/final_stable_k562_high_cage_10.srt.mrg.bed"
FILE_UNstable<-"/local1/home/ablumber/CAGE/TSS_data/final_UNstable_k562_high_CAGE_10.match.srt.mrg.bed"


setwd("/local1/home/ablumber/K562/updated_data/kmer_count/temp_files/")
final_data_3mer<-prepare_kmer_data(stable_genes=FILE_stable,unstable_genes=FILE_UNstable,pos_start=0,pos_end=600,steps=200, seqLen=400, hexamers_in=hexamers_3mer)
write.table(final_data_3mer,"final_data_from_eRNA_3mer_step1.txt",quote=F,row.names=F,col.names=T,sep="\t")
final_data_2mer<-prepare_kmer_data(stable_genes=FILE_stable,unstable_genes=FILE_UNstable,pos_start=0,pos_end=600,steps=200,seqLen=400, hexamers_in=hexamers_2mer)
write.table(final_data_2mer,"final_data_from_eRNA_2mer_step1.txt",quote=F,row.names=F,col.names=T,sep="\t")
