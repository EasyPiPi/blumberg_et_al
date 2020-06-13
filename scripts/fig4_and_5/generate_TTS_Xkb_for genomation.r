
setwd("/mnt/grid/siepel/hpc_norepl/home/data/ablumber/matched_data")
files_list<-list.files(pattern = "*.bed")
win=1000

 for (i in files_list){
 myfile<-read.table(i,sep="\t",header=F)
 plus_file<-subset(myfile,V6=="+")

 plus_file$V2<-plus_file$V3-win
 plus_file$V3<-plus_file$V3+win
 name<-paste("plus_TTS_",i,"")
   write.table(plus_file,name,quote = F, col.names = T, row.names = F,sep="\t")
    minus_file<-subset(myfile,V6=="-")

 minus_file$V3<-minus_file$V2+win
 minus_file$V2<-minus_file$V2-win
 name<-paste("minus_TTS_",i,"")
   write.table(minus_file,name,quote = F, col.names = T, row.names = F,sep="\t")
 }
