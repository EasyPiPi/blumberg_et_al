
setwd("/local1/home/ablumber/K562/updated_data/kmer_count/temp_files/")

library(ggplot2)
library(cowplot)

fontSize=60
theme_set(theme_cowplot(font_size=fontSize))

final_data<-read.table("final_data_from_mRNA_3mer_step1.txt",header=T,sep="\t")

##########################
####
#### mRNA
####
#### change final plot based on plots number and names
####
##########################

final_data$position<-as.character(final_data$position)

temp<-aggregate(log_ratio~V2, data =final_data, FUN= "mean" )
temp<-subset(final_data,position=="0")
y<- factor(temp$V2[order(temp$log_ratio)])
y_3mer<-as.vector(y)
print(y_3mer)

cbPalette <- c( "#009E73", "#F0E442", "#0072B2", "#D55E00")

plot_final<-ggplot(final_data, aes(x=V2, ymin=0,ymax=log_ratio,color=position)) + geom_linerange(position = position_dodge(width = 0.5), size = 20, alpha=0.75) +  theme(legend.key.height = unit(2, "cm"))
plot_final<-plot_final + xlab("") + ylab("Log2(Stable/Unstable)" ) + theme(axis.text.x  = element_text(angle=45, vjust=0.5))
#plot_final<-plot_final +   scale_color_manual(values=cbPalette)  + ylim(range(lim_range))
plot_final<-plot_final +scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "500", "1000","1500"),labels=c("0-1000","500-1500","1000-2000","1500-2500"))+  scale_x_discrete(limits=y_3mer)
# scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "200", "400","600"),labels=c("0-400","200-600","400-800","600-1000"))+  scale_x_discrete(limits=dinuc)
plot_final<-plot_final+ geom_hline(yintercept=0, linetype="dashed") + theme(legend.justification=c(0,1),legend.position=c(0,1))

setwd("/local1/home/ablumber/K562/updated_data/kmer_count/temp_files/")

final_data_2mer<-read.table("final_data_from_mRNA_2mer_step1.txt",header=T,sep="\t")

##########################
####
#### change final plot based on plots number and names
####
##########################

final_data_2mer$position<-as.character(final_data_2mer$position)

temp<-aggregate(log_ratio~V2, data =final_data_2mer, FUN= "mean" )
temp<-subset(final_data_2mer,position=="0")
y<- factor(temp$V2[order(temp$log_ratio)])
y_2mer<-as.vector(y)
print(y_2mer)

cbPalette <- c( "#009E73", "#F0E442", "#0072B2", "#D55E00")

plot_final_2mer<-ggplot(final_data_2mer, aes(x=V2, ymin=0,ymax=log_ratio,color=position)) + geom_linerange(position = position_dodge(width = 0.5), size = 20, alpha=0.75) +  theme(legend.key.height = unit(2, "cm"))
plot_final_2mer<-plot_final_2mer + xlab("") + ylab("Log2(Stable/Unstable)" )+ theme(axis.text.x  = element_text(angle=45, vjust=0.5))
plot_final_2mer<-plot_final_2mer +scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "500", "1000","1500"),labels=c("0-1000","500-1500","1000-2000","1500-2500"))+  scale_x_discrete(limits=y_2mer)
# scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "200", "400","600"),labels=c("0-400","200-600","400-800","600-1000"))+  scale_x_discrete(limits=dinuc)
plot_final_2mer<-plot_final_2mer+ geom_hline(yintercept=0, linetype="dashed") + theme(legend.justification=c(0,1),legend.position=c(0,1))
plot_both<-plot_grid(plot_final_2mer, plot_final, labels = c('A', 'B'),ncol=1,label_size = 40)

setwd("/local1/home/ablumber/K562/updated_data/kmer_count/plots/")

save_plot("mRNA_updated_2_3MER_base_composition.pdf",plot_both,ncol=1,nrow=1, base_height=40, base_width=60,limitsize = FALSE)



##########################
####
###   lincRNA
####
##########################
setwd("/local1/home/ablumber/K562/updated_data/kmer_count/temp_files/")

final_data<-read.table("final_data_from_lincs_3mer_step1.txt",header=T,sep="\t")

final_data$position<-as.character(final_data$position)

plot_final<-ggplot(final_data, aes(x=V2, ymin=0,ymax=log_ratio,color=position)) + geom_linerange(position = position_dodge(width = 0.5), size = 20, alpha=0.75) +  theme(legend.key.height = unit(2, "cm"))
plot_final<-plot_final + xlab("") + ylab("Log2(Stable/Unstable)" ) + theme(axis.text.x  = element_text(angle=45, vjust=0.5))
#plot_final<-plot_final +   scale_color_manual(values=cbPalette)  + ylim(range(lim_range))
plot_final<-plot_final +scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "500", "1000","1500"),labels=c("0-1000","500-1500","1000-2000","1500-2500"))+  scale_x_discrete(limits=y_3mer)
# scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "200", "400","600"),labels=c("0-400","200-600","400-800","600-1000"))+  scale_x_discrete(limits=dinuc)
plot_final<-plot_final+ geom_hline(yintercept=0, linetype="dashed") + theme(legend.justification=c(0,1),legend.position=c(0,1))

setwd("/local1/home/ablumber/K562/updated_data/kmer_count/temp_files/")

final_data_2mer<-read.table("final_data_from_lincs_2mer_step1.txt",header=T,sep="\t")

final_data_2mer$position<-as.character(final_data_2mer$position)

plot_final_2mer<-ggplot(final_data_2mer, aes(x=V2, ymin=0,ymax=log_ratio,color=position)) + geom_linerange(position = position_dodge(width = 0.5), size = 20, alpha=0.75) +  theme(legend.key.height = unit(2, "cm"))
plot_final_2mer<-plot_final_2mer + xlab("") + ylab("Log2(Stable/Unstable)" )+ theme(axis.text.x  = element_text(angle=45, vjust=0.5))
plot_final_2mer<-plot_final_2mer +scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "500", "1000","1500"),labels=c("0-1000","500-1500","1000-2000","1500-2500"))+  scale_x_discrete(limits=y_2mer)
# scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "200", "400","600"),labels=c("0-400","200-600","400-800","600-1000"))+  scale_x_discrete(limits=dinuc)
plot_final_2mer<-plot_final_2mer+ geom_hline(yintercept=0, linetype="dashed") + theme(legend.justification=c(0,1),legend.position=c(0,1))
plot_both<-plot_grid(plot_final_2mer, plot_final, labels = c('A', 'B'),ncol=1,label_size = 40)

setwd("/local1/home/ablumber/K562/updated_data/kmer_count/plots/")
save_plot("lincRNA_updated_2_3MER_base_composition.pdf",plot_both,ncol=1,nrow=1, base_height=40, base_width=60,limitsize = FALSE)




##########################
####
###   lincRNA
####
##########################
setwd("/local1/home/ablumber/K562/updated_data/kmer_count/temp_files/")

final_data<-read.table("final_data_from_eRNA_3mer_step1.txt",header=T,sep="\t")

final_data$position<-as.character(final_data$position)

plot_final<-ggplot(final_data, aes(x=V2, ymin=0,ymax=log_ratio,color=position)) + geom_linerange(position = position_dodge(width = 0.5), size = 20, alpha=0.75) +  theme(legend.key.height = unit(2, "cm"))
plot_final<-plot_final + xlab("") + ylab("Log2(Stable/Unstable)" ) + theme(axis.text.x  = element_text(angle=45, vjust=0.5))
#plot_final<-plot_final +   scale_color_manual(values=cbPalette)  + ylim(range(lim_range))
#plot_final<-plot_final +scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "500", "1000","1500"),labels=c("0-1000","500-1500","1000-2000","1500-2500"))+  scale_x_discrete(limits=y_3mer)
plot_final<-plot_final + scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "200", "400","600"),labels=c("0-400","200-600","400-800","600-1000"))+  scale_x_discrete(limits=y_3mer)
plot_final<-plot_final+ geom_hline(yintercept=0, linetype="dashed") + theme(legend.justification=c(1,1),legend.position=c(1,1))

setwd("/local1/home/ablumber/K562/updated_data/kmer_count/temp_files/")

final_data_2mer<-read.table("final_data_from_eRNA_2mer_step1.txt",header=T,sep="\t")

final_data_2mer$position<-as.character(final_data_2mer$position)

plot_final_2mer<-ggplot(final_data_2mer, aes(x=V2, ymin=0,ymax=log_ratio,color=position)) + geom_linerange(position = position_dodge(width = 0.5), size = 20, alpha=0.75) +  theme(legend.key.height = unit(2, "cm"))
plot_final_2mer<-plot_final_2mer + xlab("") + ylab("Log2(Stable/Unstable)" )+ theme(axis.text.x  = element_text(angle=45, vjust=0.5))
#scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "500", "1000","1500"),labels=c("0-1000","500-1500","1000-2000","1500-2500"))+  scale_x_discrete(limits=y_2mer)
plot_final_2mer<-plot_final_2mer + scale_color_discrete(name="Window positons\n(downstream to TSS)", limits=c("0", "200", "400","600"),labels=c("0-400","200-600","400-800","600-1000"))+  scale_x_discrete(limits=y_2mer)
plot_final_2mer<-plot_final_2mer+ geom_hline(yintercept=0, linetype="dashed") + theme(legend.justification=c(1,1),legend.position=c(1,1))
plot_both<-plot_grid(plot_final_2mer, plot_final, labels = c('A', 'B'),ncol=1,label_size = 40)

setwd("/local1/home/ablumber/K562/updated_data/kmer_count/plots/")
save_plot("eRNA_updated_2_3MER_base_composition.pdf",plot_both,ncol=1,nrow=1, base_height=40, base_width=60,limitsize = FALSE)
