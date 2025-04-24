args=commandArgs(T)
inputfile=args[1]
outpdffile=args[2]
library(ggplot2)
library(ggrepel)
df<-read.table(inputfile,header=T)
df$reference_chr<-substr(df$Rchr,6,10)
f_levels<- unique(df$reference_chr)
df$reference_chr<-factor(df$reference_chr,levels = f_levels)
reference_genome<-unique(substr(df$Rchr,1,4))
df$query_chr<-substr(df$Qchr,6,10)
query_name<-substr(df$Qchr,1,5)
df$subgenome=paste0(query_name,df$subg)
ggplot(data=df,aes(x=reference_chr,y=MashD))+geom_line(aes(group=subgenome),size=1.2,alpha=0.5)+
  geom_line(aes(group = reference_chr),linetype="dotted",size=0.3)+
  geom_point(aes(colour=subgenome,shape=subgenome),size=3,alpha=0.5)+
  geom_text_repel(aes(label=query_chr))+
  labs(x=reference_genome,y='Mash-distance')+
  scale_color_manual(values = c("red","blue"))+theme_bw()
ggsave(file=outpdffile, width = 220, height = 120, units = "mm")
