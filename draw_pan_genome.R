data=read.table(file="./depth/pan_genome51.record",header = TRUE)
#plot(x=data$file,y=data$num_seqs,"o")
library(ggplot2)
library(ggrepel)
name=data$file
#draw_span=data$sum_len
#data=as.data.frame(data)
data$file=factor(data$file,levels=name)
data$sum_len=data$sum_len/1000000
ggplot(data,aes(x=seq(1:nrow(data)),y=sum_len,group=1))+
  geom_line()+
  geom_point()+
  #geom_text_repel(aes(seq(1:nrow(data)),sum_len, label=file))+
  theme_classic()+
  #theme(text = element_text(size = 30))+
  #theme(axis.text.x = element_text(angle = 0,hjust = 0))+
  xlab("Sample Number")+
  ylab("Span(Mb)")
