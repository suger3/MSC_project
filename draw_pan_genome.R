data=read.table(file="by_distance.record",header = TRUE)
#plot(x=data$file,y=data$num_seqs,"o")
library(ggplot2)

library(tidyverse)

library(ggrepel)

library(geosphere)
info=read.csv("CoordinatesPangenome.CSV",header = TRUE)
info$distance=0
for (i in seq(1:nrow(info))){
  info[i,12]=distm(c(info[i,8], info[i,9]), c(info[7,8], info[7,9]), fun = distHaversine)
}
data=left_join(data,info)
writeLines(data[order(data$Species_name),]$Enum,"by_sn.list")


data$Enum=factor(data$Enum)
data$sum_len=data$sum_len/1000000

##by order
ggplot(data,aes(x=seq(1:nrow(data)),y=sum_len))+
  geom_line()+
  geom_point(aes(color=Species_name),size=3)+
  geom_text_repel(aes(seq(1:nrow(data)),sum_len, label=Enum))+
  theme_classic()+
  #theme(text = element_text(size = 30))+
  #theme(axis.text.x = element_text(angle = 0,hjust = 0))+
  xlab("Sample Number")+
  ylab("Span(Mb)")

###by sn
ggplot(data,aes(x=distance,y=sum_len))+
  geom_line()+
  geom_point(aes(),size=2)+
  geom_text_repel(aes(distance,sum_len, label=Enum))+
  theme_classic()+
  #theme(text = element_text(size = 30))+
  #theme(axis.text.x = element_text(angle = 0,hjust = 0))+
  xlab("Distance to reference genome")+
  ylab("Span(Mb)")+
  scale_shape_manual(values =seq(1:15) )
