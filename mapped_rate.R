library(tidyverse)
library(ggplot2)

#draw density
fl=list.files(pattern=".tsv$")
file=gsub("\\.tsv","",fl)
Enum_list=gsub("E","",file)
mapped_rate=array(1:length(Enum_list))
total_mapped_base=array(1:length(Enum_list))
total_Base=array(1:length(Enum_list))
total_mapped_read_base=array(1:length(Enum_list))
mapped_read_rate=array(1:length(Enum_list))
all_data=read.table("../all_report.csv",sep=",",header=TRUE,colClasses=c("Enum"="character"))
for (i in seq(1,length(fl))){
  start_data = read.table(fl[i],sep="\t",header=TRUE)
  total_mapped_base[i]=sum(start_data$mapped_base)
  total_Base[i]=as.numeric(all_data[all_data["Enum"]==Enum_list[i]][3])
  mapped_rate[i]=total_mapped_base[i]/total_Base[i]*100
}
Enum=file
mapped_table=data.frame(Enum,total_mapped_base,total_Base,mapped_rate)
space_table= read.csv("./Enum_name.csv",header=TRUE)
diplot_name=c("anglica","rostkoviana","vigursii","montana","rivularis")
uncertain=c("reayensis")
test=left_join(mapped_table,space_table)
test["ploids"]="tetraploids"
for (i in diplot_name){
  test["ploids"]=ifelse(grepl(i,test$Species_name),"diploids",test$ploids)
}
for (i in uncertain){
  test["ploids"]=ifelse(grepl(i,test$Species_name),"unknown",test$ploids)
}

library(ggrepel)
ggplot(test, aes(Enum, mapped_rate)) +
  geom_point(aes(colour=ploids,size=2))  + 
  theme_bw() +
  ylab("Base mapped rate % ") +
  xlab("Sample") +
  geom_text(aes(label=Enum), vjust=-0.3, size=3.5)

            