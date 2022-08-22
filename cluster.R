contigs=as.data.frame(readRDS("all_contigs.Rds"))
#info=readRDS("all_info.Rds")
info=read.csv("CoordinatesPangenome.CSV",header = TRUE)
library(tidyverse)
library(ggplot2)
library(pheatmap2)
library(fastcluster)
kmeans_cluster=read.csv("kmeans_cluster.csv")
colnames(kmeans_cluster)=c("contigs","Kc")
table(kmeans_cluster$Kc)
#     1      2      3 
#265917 449391 220124 
contigs$contigs=row.names(contigs)
alldata=left_join(contigs,kmeans_cluster)
len=readRDS("length.rds")
colnames(len)=c("contigs","length")
alldata=left_join(alldata,len)

all_0=which(rowSums(contigs[1:52]<0.2) ==52)
all_0=alldata[all_0,]
sum(all_0$length)
writeLines(all_0$contigs,con="real_all_0.txt",sep="\n")

##by kmeans
#Agenome
cluster3=alldata %>%
  filter(Kc==3)
sum(cluster3$length)
#143985929

#Bgenome
cluster1=alldata %>%
  filter(Kc==1)
sum(cluster1$length)
#437354541

#pav
cluster2=alldata %>%
  filter(Kc==2)
sum(cluster2$length)
#401535949
#the Bgenome are far away from what I predicted
#########################################################


##by specific value

pre_1000=readRDS("all_contigs.Rds")

have=0.6
nhave=0.2

#per_1000_data=read.csv("heatmap_gt4999_per_5000.csv",header = TRUE)

pre_1000=filter(alldata,length>4999)
#heatmap_order=table[test$tree_row$order,]
#heatmap_name_order=as.data.frame(row.names(heatmap_order))
#colnames(heatmap_name_order)=c("contigs")
#table$contigs=rownames(table)
row_order=readRDS("heatmap_order.rds")
table_tmp=as.data.frame(row_order)
colnames(table_tmp)=c("contigs")
pre_1000=left_join(table_tmp,pre_1000)

name=pre_1000$contigs
pre_1000=pre_1000[,2:53]
pre_1000=apply(pre_1000,2,as.numeric)
pre_1000=pre_1000[which(rowSums(pre_1000 >0) >0),]
#pre_heatmap_Data=pre_heatmap_Data[which(rowSums(pre_heatmap_Data) > 0),]#get with depth not all 0
pre_1000=as.matrix(pre_1000)
rownames(pre_1000)=name

annotation_row <- data.frame(pre_1000[,1:2])
colnames(annotation_row)=c("like_sf","nature")
annotation_row$nature="uncertain"

PAV_index1=which(rowSums(pre_1000>have) >= 1) #at least 1 samples have
PAV_index2=which(rowSums(pre_1000>have) < 15) #at most 14 samples have
#PAV_index3=which(rowSums(pre_1000>nhave) >= 1) #at least 51 samples have
#PAV_index4=which(rowSums(pre_1000<nhave) <= 8) #at least 51 samples have
#Ageome_index2=which(rowSums(pre_1000<nhave) <= 1) #at most 1 samples have not
#PAV_index5=intersect(PAV_index1,PAV_index2)
PAV_index=intersect(PAV_index1,PAV_index2)
annotation_row[PAV_index,]$nature="PAV"
#sum(pre_1000[PAV_index,]$length)
#288723125    1<=have<51 other uncertain
#372631743    1<=have<51 other PAV


annotation_row$contigs=row.names(annotation_row)
annotation_row=left_join(annotation_row,len)
sum(as.numeric(filter(annotation_row,nature=="PAV")$length))
#PAV 126659301
sum(as.numeric(filter(annotation_row,nature=="Agenome")$length))
#Agenome 168002413
sum(as.numeric(filter(annotation_row,nature=="Bgenome")$length))
#Bgenome 137183897

Ageome_index1=which(rowSums(pre_1000>have) >= 50) #at least 51 samples have
#Ageome_index2=which(rowSums(pre_1000<nhave) <= 1) #at most 1 samples have not
#Ageome_index=intersect(Ageome_index1,Ageome_index2)
annotation_row[Ageome_index1,]$nature="Agenome"


#15 only have A,  37 have A and B
Bgeome_index1=which(rowSums(pre_1000>have) <= 42) #at least 37 samples have
Bgeome_index2=which(rowSums(pre_1000>have) >= 32) #at least 37 samples have
#Bgeome_index3=which(rowSums(pre_1000<nhave) > 10) #at least 10 samples have
#Bgeome_index4=which(rowSums(pre_1000<nhave) < 20) #at least 51 samples have
Btmp=intersect(Bgeome_index1,Bgeome_index2)
#Btmp1=intersect(Btmp,Bgeome_index3)
#Bgeome_index=intersect(Btmp1,Bgeome_index4)
annotation_row[Btmp,]$nature="Bgenome"


annotation_row$like_sf="uncertain"
core_index1=which(rowSums(pre_1000>have) >= 51) #at least 51 samples have
core_index2=which(rowSums(pre_1000<nhave) <= 1) #at most 1 samples have not
core_index=intersect(core_index1,core_index2)
annotation_row[core_index,]$like_sf="core"

rare_index1=which(rowSums(pre_1000<nhave) >= 51)# at least 51 samples have not
rare_index2=which(rowSums(pre_1000>have) >= 1 )# at most 1 samples have
rare_index=intersect(rare_index1,rare_index2)
annotation_row[rare_index,]$like_sf="rare"

dispensable_index1=which(rowSums(pre_1000>have) < 51) #less 51 samples have
dispensable_index2=which(rowSums(pre_1000>have) > 1 ) #over  1 samples have
dispensable_index3=which(rowSums(pre_1000<nhave) < 51 )#less 51 samples have not
dispensable_index4=which(rowSums(pre_1000<nhave) > 1 ) #over 1 samples have not
temp_index=intersect(dispensable_index1,dispensable_index2)
temp_index=intersect(temp_index,dispensable_index3)
temp_index=intersect(temp_index,dispensable_index4)
annotation_row[temp_index,]$like_sf="dispensable"

annotation_col <- data.frame(info[c("uncertain_ploids")])
colnames(annotation_col) = c("ploids")
rownames(annotation_col) <- colnames(pre_1000)
pheatmap2(pre_1000,scale = "none",show_rownames=FALSE,treeheight_row = 0,angle_col ="0",annotation_col = annotation_col,annotation_row = annotation_row)

##############################################################################################
#pan-genome
info=as.data.frame(read.csv("CoordinatesPangenome.CSV",header = TRUE)[1:11])
library(geosphere)
sn_order=info[order(info[,5]),]
info$distance=0
for (i in seq(1:nrow(info))){
  info[i,12]=distm(c(info[i,8], info[i,9]), c(info[7,8], info[7,9]), fun = distHaversine)
}
distance_order=info[order(info[,12]),]

writeLines(distance_order$E.number)

#############################################################################################
#heatmap
contigs_nature=readRDS("contigs_nature.rds")
draw_heatmap= function(table,picture_name="plot.png"){
  col_num=ncol(table)
  table=as.data.frame(table)
  table$contigs=row.names(table)
  table=left_join(table,contigs_nature)
  table=left_join(table,len)
  #table = filter(table,nature=="PAV"|nature=="Bgenome")
  table=table[which(rowSums(table[,1:col_num] >0.6) >0),]
  
  name=table$contigs
  sum_length=sum(table$length)
  #table=table[,1:52]
  table=apply(table[,1:col_num],2,as.numeric)
  #pre_heatmap_Data=pre_heatmap_Data[which(rowSums(pre_heatmap_Data) > 0),]#get with depth not all 0
  table=as.matrix(table)
  rownames(table)=name
  
  annotation_row <- data.frame(table[,1])
  colnames(annotation_row)=c("contigs")
  annotation_row$contigs=row.names(annotation_row)
  annotation_row=left_join(annotation_row,contigs_nature)
  name=annotation_row$contigs
  annotation_row=as.data.frame(annotation_row[,2])
  colnames(annotation_row)=c("type")
  row.names(annotation_row)=name
  bk <- c(seq(0,0.59,by=0.01),seq(0.6,1,by=0.01),seq(1.01,4,by=0.01))
  annotation_col <- data.frame(info[c("Species_name")])
  rownames(annotation_col) <- info$Enum
  ann_colors = list(
    type = c(Agenome="green", PAV="black",Bgenome="red",uncertain="lightblue")
  )
  png(paste0(picture_name,"_ALL_in_paper.png"),width = 800,height=800)
  pheatmap2(table,scale = "none",show_rownames=FALSE,treeheight_row = 0,angle_col ="0",
            annotation_row = annotation_row,cluster_rows = F,
            annotation_col = annotation_col,
            #color = c(rep("#000000",60),colorRampPalette(colors = c("white","yellow"))(41),colorRampPalette(colors = c("yellow","red"))(300)),
            #legend_breaks=seq(0,4,0.5),
            #breaks=bk,
            fontsize = 12,
            annotation_names_row = FALSE,
            annotation_colors = ann_colors)
  dev.off()
}
############################save heatmap_order



annotation_col <- data.frame(info["genetype"])
rownames(annotation_col) <- info$Enum
#test = pheatmap2(pre_1000,scale = "none",treeheight_row = 0,show_rownames=FALSE,
#          angle_col ="0",annotation_col = annotation_col,annotation_row=annotation_row,
#          annotation_colors = ann_colors)


pheatmap2(pre_1000,scale = "none",treeheight_row = 0,show_rownames=FALSE, cluster_rows = F,
                 angle_col ="0",annotation_col = annotation_col,annotation_row=annotation_row,
                 annotation_colors = ann_colors)
draw_PAV_heatmap= function(table,picture_name="plot.png"){
  col_num=ncol(table)
  table=as.data.frame(table)
  table$contigs=row.names(table)
  table=left_join(table,contigs_nature)
  table=left_join(table,len)
  table = filter(table,nature=="PAV")
  table=table[which(rowSums(table[,1:col_num] >0.6) >0),]
  
  name=table$contigs
  sum_length=sum(table$length)
  #table=table[,1:52]
  table=apply(table[,1:col_num],2,as.numeric)
  
  
  #pre_heatmap_Data=pre_heatmap_Data[which(rowSums(pre_heatmap_Data) > 0),]#get with depth not all 0
  table=as.matrix(table)
  rownames(table)=name
  
  #colnames(annotation_row)=c("contigs")
  #annotation_row$contigs=row.names(annotation_row)
  #annotation_row=left_join(annotation_row,contigs_nature)
  #name=annotation_row$contigs
  #annotation_row=as.data.frame(annotation_row[,2])
  #colnames(annotation_row)=c("nature")
  #row.names(annotation_row)=name
  annotation_col =info[info$Enum%in%colnames(table),]
  name=annotation_col$Enum
  annotation_col$distance=-1
  for (i in seq(1:nrow(annotation_col))){
    annotation_col[i,12]=distm(c(annotation_col[i,8], annotation_col[i,9]), c(annotation_col[1,8], annotation_col[1,9]), fun = distHaversine)
  }
  annotation_col=as.data.frame(annotation_col[,c(5)])
  colnames(annotation_col)=c("Species_name")
  rownames(annotation_col) <- name
  bk <- c(seq(0,0.59,by=0.01),seq(0.6,1,by=0.01),seq(1.01,4,by=0.01))
  png(paste0(picture_name,"_PAV.png"),width = 1500,height=800)
  pheatmap2(table,scale = "none",show_rownames=FALSE,treeheight_row = 0,
            angle_col ="0",
            annotation_col = annotation_col,
            #color = c(rep("#000000",60),colorRampPalette(colors = c("white","yellow"))(41),colorRampPalette(colors = c("yellow","red"))(300)),
            #legend_breaks=c(0,0.6,1,seq(1.5,4,0.5)),
            #breaks=bk,
            #border_color="white",
            fontsize = 20,
            cluster_rows = F,
            annotation_names_row = FALSE,
            main = paste0("heatmap length:",sum_length)
            )
  #mtext(,side = 4,outer=TRUE)
  dev.off()
}

#sn_order
anglica=filter(sn_order,Species_name=="E. anglica")
table=pre_1000[,anglica$Enum]
#draw_PAV_heatmap(table)
draw_heatmap(table,"anglica")
#pheatmap2(pre_1000[,anglica$E.number],scale = "none",show_rownames=FALSE,treeheight_row = 0,angle_col ="0",annotation_col = annotation_col,annotation_row = annotation_row)

arctica=filter(sn_order,Species_name=="E. arctica")
table=pre_1000[,arctica$Enum]
#draw_PAV_heatmap(table)
draw_heatmap(table,"arctica")

arctica_septentrionalis=filter(sn_order,Species_name=="E. arctica"|Species_name=="E. septentrionalis")
table=pre_1000[,arctica_septentrionalis$Enum]
draw_PAV_heatmap(table,"arctica_septentrionalis")
draw_heatmap(table,"arctica_septentrionalis")


foulaensis=filter(sn_order,Species_name=="E. foulaensis"|Species_name=="E. confusa/foulaensis"|Species_name=="E. foulaensis/marshallii")
table=pre_1000[,foulaensis$Enum]
#draw_PAV_heatmap(table)
draw_heatmap(table,"foulaensis")

micrantha=filter(info,Species_name=="E. micrantha")
table=pre_1000[,micrantha$Enum]
draw_heatmap(table,"micrantha")
#draw_PAV_heatmap(table,"micrantha")


glanduligera_rostkoviana_anglica=filter(sn_order,Species_name=="E. glanduligera"|Species_name=="E. rostkoviana"|Species_name=="E. arctica")
draw_PAV_heatmap(pre_1000[,glanduligera_rostkoviana_anglica$Enum],"glanduligera_rostkoviana_anglica")
draw_heatmap(pre_1000[,glanduligera_rostkoviana_anglica$Enum],"glanduligera_rostkoviana_anglica")

rostkoviana=filter(sn_order,Species_name=="E. rostkoviana")
draw_heatmap(pre_1000[,rostkoviana$Enum],"rostkoviana")
draw_PAV_heatmap(pre_1000[,rostkoviana$Enum],"rostkoviana")

rivularis=filter(sn_order,Species_name=="E. rivularis")
draw_heatmap(pre_1000[,rivularis$Enum],"rivularis")
#draw_PAV_heatmap(pre_1000[,rivularis$Enum])

other=filter(sn_order,Species_name=="E. scottica"|Species_name=="E. montana"|Species_name=="E. septentrionalis"|Species_name=="E. vigursii"|Species_name=="unknow"|Species_name=="E. reayensis")
draw_heatmap(pre_1000[,other$Enum],"other")
#draw_PAV_heatmap(pre_1000[,other$Enum])

vigursii_rivularis=filter(sn_order,Species_name=="E. rivularis"|Species_name=="E. vigursii")
draw_heatmap(pre_1000[,vigursii_rivularis$Enum],"vigursii_rivularis")
draw_PAV_heatmap(pre_1000[,vigursii_rivularis$Enum],"vigursii_rivularis")
draw_PAV_heatmap(pre_1000,"all_sample")
#distance order
draw_by_distance=function(table,E_num){
  position=which(info[["Enum"]]==E_num)
  info$distance=0
  if (is.na(info[position,8])){
    return(NULL)
  }
  if (is.na(info[position,9])){
    return(NULL)
  }
  for (i in seq(1:nrow(info))){
    info[i,12]=distm(c(info[i,8], info[i,9]), c(info[position,8], info[position,9]), fun = distHaversine)
  }
  ordered_info=info[order(info$distance),][1:5,]
  draw_heatmap(pre_1000[,ordered_info$Enum],E_num)
  draw_PAV_heatmap(pre_1000[,ordered_info$Enum],E_num)
}
for (i in colnames(pre_1000)){
  print(i)
  draw_by_distance(table,i)
}
