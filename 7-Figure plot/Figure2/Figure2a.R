library(ggplot2)
library(gtrellis)
library(RColorBrewer)
library(circlize)
library(pheatmap)

df.barplot <- read.csv("xiangy.barplot.csv",head=T)
df.A.barplot <- read.csv("xiangy.A.barplot.csv",head=T)
df.B.barplot <- read.csv("xiangy.B.barplot.csv",head=T)

### all stats barplot
p1 = ggplot(df.barplot, aes(x=species,y=number,fill=type))+
  geom_bar(position="fill", stat = "identity",width=0.6)+theme_bw()+
  theme(legend.title=element_blank(),legend.text =element_text(face = "bold",size = 11),legend.position = "top")+
  theme(axis.text.x= element_text(size = 12, color = "black", face = "bold.italic", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 14, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_y_continuous(labels = scales::percent)
ggsave(p1,filename="1.all_stats.pdf",width=2,height=5)

###### changed A heatmap
df.changed.A <- read.csv("changed.homologous.bin.A.species_sort_ks.matrix.csv",head=T,row.names = 1)
dfchangedall.A<-as.matrix(df.changed.A)
p2 = pheatmap(dfchangedall.A,
            scale = "none",
            cluster_cols=F,
            cluster_rows=F,
            border_color=NA,
            show_rownames = F,
            color = colorRampPalette(colors = c("#67A9CF","white","#FCD383"))(100))
ggsave(p2,filename="2.A.heatmap.pdf",width=4,height=5)

###### changed B heatmap
df.changed.B <- read.csv("changed.homologous.bin.B.species_sort_ks.matrix.csv",head=T,row.names = 1)
dfchangedall.B<-as.matrix(df.changed.B)
p3 = pheatmap(dfchangedall.B,
           scale = "none",
           cluster_cols=F,
           cluster_rows=F,
           border_color=NA,
           show_rownames = F,
           color = colorRampPalette(colors = c("#67A9CF","white","#FCD383"))(100))
ggsave(p3,filename="3.B.heatmap.pdf",width=4,height=5)

###### changed A stats barplot
df.A.barplot$pan <- factor(df.A.barplot$pan, levels=c("1","2","3","4","5","6","7","8"), ordered=TRUE)
p4  = ggplot(df.A.barplot, aes(x=pan,y=number,fill=pan))+
  geom_bar(stat = "identity",width=0.6)+theme_bw()+
  theme(legend.title=element_blank(),legend.text =element_text(face = "bold",size = 11),legend.position = "top")+
  theme(axis.text.x= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 14, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())
ggsave(p4,filename="4.A.stats.pdf",width=4,height=5)

###### changed B stats barplot
df.B.barplot$pan <- factor(df.B.barplot$pan, levels=c("1","2","3","4","5","6","7","8","9","10"), ordered=TRUE)
p5 = ggplot(df.B.barplot, aes(x=pan,y=number,fill=pan))+
  geom_bar(stat = "identity",width=0.6)+theme_bw()+
  theme(legend.title=element_blank(),legend.text =element_text(face = "bold",size = 11),legend.position = "top")+
  theme(axis.text.x= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 14, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())
ggsave(p5,filename="5.B.stats.pdf",width=4,height=5)





