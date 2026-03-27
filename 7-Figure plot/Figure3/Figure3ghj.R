library(ggplot2)

########### Figure 3g - methylation ##########
all.cpg.chg.chh.20bin <- read.csv("all.uptaddown.20bin.cpg.chg.chh.csv")
all.cpg.chg.chh.20bin$methylation <- factor(all.cpg.chg.chh.20bin$methylation, levels=c("cpg", "chg", "chh"))
p1=ggplot(all.cpg.chg.chh.20bin,aes(x=bin,y=value,colour=methylation))+
  geom_line(size=1)+theme_bw()+
  scale_colour_manual(values= c("#fc8d62","#1f78b4","#a6cee3"),
                      labels=c("CG", "CHG", "CHH"))+
  labs(x = "",y="Methylation levels (%)")+theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,21,40,60),labels = c("-50k", "start","end","+50k"))
ggsave(p1,filename="all.uptaddown.20bin.cpg.chg.chh.pdf",width=3.7,height=2.7)

########### Figure 3h - pangene density ##########
gene.density <- read.csv("all.pangene.density.csv")
p2 = ggplot(gene.density,aes(x=bin,y=gene.density,group=gene.class,colour=gene.class))+
  geom_line(size=0.8)+theme_bw()+
  #scale_colour_manual(values=c("#1f78b4","#b2df8a","#a6cee3"))+
  labs(x = "",y="Gene density")+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,11,20,30),labels = c("-50k", "start","end","+50k"))
ggsave(p2,filename="all.pangene.density.pdf",width=3.7,height=2.7)


########### Figure 3i - CNS enrichment ##########
library(ggplot2)
library(gtrellis)
library(RColorBrewer)
library(circlize)
library("ComplexHeatmap")

df=read.csv("tadup.tad.taddown.50kb.20bin.observed.random.cns.csv",head=T,row.names = 1)
df.mat<-as.matrix(df)
p3=pheatmap(df.mat,
            scale = "none",
            cluster_cols=F,
            cluster_rows=F,
            border_color=NA,
            show_rownames = T)
ggsave(p3,filename="CNS_enrichment.pdf",width=3.7,height=2.7)







