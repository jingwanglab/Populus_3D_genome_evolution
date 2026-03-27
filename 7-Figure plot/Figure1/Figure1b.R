library(ggplot2)
library(ggpubr)

##### Gene expression in A and B compartments
ab.gene_exp <- read.csv("all.ab.gene.expression.csv")
p1=ggplot(ab.gene_exp,aes(x=compartment,y=log2(FPKM+1),fill=compartment))+
  geom_boxplot(outlier.shape=NA,width=0.5)+theme_bw()+
  labs(x = "",y="Gene expression")+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#fee090","#91bfdb"))+
  guides(fill="none")+
ggsave(p1,filename="all.ab.gene_exp.pdf",width=2.6,height=3.4)

##### TE expression in A and B compartments
ab.TE_exp <- read.csv("all.ab.TE.expression.csv")
p2=ggplot(ab.TE_exp,aes(x=compartment,y=log2(FPKM+1),fill=compartment))+
  geom_boxplot(outlier.shape=NA,width=0.5)+theme_bw()+
  labs(x = "",y="TE expression")+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#fee090","#91bfdb"))+
  guides(fill="none")+
ggsave(p2,filename="all.ab.TE_exp.pdf",width=2.6,height=3.4)
