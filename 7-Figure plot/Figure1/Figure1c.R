
#### Comparison of methylation levels (CG, CHG and CHH) of A/B compartments across species
library(ggplot2)
library(ggpubr)

all.ab.chg <-read.csv("all.ab.chg.csv",head=T)
all.ab.chh <-read.csv("all.ab.chh.csv",head=T)
all.ab.cpg <-read.csv("all.ab.cpg.csv",head=T)

############chg
p1=ggplot(all.ab.chg,aes(x=status,y=methylation,fill=status))+
  geom_boxplot(lwd=0.6,outlier.shape = NA,width=0.6)+theme_bw()+
  theme(legend.position="none")+
  labs(x = " ",y="CHG")+
  scale_fill_manual(values=c("#fee090","#91bfdb"))+
  theme(axis.text= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x= element_text(size = 15, color = "black", face = "italic", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())
ggsave(p1,filename="all.ab.chg.pdf",width=2.6,height=3)

############chh
p2=ggplot(all.ab.chh,aes(x=status,y=methylation/100,fill=status))+
  geom_boxplot(lwd=0.6,outlier.shape = NA,width=0.6)+theme_bw()+
  theme(legend.position="none")+
  labs(x = " ",y="CHH")+
  scale_fill_manual(values=c("#fee090","#91bfdb"))+
  theme(axis.text= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x= element_text(size = 15, color = "black", face = "italic", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())
ggsave(p2,filename="all.ab.chh.pdf",width=2.6,height=3)

###############cpg
p3=ggplot(all.ab.cpg,aes(x=status,y=methylation,fill=status))+
  geom_boxplot(lwd=0.6,outlier.shape = NA,width=0.6)+theme_bw()+
  theme(legend.position="none")+
  labs(x = " ",y="CG")+
  scale_fill_manual(values=c("#fee090","#91bfdb"))+
  theme(axis.text= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x= element_text(size = 15, color = "black", face = "italic", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())
ggsave(p3,filename="all.ab.cpg.pdf",width=2.6,height=3)


