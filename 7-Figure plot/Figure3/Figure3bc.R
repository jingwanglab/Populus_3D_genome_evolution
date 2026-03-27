library(ggplot2)

########### Figure 3b - gene density ##########
gene <- read.csv("all.gene.density.10bin.csv")
p1 = ggplot(gene,aes(x=bin,y=gene.density))+
  geom_line(size=0.8,color="#fc8d62")+theme_bw()+
  labs(x = "",y="Gene density")+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold.italic"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,11,20,30),labels = c("-50k", "start","end","+50k"))
ggsave(p1,filename="all.gene.density.pdf",width=3.7,height=2.7)

########### Figure 3c - te density ##########
te <- read.csv("all.te.density.10bin.csv")
p2 = ggplot(te,aes(x=bin,y=te.density))+
  geom_line(size=0.8,color="#fc8d62")+theme_bw()+
  labs(x = "",y="TE density")+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold.italic"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,11,20,30),labels = c("-50k", "start","end","+50k"))
ggsave(p2,filename="all.te.density.pdf",width=3.7,height=2.7)
