library(ggplot2)


tadsize=read.csv("tadsize.csv",header = F) 
tadsize$V2=factor(tadsize$V2)

p=ggplot(tadsize,aes(x=V1/1000,group=V2,colour=V2,fill=V2))+theme_bw()+theme(legend.position=c(0.85,0.53))+
  geom_line(stat ="density",adjust=1,size=1)+scale_x_continuous(limits=c(0,1500),breaks=seq(0,1500,250))+
  scale_colour_manual(values= c("#6a51a3","#807dba","#9e9ac8","#bcbddc","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c"),
                      labels=c("P.pseudoglauca", "P.wuana", "P.szechuanica","P.lasiocarpa","P.yunnanensis","P.rotundifolia","P.qiongdaoensis","P.davidiana","P.koreana","P.adenopoda","P.simonii"))+
  theme(legend.title=element_blank(),legend.text =element_text(face = "italic",size = 11))+
  theme(axis.text = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+ylab("Density")+
  theme(axis.title= element_text(size = 14, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+xlab("TAD size(kb)")+
  theme(panel.grid =element_blank(),panel.background = element_blank())
ggsave(p,filename="tadsize.pdf",width=8,height=5)
