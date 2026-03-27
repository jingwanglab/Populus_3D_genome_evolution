#####ATAC peak in A/B 
library(ggplot2)

df= read.csv("all.ab.peak2gene_dis.csv",head=T) 

p=ggplot(df, aes(x=peak2gene_dis+0.001, color=compartment))+
  geom_density(size=1)+
  theme_bw() +
  theme(legend.title = element_blank())+
  scale_x_log10(breaks=c(0.001,1,100,1000,10000,100000,1000000),labels=c(0,1,100,1000,10000,100000,1000000)) +
  xlab('Distance between Peak and nearest gene (bp)')+ylab('Density')+
  scale_color_manual(values=c("#fee090","#91bfdb"))+
  theme(axis.text= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())
ggsave(p,filename="all.ab.peak2gene_dis.pdf",width=6.8,height=3.7)

