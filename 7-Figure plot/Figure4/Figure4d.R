library(ggplot2)
library(ggpubr)

all.changed_tad.kaks <-read.csv("all.changed_tad.kaks.csv",head=T)
all.changed_tad.kaks$conservation <- factor(all.changed_tad.kaks$conservation, levels=c("H-conserved", "Conserved","Diverged"))

###########all.ka
p1 = ggplot(all.changed_tad.kaks,aes(x=conservation,y=ka))+
  geom_boxplot(aes(fill=conservation),outlier.shape=NA,width=0.6)+
  facet_wrap(~region,nrow=1)+theme_bw()+
  labs(x = "",y="ka")+
  theme(strip.text = element_text(size = 15))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 9, face = "bold"),legend.position = 'none')+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#66c2a5","#80b1d3","#fc8d62"))+
  stat_compare_means(method = 'wilcox.test', comparisons = list(c("H-conserved", "Conserved"),c("H-conserved","Diverged"),c("Conserved","Diverged")),label = 'p.format')

ka_p_values <- compare_means(ka ~ conservation, 
                          data = all.changed_tad.kaks, 
                          group.by = "region",
                          method = "wilcox.test")

###########all.ks
p2 = ggplot(all.changed_tad.kaks,aes(x=conservation,y=ks))+
  geom_boxplot(aes(fill=conservation),outlier.shape=NA,width=0.6)+
  facet_wrap(~region,nrow=1)+theme_bw()+
  labs(x = "",y="ks")+
  theme(strip.text = element_text(size = 15))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 9, face = "bold"),legend.position = 'none')+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#66c2a5","#80b1d3","#fc8d62"))+
  stat_compare_means(method = 'wilcox.test', comparisons = list(c("H-conserved", "Conserved"),c("H-conserved","Diverged"),c("Conserved","Diverged")),label = 'p.format')

ks_p_values <- compare_means(ks ~ conservation, 
                             data = all.changed_tad.kaks, 
                             group.by = "region",
                             method = "wilcox.test")

###########all.ka/ks
all.changed_tad.kaks$kaks = all.changed_tad.kaks$ka/all.changed_tad.kaks$ks
p3 = ggplot(all.changed_tad.kaks,aes(x=conservation,y=kaks))+
  geom_boxplot(aes(fill=conservation),outlier.shape=NA,width=0.6)+
  facet_wrap(~region,nrow=1)+theme_bw()+
  labs(x = "",y="ka/ks")+
  theme(strip.text = element_text(size = 15))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 9, face = "bold"),legend.position = 'none')+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#66c2a5","#80b1d3","#fc8d62"))+
  stat_compare_means(method = 'wilcox.test', comparisons = list(c("H-conserved", "Conserved"),c("H-conserved","Diverged"),c("Conserved","Diverged")),label = 'p.format')

kaks_p_values <- compare_means(kaks ~ conservation, 
                             data = all.changed_tad.kaks, 
                             group.by = "region",
                             method = "wilcox.test")


ggsave(p1,filename="all.changed_tad.ka.pdf",width=8,height=3.7)
ggsave(p2,filename="all.changed_tad.ks.pdf",width=8,height=3.7)
ggsave(p3,filename="all.changed_tad.kaks.pdf",width=8,height=3.7)
