library(ggplot2)
library(ggpubr)
#################### Figure2e - TE coverage #############
te.coverage<-read.csv("all.changed.conserved.homologous.bin.te.coverage.csv",head=T)
te.coverage$status <- factor(te.coverage$status, levels=c("Conserved","Changed"))

p1=ggplot(te.coverage,aes(x=status,y=te.coverage,fill=status))+
  geom_boxplot(outlier.shape=NA,width=0.5)+theme_bw()+
  labs(x = " ",y="TE coverage")+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 8, face = "bold"),legend.position = 'none')+
  theme(axis.text.y= element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 8,vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x= element_text(size = 10, color = "black", face = "bold.italic", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y= element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#91bfdb","#fee090"))+
  stat_compare_means(method = 'wilcox.test', comparisons = list(c("Conserved","Changed")))+
  scale_y_continuous(limits = c(0,1.1),breaks=seq(0, 1.1, 0.5))
ggsave(p1,filename="all.changed.conserved.homologous.bin.te.coverage.pdf",width=3.2,height=4.2)

#################### Figure2f - SV coverage #############
sv.coverage<-read.csv("all.changed.conserved.homologous.bin.sv.coverage.csv",head=T)
sv.coverage$status <- factor(sv.coverage$status, levels=c("Conserved","Changed"))

p2=ggplot(sv.coverage,aes(x=status,y=sv.coverage,fill=status))+
  geom_boxplot(outlier.shape=NA,width=0.5)+theme_bw()+
  labs(x = " ",y="SV coverage")+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 8, face = "bold"),legend.position = 'none')+
  theme(axis.text.y= element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 8,vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x= element_text(size = 10, color = "black", face = "bold.italic", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y= element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#91bfdb","#fee090"))+
  stat_compare_means(method = 'wilcox.test', comparisons = list(c("Conserved","Changed")))+
  scale_y_continuous(limits = c(0,0.5),breaks=seq(0, 0.5, 0.1))
ggsave(p2,filename="all.changed.conserved.homologous.bin.sv.coverage.pdf",width=2,height=4.8)


#################### Figure2g - gene expression #############

df <-read.csv("all.ab.ba.gene.alltra_tpm.csv",head=T)
df$species1 <- factor(df$species1, levels=c("Pkor", "Others"))
df$status <- factor(df$status, levels=c("A-B","B-A"))

p3=ggplot(df,aes(x=status,y=log2(TPM+1),fill=species1))+
  geom_boxplot(outlier.shape=NA,outlier.size=0.8,width=0.5)+theme_bw()+
  labs(x = "",y="Gene expression")+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold.italic"),legend.position = "top")+
  theme(axis.text.y= element_text(size = 11, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 12, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#67A9CF","#FCD383"),labels=c("P. kor","Others"))+
  stat_compare_means(method = 'wilcox.test', label.y = c(10,10))+
  guides(fill=guide_legend(keywidth = 0.8,keyheight = 0.8))
ggsave(p3,filename="all.ab.ba.gene.alltra_tpm.pdf",width=3,height=4.5)


#################### Figure2h - methylation #############

chg <-read.csv("all.Pkor.chg.csv",head=T)
chh <-read.csv("all.Pkor.chh.csv",head=T)
cpg <-read.csv("all.Pkor.cpg.csv",head=T)

##chg
chg$status <- factor(chg$status, levels=c("A-A","A-B","B-B","B-A"))
c1=ggplot(chg,aes(x=bin,y=methylation,colour=status))+
  geom_line(size=0.8)+theme_bw()+
  scale_colour_manual(values= c("#cccccc","#91bfdb","#999999","#fee090"),labels=c("A-A","A-B","B-B","B-A"))+
  labs(x = "",y="")+ggtitle("CHG")+theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,21,50,70),labels = c("-2k", "start","end","+2k"))
ggsave(c1,filename="chg.pdf",width=4,height=2.4)
##chh
chh$status <- factor(chh$status, levels=c("A-A","A-B","B-B","B-A"))
c2=ggplot(chh,aes(x=bin,y=methylation,colour=status))+
  geom_line(size=0.8)+theme_bw()+
  scale_colour_manual(values= c("#cccccc","#91bfdb","#999999","#fee090"),labels=c("A-A","A-B","B-B","B-A"))+
  labs(x = "",y="")+ggtitle("CHH")+theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,21,50,70),labels = c("-2k", "start","end","+2k"))
ggsave(c2,filename="all.chh.pdf",width=4,height=2.4)
##cpg
cpg$status <- factor(cpg$status, levels=c("A-A","A-B","B-B","B-A"))
c3=ggplot(cpg,aes(x=bin,y=methylation,colour=status))+
  geom_line(size=0.8)+theme_bw()+
  scale_colour_manual(values= c("#cccccc","#91bfdb","#999999","#fee090"),labels=c("A-A","A-B","B-B","B-A"))+
  labs(x = "",y="")+ggtitle("CG")+theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,21,50,70),labels = c("-2k", "start","end","+2k"))
ggsave(c3,filename="all.cpg.pdf",width=4,height=2.4)



