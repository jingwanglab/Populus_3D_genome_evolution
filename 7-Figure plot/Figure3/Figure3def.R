library(ggplot2)
library(ggpubr)

########### Figure 3d - gene expression ##########
gene <- read.csv("all.tad.boundary.20kb.tpm.csv")
p1=ggplot(gene,aes(x=region,y=log2(TPM+1),fill=region))+
  geom_boxplot(outlier.shape=NA,width=0.5)+theme_bw()+
  labs(x = "",y="Gene expression")+
  scale_fill_manual(values = c("#a1d99b", "#fdae6b"))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 12,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  guides(fill="none")+
  scale_x_discrete(labels=c("boundary","interior"))+
#  stat_compare_means(method = 'wilcox.test', comparisons = list(c("boundary", "tad")),label = 'p.signif',label.y = c(9.5))
ggsave(p1,filename="all.tad.boundary.20kb.tpm.pdf",width=2.6,height=3.4)

########### Figure 3e - TE expression ##########
te <- read.csv("all.tad.boundary.20kb.te.tpm.csv")
p2=ggplot(te,aes(x=region,y=log2(TPM+1),fill=region))+
  geom_boxplot(outlier.shape=NA,width=0.5)+theme_bw()+
  labs(x = "",y="TE expression")+
  scale_fill_manual(values = c("#a1d99b", "#fdae6b"))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 12,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  guides(fill="none")+
  scale_x_discrete(labels=c("boundary","interior"))+
#  stat_compare_means(method = 'wilcox.test', comparisons = list(c("boundary", "tad")),label = 'p.signif',label.y = c(9.5))
ggsave(p2,filename="all.tad.boundary.20kb.te.tpm.pdf",width=2.6,height=3.4)

########### Figure 3f - ATAC peak ##########
df= read.csv("all.boundary.interior.peak2gene_dis.csv",head=T) 
df <- subset(df,region=="boundary" | region=="interior")
p3=ggplot(df, aes(x=region,y=log10(peak2gene_dis+1),fill=region))+
  geom_boxplot(lwd=0.6,outlier.shape=NA,width=0.5)+theme_bw()+
  theme(legend.position="none")+
  labs(x = " ",y="Distance between peak and nearest gene")+
  scale_fill_manual(values=c("#67a9cf","#fc8d62"))+
  theme(axis.text= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
#  stat_compare_means(method = 'wilcox.test', comparisons = list(c('boundary','interior')),label = 'p.signif',size=4)
ggsave(p3,filename="all.boundary.interior.peak2gene_dis.pdf",width=2.3,height=4.2)
