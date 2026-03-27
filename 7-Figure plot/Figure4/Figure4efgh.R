library(ggplot2)
library(ggpubr)
library(cowplot)

################# Figure 4e - CNS enrichment ###################

df=read.csv("all.changed_tad.flank50kb.20bin.cns.enrichment.csv",head=T)
df$conservation<- factor(df$conservation, levels=c("H-conserved", "Conserved","Diverged"))
p1 =ggplot(df,aes(x=bin,y=o.r,color=conservation))+
  geom_line(size=0.8)+theme_bw()+
  labs(x = " ",y="CNS enrichment")+
  scale_colour_manual(values= c("#a6cee3","#67a9cf","#fc8d62"))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 8),legend.position = "top")+
  theme(axis.text.y= element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold",size = 8,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 10, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,21,40,60),labels = c("-50k", "start","end","+50k"))
ggsave(p1,filename="all.changed_tad.flank50kb.20bin.cns.enrichment.pdf",width=5.7,height=4.0)


################# Figure 4f - gene expression ###################

dt <-read.csv("all.changed_tad.gene_expression.csv",head=T)
dt$conservation <- factor(dt$conservation, levels=c("H-conserved", "Conserved","Diverged"))

p2 = ggplot(dt,aes(x=conservation,y=log2(TPM+1)))+
  geom_boxplot(aes(fill=conservation),outlier.shape=NA,width=0.6)+
  facet_wrap(~region,nrow=1)+theme_bw()+
  labs(x = "",y="Gene expression")+
  theme(strip.text = element_text(size = 15))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 9, face = "bold"),legend.position = 'none')+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#66c2a5","#80b1d3","#fc8d62"))+
  stat_compare_means(method = 'wilcox.test', comparisons = list(c("H-conserved", "Conserved"),c("H-conserved","Diverged"),c("Conserved","Diverged")),label = 'p.format')

p_values_tpm <- compare_means(TPM ~ conservation, 
                              data = dt, 
                              group.by = "region",
                              method = "wilcox.test")
ggsave(p2,filename="all.changed_tad.gene_expression.pdf",width=10,height=4.0)


################# Figure 4h - methylation ###################

all.chg <-read.csv("all.high.con.spe.uptaddown.20bin.chg.csv",head=T)
all.chh <-read.csv("all.high.con.spe.uptaddown.20bin.chh.csv",head=T)
all.cpg <-read.csv("all.high.con.spe.uptaddown.20bin.cpg.csv",head=T)

## chg
all.chg$conservation <- factor(all.chg$conservation, levels=c("H-conserved", "Conserved","Diverged"))
c1 = ggplot(all.chg,aes(x=bin,y=methylation,colour=conservation))+
  geom_line(size=0.8)+theme_bw()+
  scale_colour_manual(values= c("#66c2a5","#80b1d3","#fc8d62"))+
  labs(x = "",y="")+ggtitle("CHG")+theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.text= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,21,40,60),labels = c("-50k", "start","end","+50k"))

## chh
all.chh$conservation <- factor(all.chh$conservation, levels=c("H-conserved", "Conserved","Diverged"))
c2 = ggplot(all.chh,aes(x=bin,y=methylation,colour=conservation))+
  geom_line(size=0.8)+theme_bw()+
  scale_colour_manual(values= c("#66c2a5","#80b1d3","#fc8d62"))+
  labs(x = "",y="")+ggtitle("CHH")+theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.text= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,21,40,60),labels = c("-50k", "start","end","+50k"))
## cpg
all.cpg$conservation <- factor(all.cpg$conservation, levels=c("H-conserved", "Conserved","Diverged"))
c3 = ggplot(all.cpg,aes(x=bin,y=methylation,colour=conservation))+
  geom_line(size=0.8)+theme_bw()+
  scale_colour_manual(values= c("#66c2a5","#80b1d3","#fc8d62"))+
  labs(x = "",y="")+ggtitle("CG")+theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.text= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_x_discrete(limits = c(1,21,40,60),labels = c("-50k", "start","end","+50k"))

c =plot_grid(c3,c1,c2,nrow = 1)
ggsave(c,filename="all.changed_tad.flank50kb.20bin.methy.pdf",width=10,height=4.0)


#################### Figure 4g - ATAC peak ###################

df.all= read.csv("all.changed_tad.atac.csv",head=T) 
df.all$conservation <- factor(df.all$conservation, levels=c("H-conserved", "Conserved","Diverged"))

p3 =ggplot(df.all,aes(x=conservation,y=log10(peak2gene_dis+1)))+
  geom_boxplot(aes(fill=conservation),outlier.shape=NA,width=0.5)+
  facet_wrap(~region,nrow=1)+theme_bw()+
  labs(x = " ",y="Distance between peak and nearest gene")+
  theme(strip.text = element_text(size = 10))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 8, face = "bold"),legend.position = 'none')+
  theme(axis.text.y= element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 8,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 10, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#a6cee3","#67a9cf","#fc8d62"))+
  stat_compare_means(method = 'wilcox.test', comparisons = list(c("H-conserved", "Conserved"),c("H-conserved","Diverged"),c("Conserved","Diverged")),label = 'p.signif')

p_values_atac <- compare_means(peak2gene_dis ~ conservation, 
                             data = df.all, 
                             group.by = "region",
                             method = "wilcox.test")
ggsave(p3,filename="all.changed_tad.atac.pdf",width=10,height=4.4)

