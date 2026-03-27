library(ggplot2)
library(ggpubr)

################# Figure 5a - sv coverage ###################

sv <-read.csv("all.high.con.spe.tad.boundary.interior.sv.coverage.csv",head=T)
sv = subset(sv, region!='tad')
sv$conservation <- factor(sv$conservation, levels=c("H-conserved", "Conserved","Diverged"))

p1 = ggplot(sv,aes(x=conservation,y=sv.coverage))+
  geom_boxplot(aes(fill=region),outlier.shape=NA,width=0.6)+
  labs(x = "",y="SV coverage")+
  theme(strip.text = element_text(size = 15))+
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 9, face = "bold"),legend.position = 'none')+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(face = "bold",color = "black",size = 10,vjust = 0.5, hjust = 0.5,angle = 30))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#66c2a5","#80b1d3","#fc8d62"))

p_values_sv <- compare_means(sv.coverage ~ region, 
                               data = sv, 
                               group.by = "conservation",
                               method = "wilcox.test")
ggsave(p1,filename="all.high.con.spe.boundary.interior.sv.coverage.pdf",width=4.5,height=3.5)



################# Figure 5b - sv enrichment ###################
library(ggplot2)
library(gtrellis)
library(RColorBrewer)
library(circlize)
library("ComplexHeatmap")

df <- read.csv("all.spe.uptaddown.50kb.20bin.sv.breakpoint.observed.random.csv",head=T)
df.mat<-as.matrix(df)

labels_row = c("Ppse", "Pwua", "Psze","Plas","Pyun","Prot","Pqio","Pdav","Pade","Psim")
p2=pheatmap(df.mat,
            scale = "none",
            cluster_cols=F,
            cluster_rows=F,
            border_color=NA,
            show_rownames = T,
            labels_row =labels_row )
ggsave(p2,filename="all.diverged_tad.sv.enrichment.pdf",width=4.5,height=3.5)

################# Figure 5c - TE coverage ###################

dt= read.csv('all_spe_bound_SV.all_TE.coverage.csv')
dt$type=factor(d$type,levels = c("H-conserved", "Conserved","Diverged"))

p3 = ggplot(dt, aes(x=type,y=te.coverage,fill=sv))+
  geom_boxplot(outlier.colour = NA,width=0.6)+
  theme_bw()+
  theme(
    #panel.grid =element_blank(),
    #panel.background=element_blank()
    axis.text.y= element_text(size = 9, color = "black", face = "bold"),
    axis.title.y= element_text(size = 12, color = "black", face = "bold"),
    axis.title.x= element_text(size = 12, color = "black", face = "bold"),
    axis.text.x= element_text(face = "bold",color = "black",size = 10,vjust = 0.5, hjust = 0.5,angle = 30))+
  ylab('TE coverage')+xlab('')

p_values_te <- compare_means(te.coverage ~ sv, 
                          data = dt, 
                          group.by = "type",
                          method = "wilcox.test")
ggsave(p3, filename="all_spe_bound_SV.all_TE.coverage.pdf", width=4.5, height=3.5, units=c("cm"),colormodel="srgb")

