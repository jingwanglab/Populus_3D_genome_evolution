
library(readxl)
library(ggplot2)
library(cowplot)

##### The proportion of A/B compartment across the genomes of the eleven Populus species
ab.ratio <- read_excel("all_species_AB_proportion.xlsx")
ab.ratio$species <- factor(ab.ratio$species, levels=c("Ppse", "Pwua", "Psze","Pyun","Pkor","Psim","Plas","Pdav","Prot","Pqio","Pade"), ordered=TRUE)
p1=ggplot(ab.ratio,aes(x=species,y=ratio,fill=compartment))+
  geom_bar(stat="identity",position="stack",width=0.7)+theme_bw()+
  labs(x = "",y="Proportion of compartment(%)")+theme(legend.position = "top") +
  theme(legend.title=element_blank(),legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black",face = "bold.italic",size = 11,vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_fill_manual(values=c("#fee090","#91bfdb"),labels=c("A compartment","B compartment"))
ggsave(p1,filename="all.ab.proportion.pdf",width=8,height=4.4)

###### The proportions of genes identified in A and B compartments
ab.gene.ratio <- read_excel("all_species_AB_genes_stats.xlsx")
p2 = ggplot(ab.gene.ratio, aes(x='',value,fill=variable))+
  geom_bar(stat='identity',position = 'stack')+
  coord_polar(theta = 'y')+
  labs(x = "", y = "", title = "") +
  theme(axis.ticks = element_blank(),
        legend.position = 'top') +
  facet_wrap(~num, ncol=1) +
  scale_fill_manual(values=c("#fee090","#91bfdb"))+
  theme_nothing()
ggsave(p2, filename="ab.gene.ratio.pdf", width=0.6, height=8)

###### The proportions of TEs identified in A and B compartments
ab.te.ratio <- read_excel("all_species_AB_TE_stats.xlsx")
p3 = ggplot(ab.te.ratio, aes(x='',value,fill=variable))+
  geom_bar(stat='identity',position = 'stack')+
  coord_polar(theta = 'y')+
  labs(x = "", y = "", title = "") +
  theme(axis.ticks = element_blank(),
        legend.position = 'top') +
  facet_wrap(~num, ncol=1) +
  scale_fill_manual(values=c("#fee090","#91bfdb"))+
  theme_nothing()
ggsave(p3, filename="ab.te.ratio.pdf", width=0.6, height=8)
