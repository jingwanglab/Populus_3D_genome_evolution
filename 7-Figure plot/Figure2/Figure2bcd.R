#################Figure 2b###################
library(ggplot2)
all.ab.change <- read.csv("all.ab.change.homologous.bin.csv",head=T)
all.ab.change$species <- factor(all.ab.change$species, levels=c("Ppse","Pwua","Pyun","Psze","Plas","Psim","Pade","Pdav","Prot","Pqio"), ordered=TRUE)
all.ab.change$status <- factor(all.ab.change$status, levels=c("A-B","B-A","A-A","B-B"))
p1=ggplot(all.ab.change, aes(x=species,y=number,fill=status))+
  geom_bar(position="fill", stat = "identity",width=0.6)+theme_bw()+
  theme(legend.title=element_blank(),legend.text =element_text(face = "bold",size = 11),legend.position = "top")+
  labs(x = " ",y="Proportion")+scale_fill_manual(values=c("#91bfdb","#fee090","#cccccc","#999999"))+
  theme(axis.text.x= element_text(size = 12, color = "black", face = "bold.italic", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y= element_text(size = 12, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 14, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  scale_y_continuous(labels = scales::percent)
ggsave(p1,filename="all.ab.change.homologous.bin.pdf",width=6.5,height=3.5)

#################Figure 2c###################
library(ggpubr)
library(ggpmisc)

df= read.csv("ab.change.ks.csv",head=T)

p2=ggplot(df, aes(x=ks, y = ab.conserved))+
  geom_point(color="#fee090",size=3)+theme_bw()+
  geom_smooth(method = "lm", color = "#91bfdb",size=1.2, fill = "lightgray")+
  stat_cor(method = "pearson", 
           label.x =0.04, label.y =0.3)+
  theme(axis.text= element_text(size = 12, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  labs(x = "Ks",y="The ratio of conserved compartment")
ggsave(p2,filename="ab.conserved.ks.cor.pdf",width=4.9,height=3.7)

#################Figure 2d###################
dt = read.csv("changed_AB_pan_stats.csv") 
p3 = ggplot(changeA,aes(x="",y=value, fill=pan_type))+
  geom_bar(stat='identity')+
  coord_polar(theta = 'y')+
  facet_wrap(~type)
ggsave(p3, filename="changed_AB_pan_stats.pdf", width=8, height=8, units=c("cm"),colormodel="srgb")
