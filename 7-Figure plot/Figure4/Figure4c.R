library(ggplot2)

#################Figure 4c###################

df= read.csv("tad.change.ks.csv",head=T)
p1=ggplot(df, aes(x=ks, y = tad.conserved))+
  geom_point(color="#fc8d62",size=3)+theme_bw()+
  geom_smooth(method = "lm", color = "#91bfdb", size=1.2, fill = "lightgray")+
  stat_cor(method = "pearson", 
           label.x =0.062, label.y =0.43)+
  theme(axis.text= element_text(size = 12, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.title= element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid =element_blank(),panel.background=element_blank())+
  labs(x = "Ks",y="The ratio of conserved TAD")
ggsave(p1,filename="tad.conserved.ks.cor.pdf",width=4.9,height=3.7)
