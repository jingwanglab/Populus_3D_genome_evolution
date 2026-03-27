##### Pangene in A/B
library(ggplot2)
library(cowplot)
pan=c('core','core','dis','dis','pri','pri')
cmp = c('A','B','A','B','A','B')
ratio = c(0.684,0.316,0.604,0.396,0.568,0.432)
df =data.frame(pan,cmp,ratio)
a1 = ggplot(df, aes(x=pan,ratio,fill=cmp))+
  geom_bar(width=0.5,stat='identity',position = 'stack')+
  coord_polar(theta = 'y')+
  labs(x = "", y = "", title = "") +
  theme(axis.ticks = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values=c("#41B3E2", "#F1786B"))+
  theme_nothing()
ggsave(a1, filename="pan.AB.ratio.pdf", width=6, height=6, units=c("cm"),colormodel="srgb")

##### CNS in A/B
dt = read.csv('ab.cns.csv',sep='\t',header=FALSE)
a2 = ggplot(dt, aes(AB,count,fill=AB))+
  geom_violin(width=0.9)+
  geom_boxplot(width=0.2,outlier.colour = NA)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none')+
  scale_fill_manual(values=c("#41B3E2", "#F1786B"))+
  xlab('')+ylab('CNS count')

ggsave(a2, filename="CNS.AB.count.pdf", width=4, height=4, units=c("cm"),colormodel="srgb")

