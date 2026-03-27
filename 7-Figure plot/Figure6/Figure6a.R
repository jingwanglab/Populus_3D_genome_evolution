library(ggplot2)
library(gridExtra)

dup = read.csv('./Pkor/bound_sv/dup_stats_count.txt',header=F)
p1=ggplot(dup,aes(V1))+
  geom_density(fill='#a6cee3',alpha=0.9)+
  geom_vline(xintercept = 1693)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab('Density')+xlab('DUP count')

indel = read.csv('./Pkor/bound_sv/indel_stats_count.txt',header=F)
p2=ggplot(indel,aes(V1))+
  geom_density(fill='#2166ac',alpha=0.8)+
  geom_vline(xintercept = 4637)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab('Density')+xlab('InDel count')

inv = read.csv('./Pkor/bound_sv/inv_stats_count.txt',header=F)
p3=ggplot(inv,aes(V1))+
  geom_density(fill='#fc8d62',alpha=0.8)+
  geom_vline(xintercept = 364)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab('Density')+xlab('INV count')

tra = read.csv('./Pkor/bound_sv/trans_stats_count.txt',header=F)
p4=ggplot(tra,aes(V1))+
  geom_density(fill='#d73027',alpha=0.7)+
  geom_vline(xintercept = 1126)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab('Density')+xlab('TRANS count')

p = grid.arrange(p2,p1,p3,p4,nrow=1)
ggsave(p, filename="./Pkor/bound_sv/SV_random_count.pdf", width=15, height=5, units=c("cm"),colormodel="srgb")


