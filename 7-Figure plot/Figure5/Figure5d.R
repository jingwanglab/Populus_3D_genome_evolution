library(ggplot2)
library(ggsignif)

####inv
dt1 = read.csv('Pkor.spe.tad.inv.overlap',sep='\t',header = FALSE)
dt1$len = dt1$V6-dt1$V5
dt2 = read.csv('Pkor.hcon.con.tad.inv.overlap',sep='\t',header = FALSE)
dt2$len = dt2$V6-dt2$V5
dt1$type='diverged'
dt2$type='conserved'
dt = rbind(dt1,dt2)

p1=ggplot(dt, aes(log10(len),fill=type))+
  geom_density(alpha=0.3)+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.background=element_blank())

p2=ggplot(dt, aes(type,log10(len),fill=type))+
  geom_boxplot(width=0.5,outlier.colour = NA)+
  geom_signif(comparisons = list(c("diverged", "conserved")))+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.backgrSound=element_blank())
ggsave(p1, filename="Pkor_tad_INV.overlap1.pdf", width=8, height=5, units=c("cm"),colormodel="srgb")
ggsave(p2, filename="Pkor_tad_INV.overlap2.pdf", width=8, height=5, units=c("cm"),colormodel="srgb")

####indel
dt1 = read.csv('Pkor.spe.tad.indel.overlap',sep='\t',header = FALSE)
dt1$len = dt1$V6-dt1$V5
dt2 = read.csv('Pkor.hcon.con.tad.indel.overlap',sep='\t',header = FALSE)
dt2$len = dt2$V6-dt2$V5
dt1$type='diverged'
dt2$type='conserved'
dt = rbind(dt1,dt2)

p1=ggplot(dt, aes(log10(len),fill=type))+
  geom_density(alpha=0.3)+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.background=element_blank())

p2=ggplot(dt, aes(type,log10(len),fill=type))+
  geom_boxplot(width=0.5,outlier.colour = NA)+
  geom_signif(comparisons = list(c("diverged", "conserved")))+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.background=element_blank())

ggsave(p1, filename="Pkor_tad_indel.overlap1.pdf", width=8, height=5, units=c("cm"),colormodel="srgb")
ggsave(p2, filename="Pkor_tad_indel.overlap2.pdf", width=8, height=5, units=c("cm"),colormodel="srgb")


####dup
dt1 = read.csv('Pkor.spe.tad.dup.overlap',sep='\t',header = FALSE)
dt1$len = dt1$V6-dt1$V5
dt2 = read.csv('Pkor.hcon.con.tad.dup.overlap',sep='\t',header = FALSE)
dt2$len = dt2$V6-dt2$V5
dt1$type='diverged'
dt2$type='conserved'
dt = rbind(dt1,dt2)

p1=ggplot(dt, aes(log10(len),fill=type))+
  geom_density(alpha=0.3)+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.background=element_blank())

p2=ggplot(dt, aes(type,log10(len),fill=type))+
  geom_boxplot(width=0.5,outlier.colour = NA)+
  geom_signif(comparisons = list(c("diverged", "conserved")))+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.background=element_blank())

ggsave(p1, filename="Pkor_tad_dup.overlap1.pdf", width=8, height=5, units=c("cm"),colormodel="srgb")
ggsave(p2, filename="Pkor_tad_dup.overlap2.pdf", width=8, height=5, units=c("cm"),colormodel="srgb")

####TRANS
dt1 = read.csv('Pkor.spe.tad.trans.overlap',sep='\t',header = FALSE)
dt2 = read.csv('Pkor.hcon.con.tad.trans.overlap',sep='\t',header = FALSE)
dt1$len = dt1$V6-dt1$V5
dt2$len = dt2$V6-dt2$V5
dt1$type='diverged'
dt2$type='conserved'
dt = rbind(dt1,dt2)

p1 = ggplot(dt, aes(log10(len),fill=type))+
  geom_density(alpha=0.3)+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.background=element_blank())

p2=ggplot(dt, aes(type,log10(len),fill=type))+
  geom_boxplot(width=0.5,outlier.colour = NA)+
  geom_signif(comparisons = list(c("diverged", "conserved")))+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.background=element_blank())

ggsave(p1, filename="Pkor_con_tad_trans.overlap1.pdf", width=8, height=5, units=c("cm"),colormodel="srgb")
ggsave(p2, filename="Pkor_tad_trans.overlap2.pdf", width=8, height=5, units=c("cm"),colormodel="srgb")

