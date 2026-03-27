library(topGO)
library(dplyr)
library(data.table)
library(GO.db)  

d <- read.table(file="all_spe_pep_removeStop.fa.GO", header=FALSE)
colnames(d) <- c("gene_id","go_id")
d$gene_id <- as.character(d$gene_id)
d$go_id   <- as.character(d$go_id)

# gene2GO: GO list corresponding to each gene
all_go <- lapply(split(d, sub("\\.\\d+$", "", d[, 1])), function(x) unique(x[, 2]))
geneNames <- names(all_go)

gene <- read.csv("all.high.con.tad.geneid.csv", header=FALSE)
colnames(gene) <- c("geneid")
outlier_gene <- as.vector(gene$geneid)

geneList <- factor(as.integer(geneNames %in% outlier_gene))
names(geneList) <- geneNames

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = all_go)

restRes <- runTest(GOdata, algorithm="weight01", statistic="fisher")

# 1) All GO terms used
go_ids <- usedGO(GOdata)

# 2) Get raw p-values for each term (full set)
pv_all <- score(restRes, whichGO = go_ids)
pv_all <- pv_all[!is.na(pv_all)]

# 3) Statistics for Annotated / Significant / Expected (full set)
ts_all <- termStat(GOdata, whichGO = names(pv_all))

# 4) BH correction for all terms (all terms corrected together)
padj_all <- p.adjust(pv_all, method="BH")

# 5) Create results table
gene_table_all <- data.frame(
  GO.ID       = names(pv_all),
  Term        = Term(GOTERM[names(pv_all)]),
  Annotated   = ts_all$Annotated,
  Significant = ts_all$Significant,
  Expected    = ts_all$Expected,
  P.Value     = unname(pv_all),
  FDR         = unname(padj_all),
  RichFactor  = ts_all$Significant / ts_all$Annotated,
  stringsAsFactors = FALSE
)

gene_table_all <- gene_table_all[order(gene_table_all$FDR, gene_table_all$P.Value), ]
write.csv(gene_table_all,
          file = "all.high.con.goenrichment.csv",
          row.names = FALSE,
          quote = FALSE)


##### Plot
#install.packages("forcats")
library(forcats)
library(ggplot2)

d1 = read.csv('all.high.con.goenrichment.csv')
d2 = read.csv('all.con.goenrichment.csv')
d3 = read.csv('all.spe.tad.goenrichment.csv')

d1$type='H-conserved'
d2$type='Conserved'
d3$type='Diverged'
d1_10 <- d1[1:10, ]
d2_10 <- d2[1:10, ]
d3_10 <- d3[1:10, ]
#d1_sorted <- d1_10[order(d1_10$P.Value, decreasing = TRUE), ]
#d2_sorted <- d2_10[order(d2_10$P.Value, decreasing = TRUE), ]
#d3_sorted <- d3_10[order(d3_10$P.Value, decreasing = TRUE), ]

dt = rbind(d1,d2,d3)
dt=rbind(d1_10,d2_10,d3_10)

#dt=subset(dt,FDR<0.05)
#dt=subset(dt,P.Value<0.05)
#write.csv(dt, 'all.diverged_TAD.GO.csv',row.names = FALSE)
dt$conservation = factor(dt$conservation, levels=c('H-conserved','Conserved','Diverged'))

p = ggplot(dt,aes(conservation,Term))+
  geom_point(aes(size=RichFactor,color=-log10(as.numeric(FDR))))+
  scale_y_discrete(limits=dt$Term)  +
  #scale_color_distiller(palette = "Spectral")+
  scale_color_gradientn(colours = c("#2166AC", "#A6CEE3", "white", "orange", "red"))+ 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 7 , color='black',face="bold",angle=90),
        axis.text.y = element_text(size = 7 , color='black',face="bold"),
        panel.grid.major = element_line(size=0.1,colour = "#F4F2F3"),
        panel.grid.minor = element_line(size=0.1,colour = "#F4F2F3"),
        axis.title.x  = element_text( size= 9,color='black',face="bold" ))+
  labs(y="",x="",size="Rich factor")

ggsave(p,filename = 'all.TAD.Go_plot.pdf',width=20,height=16,units=c("cm"),dpi=300)

