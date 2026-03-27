library(DESeq2)
library(tidyr)

### Plas
coldata2 <- read.csv("/Plas/coldata.csv",header = T)
combat_count <- read.csv("/Plas/gene_count_matrix.csv",row.names = 1,header = T)
dds <- DESeqDataSetFromMatrix(countData = combat_count, 
                              colData = coldata2, 
                              design = ~ condition) 

dds <- DESeq(dds)
res = results(dds, contrast=c("condition", "hs","control"))
res = res[order(res$pvalue),]
summary(res)
write.csv(res, file = "/Plas/ck_hs.deseq.csv")

### Prot
coldata2 <- read.csv("/Prot/coldata.csv",header = T)
combat_count <- read.csv("/Prot/gene_count_matrix.csv",row.names = 1,header = T)
dds <- DESeqDataSetFromMatrix(countData = combat_count, 
                              colData = coldata2, 
                              design = ~ condition) 
dds <- DESeq(dds)
res = results(dds, contrast=c("condition", "hs","control"))
res = res[order(res$pvalue),]
summary(res)
write.csv(res, file = "/Prot/ck_hs.deseq.csv")


### Plot

dt = read.csv('/Plas/ck_hs.deseq.csv')
dt2 = dt[complete.cases(dt),]
dt2$significant <- "unchanged"
dt2$significant[dt2$padj < 0.05 & dt2$log2FoldChange >= 1 ] <- "up"
dt2$significant[dt2$padj < 0.05 & dt2$log2FoldChange <= -1 ] <- "down"
sig_genes <- subset(dt2,gene=="Polas06125")  ##Prot16230

p = ggplot(dt2,aes(log2FoldChange, -log10(padj),color = significant)) +
  geom_point(size=0.8) +
  geom_point(data = sig_genes, shape = 21, size = 2, fill = "red", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-10, 10, 2)), limits = c(-10, 10)) +
  scale_y_continuous(breaks = c(seq(0, 10, 2.5)), limits = c(0, 10)) +
  labs(x = "log2(fold change)", y = "-log10(adjusted P-value)", colour = "Expression change") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(face = "bold", color = "black", size = 10),
        axis.text = element_text(color = "black", size = 9, face = "bold"),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", color = "black", size = 10),
        legend.text = element_text(face = "bold", color = "black", size = 9),
        legend.spacing.x = unit(0, "cm")  )

ggsave(p,filename = "/Plas/Plas.deseq.pdf",width=4,height=3)


