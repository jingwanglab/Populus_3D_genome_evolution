#install.packages('geomorph')

library(geomorph)
library(ape)
library(dplyr)

tree = read.tree('/PGLS/supergene.contree')

######### Fig 1b #########
gene_exp=read.csv('/all.ab.gene.expression.csv')
fit_gene <- extended.pgls(f1 = FPKM~compartment, 
                     data = gene_exp,
                     species = "species",
                     phy = tree) 
anova(fit_gene) 

te_exp=read.csv('/all.ab.TE.expression.csv')
fit_te <- extended.pgls(f1 = FPKM~compartment, 
                          data = te_exp,
                          species = "species",
                          phy = tree) 
anova(fit_te) 


######### Fig 1c #########
#chh
chh=read.csv('/all.ab.chh.csv')
fit_chh <- extended.pgls(f1 = methylation~status, 
                          data = chh,
                          species = "species",
                          phy = tree) 
anova(fit_chh) 

#chg
chg=read.csv('/all.ab.chg.csv')
fit_chh <- extended.pgls(f1 = methylation~status, 
                         data = chg,
                         species = "species",
                         phy = tree) 
anova(fit_chg) 

#cpg
cpg=read.csv('/all.ab.cpg.csv')
fit_chh <- extended.pgls(f1 = methylation~status, 
                         data = cpg,
                         species = "species",
                         phy = tree) 
anova(fit_cpg) 
