# Populus_3D_genome_evolution

Script for Shi et al. (2026) Dynamic reorganization of three-dimensional genome architecture during Populus diversification. Nature Ecology & Evolution. accepted

1. Characteristics of A/B compartment
This section contains a series of shell scripts to quantify multi-omics and genomic feature enrichment between A and B chromatin compartments across poplar genomes, including gene/TE coverage, expression, DNA methylation, accessible chromatin peaks, pan-gene distribution and conserved non-coding sequences (CNS).

1.1 gene_TE_coverage_in_AB.sh - Script to calculate the coverage ratios of protein-coding genes and transposable elements (TEs) within A and B compartments.

1.2 Gene & TE expression in A/B compartments
gene_expression_in_AB.sh - Script to extract gene expression values and annotate genes with their corresponding A/B compartment types.
TE_expression_in_AB.sh - Script to extract TE expression values and label TEs according to their A/B compartment localization.

1.3 methylation_in_AB.sh - Script to compute average DNA methylation levels in A and B compartments.

1.4 ATACpeak_in_AB.sh - Script to characterize peak-to-gene distances for accessible chromatin peaks located in A and B compartments.

1.5 pan-gene_in_AB.sh - Script to count the number of each pan-gene subtype distributed in A and B compartments.

1.6 CNS_in_AB.sh - Script to count the number of CNS for each compartment bin and annotate the compartment type (A/B) for each record.

2. Evolution of A/B compartment
This section contains shell scripts to systematically analyze the evolutionary dynamics of A/B chromatin compartments across poplar species, including homologous compartment bin identification, classification of conserved/changed compartments, and multi-omics feature comparison between conserved and changed compartments.

2.1 homologous_bin_identify.sh - Script to identify cross-species homologous 40 kb compartment bins and classify conserved and changed A/B compartment status between Populus koreana and other poplar species.

2.2 obtain_AB_changed_file.sh - Script to classify homologous 40 kb bins into four compartment transition types (A-to-B, B-to-A, A-to-A, B-to-B) using paired compartment status files from P. koreana and each poplar species.

2.3 pangene_in_changed_AB.sh - Script to count the number of core, dispensable and private pan-genes distributed across four compartment transition types for P. koreana and each paired poplar species.

2.4 TEandSV_coverage_in_changed_AB.sh - Script to calculate TE and SV coverage in evolutionarily conserved and dynamically changed A/B compartment regions.

2.5 gene_exp_in_changed_AB.sh - Script to extract gene expression values of genes fully embedded within four types of compartment transition regions (A-to-B, B-to-A, A-to-A, B-to-B) for P. koreana and each comparative poplar species.

2.6 methylation_in_changed_AB.sh - Script to calculate average DNA methylation levels over gene bodies and 2-kb flanking regions across four A/B compartment transition states.

3. Characteristics of TADs
This section contains shell scripts for multi-omics and genomic profiling of topologically associating domains (TADs), including gene/TE density, expression, chromatin accessibility, DNA methylation, pan-gene distribution and CNS enrichment across TAD bodies and flanking regions.

3.1 Gene & TE density around TADs
gene_density_around_TAD.sh - Script to quantify gene density across 10 equal bins of TAD bodies and 50-kb flanking regions.
TE_density_around_TAD.sh - Script to quantify TE density across 10 equal bins of TAD bodies and 50-kb flanking regions.

3.2 gene_TE_expression_in_TADs.sh - Script to extract gene and TE expression levels within TAD interior and 20-kb boundary regions, retaining only fully embedded genes and TEs for comparison.

3.3 ATACpeak_in_TADs.sh - Script to characterize peak-to-gene distances of accessible chromatin peaks in TAD interior and boundary regions.

3.4 methylation_around_TADs.sh - Script to generate bin-level DNA methylation metaprofiles (CpG, CHG, CHH) across TAD upstream 50-kb flanks, TAD bodies, and downstream 50-kb flanks, with each region divided into 20 equal bins.

3.5 pangene_density_around_TADs.sh - Script to profile pan-gene distribution across TAD upstream 50-kb flanks, TAD bodies, and downstream 50-kb flanks by counting core, dispensable, and private pan-genes in 10 equal bins.

3.6 CNS_enrichment_around_TAD.sh - Script to quantify CNS enrichment across TAD flanking and body regions by comparing observed CNS counts with randomized background signals.

4. Evolution of TADs 
This section systematically explores the multi-omics evolutionary characteristics of TADs by classifying TADs into three evolutionary grades: highly conserved, conserved, and diverged TADs. Combined with CNSs, DNA methylation, gene expression, and chromatin accessibility data, it compares the differences in genomic and epigenetic landscapes among TADs with different evolutionary stabilities.

4.1 CNS.sh - Script to compare CNS enrichment across three evolutionary TAD categories and their 50-kb upstream and downstream flanking regions, with random CNS shuffling for background normalization.

4.2 methylation.sh - Script to generate bin-level DNA methylation metaprofiles for three types of evolutionary TADs and their 50-kb flanking regions, calculating average CpG, CHG, and CHH methylation levels per bin.

4.3 gene_expression.sh - Script to extract gene expression levels for three evolutionarily TAD categories, separating signals from TAD interior and TAD boundary regions.

4.4 ATAC.sh - Script to extract peak-to-gene distances for accessible chromatin peaks in full TADs, TAD boundaries, and TAD interior regions across three TAD evolutionary categories.

5. SV and TE in TADs
This section focuses on the distribution and enrichment patterns of SVs and TEs within TADs with different evolutionary stabilities.

5.1 sv_coverage_in_changedTAD.sh - Script to calculate and summarize SV coverage ratios in full TAD, TAD boundary, and TAD interior regions across highly-conserved, conserved, and diverged TADs.

5.2 sv_enrichment_in_changedTAD.sh - Script to calculate SV breakpoint fold enrichment across 20 equal bins of TAD bodies and 50-kb flanks for three TAD evolutionary types using shuffled random breakpoints as background control.

5.3 te_coverage_in_changedTAD.sh - Script to calculate TE coverage ratios in full TAD, TAD boundary and TAD interior regions across three evolutionary TAD categories.

5.4 with.without_SV_TE_coverage.py - Python script to stratify TAD boundaries into SV-present and SV-absent groups and quantify coverage of different TE subfamilies for comparative analysis.

5.5 diff_sv_in_changedTAD.sh - Script to overlap four types of structural variants (indel, inversion, duplication, translocation) with conserved and diverged TADs to characterize differential SV distribution patterns.

5.6.1 random_boundary_SV_overlap.sh - Script to generate 1000 sets of shuffled random TAD boundaries and calculate overlaps with four categories of SVs as null background.
5.6.2 random_boundary_SV_stats.py - Python script to statistically summarize SV overlapping counts across all random boundary permutations.

6. Other analysis
This section integrates the workflow scripts for identifying conserved noncoding sequences (CNSs) and auxiliary R scripts for transcriptome differential expression, gene functional enrichment and cross-species phylogenetic comparative statistical analysis.

CNS identify - Scripts to conduct multi-species whole-genome alignment, identify, filter and preprocess CNS, and generate standardized CNS BED files for downstream enrichment analyses of A/B compartments and TADs.

Deseq.R - R script performing DESeq2-based differential gene expression analysis between treatment and control groups, and generating volcano plots for DEG visualization.

GO_enrichment.R - R script utilizing topGO to conduct biological process GO enrichment for genes derived from three TAD evolutionary categories, perform BH P-value correction, and visualize functional enrichment via bubble plots.

PGLS.R - R script implementing phylogenetic generalized least squares (PGLS) regression to test interspecies correlations between genomic traits while controlling phylogenetic relatedness.

7. Figure plot
This section contains all customized R plotting scripts used for generating formal manuscript figures, including Figure 1a–f, Figure 2a–h, Figure 3a–i, Figure 4c–h, Figure 5a–d, and Figure 6a.
