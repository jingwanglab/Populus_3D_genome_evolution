#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim;
do
   cd /w/00/g/g00/hic/tad.analysis/changed_tad.gene_expression/;
   mkdir "$i"; cd "$i";
   ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/"$i".high.con.tad.bed;
   ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/"$i".con.tad.bed;
   ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/"$i".spe.tad.bed;
   ln -s /w/00/g/g00/hic/compartment/gene.expression/"$i"/"$i".gene.tpm.sort.bed ./;
   bedtools intersect -a "$i".gene.tpm.sort.bed -b "$i".high.con.tad.bed -wa -wb -nonamecheck | awk '{if($2>$6 && $3<$7) print $4"\t"$5"\t"$6"\t"$7"\t""H-conserved"}'  > "$i".high.con.tad.gene_exp;
   bedtools intersect -a "$i".gene.tpm.sort.bed -b "$i".con.tad.bed -wa -wb -nonamecheck | awk '{if($2>$6 && $3<$7) print $4"\t"$5"\t"$6"\t"$7"\t""Conserved"}'  > "$i".con.tad.gene_exp;
   bedtools intersect -a "$i".gene.tpm.sort.bed -b "$i".spe.tad.bed -wa -wb -nonamecheck | awk '{if($2>$6 && $3<$7) print $4"\t"$5"\t"$6"\t"$7"\t""Diverged"}'  > "$i".spe.tad.gene_exp;
   cat "$i".high.con.tad.gene_exp "$i".con.tad.gene_exp "$i".spe.tad.gene_exp > "$i".changed_tad.gene_expression.csv;
done
