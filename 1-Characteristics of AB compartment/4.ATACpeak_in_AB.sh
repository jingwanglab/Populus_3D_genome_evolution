#!/usr/bin/bash

for i in Ppse Pwua Psze Plas Pyun Prot Pdav Pade Pkor Psim;
do
  cd /usr_storage/atac/atac.ab/"$i";
  bedtools intersect -a "$i".a.bed -b peak2gene_dis.bed -nonamecheck -wa -wb |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t""A"}'  > "$i".a.peak2gene_dis.bed;
  bedtools intersect -a "$i".b.bed -b peak2gene_dis.bed -nonamecheck -wa -wb |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t""B"}'  > "$i".b.peak2gene_dis.bed;
  cat "$i".a.peak2gene_dis.bed "$i".b.peak2gene_dis.bed > "$i".ab.peak2gene_dis.csv;
done
  
