#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
  cd /usr_storage/atac/atac.tad/"$i";
  ln -s ../../atac.ab/"$i"/peak2gene_dis.bed ./;
  bedtools intersect -a boundary.bed -b peak2gene_dis.bed -nonamecheck -wa -wb |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t""boundary"}'  > "$i".boundary.peak2gene_dis.bed;
  bedtools intersect -a interior.bed -b peak2gene_dis.bed -nonamecheck -wa -wb |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t""interior"}'  > "$i".interior.peak2gene_dis.bed;
  cat "$i".boundary.peak2gene_dis.bed "$i".interior.peak2gene_dis.bed > "$i".boundary.interior.peak2gene_dis.csv;
done
  
