#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   cd /w/00/g/g00/hic/compartment/pangene/"$i";
   bedtools intersect -nonamecheck -a "$i"_core_gene.bed -b "$i".ab.sort.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > ab.core.bed;
   bedtools intersect -nonamecheck -a "$i"_dis_gene.bed -b "$i".ab.sort.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > ab.dis.bed;
   bedtools intersect -nonamecheck -a "$i"_pri_gene.bed -b "$i".ab.sort.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > ab.pri.bed;
   less ab.core.bed|awk '{if($5=="A") print $0}'|wc -l >>stat.txt;
   less ab.core.bed|awk '{if($5=="B") print $0}'|wc -l >>stat.txt;
   less ab.dis.bed|awk '{if($5=="A") print $0}'|wc -l >>stat.txt;
   less ab.dis.bed|awk '{if($5=="B") print $0}'|wc -l >>stat.txt;
   less ab.pri.bed|awk '{if($5=="A") print $0}'|wc -l >>stat.txt;
   less ab.pri.bed|awk '{if($5=="B") print $0}'|wc -l >>stat.txt;
done
