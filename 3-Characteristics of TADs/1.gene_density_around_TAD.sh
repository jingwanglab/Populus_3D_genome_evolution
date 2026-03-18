#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   cd /w/00/g/g00/hic/tad.analysis/"$i"/;
   mkdir gene.density;cd gene.density;
   cp /w/00/g/g00/hic/compartment/gene.class/"$i"/geneid.bed ./;
   ####slide window into 10 bins for TAD up, body, down regions respectively
   bedtools makewindows -b ../tad.bed/"$i".tad.raw.bed -n 10 -i winnum > tad.10bin.bed;
   bedtools makewindows -b ../tad.tad/"$i".tadup.50kb.bed -n 10 -i winnum > tadup.50kb.10bin.bed;
   bedtools makewindows -b ../tad.tad/"$i".taddown.50kb.bed -n 10 -i winnum > taddown.50kb.10bin.bed;
   bedtools intersect -nonamecheck -a geneid.bed -b tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.gene.density.bed;
   bedtools intersect -nonamecheck -a geneid.bed -b tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.gene.density.bed;
   bedtools intersect -nonamecheck -a geneid.bed -b taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.gene.density.bed;
   cp /w/00/g/g00/hic/tad.analysis/changxuy/gene.density/10bin.sh ./;
   bash 10bin.sh;
done
