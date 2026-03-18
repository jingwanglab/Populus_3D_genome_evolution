#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   cd /w/00/g/g00/hic/tad.analysis/"$i"/;
   mkdir te.density;cd te.density;
   bedtools intersect -nonamecheck -a te.bed -b ../gene.density/tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.te.density.bed;
   bedtools intersect -nonamecheck -a te.bed -b ../gene.density/tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.te.density.bed;
   bedtools intersect -nonamecheck -a te.bed -b ../gene.density/taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.te.density.bed;
   cp /w/00/g/g00/hic/tad.analysis/changxuy/te.density/10bin.sh ./;
   bash 10bin.sh;
done