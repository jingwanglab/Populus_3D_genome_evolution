#!/usr/bin/bash

for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   cd /w/00/g/g00/hic/compartment/te.expression/"$i";
   less te.fpkm_tracking |grep -v Contig |sed 1d |awk '{if($10>0) print $7"\t"$10}'|sed 's/:/\t/g'|less |awk -F '\t' '{split($2,a,"-")}{print $1"\t"a[1]+1"\t"a[2]"\t"$3}'|sort -k1,1V -k2,2n -k3,3n  > "$i".te.fpkm.bed;
   sort -k1,1V -k2,2n -k3,3n "$i".te.fpkm.bed >"$i".te.fpkm.sort.bed;
   bedtools intersect -nonamecheck -a "$i".te.fpkm.sort.bed -b "$i".ab.sort.bed -wa > 1.txt;
   sort 1.txt | uniq -d >2.txt; ###save duplicate regions to a new file, check which TEs exist in both A and B compartments
   less 2.txt|awk '{print $1";"$2";"$3}'|while read id;do echo "sed -i '/$id/d' "$i".te.fpkm.sort.bed">>tp.sh;done;###delete TEs that exist in both A and B compartments
   bash tp.sh;
   less "$i".te.fpkm.sort.bed|sed 's/;/\t/g' > "$i".te.fpkm.sort.uniq.bed;
   bedtools intersect -nonamecheck -a "$i".te.fpkm.sort.uniq.bed -b "$i".ab.sort.bed -wa -wb|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > "$i".ab.te.fpkm.csv;
   less "$i".ab.te.fpkm.csv| awk '$0=$0"\t'$i'"' >"$i".ab.te.fpkm.name.csv
done
