#!/usr/bin/bash 

for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do 
   cd /w/00/g/g00/hic/compartment/gene.expression/"$i";
   less genes.fpkm_tracking |grep -v Contig |awk '{print $7"\t"$10}'|sed 1d |sed 's/:/\t/g'|less |awk -F '\t' '{split($2,a,"-")}{print $1"\t"a[1]+1"\t"a[2]"\t"$3}'|sort -k1,1V -k2,2n -k3,3n  > "$i".gene.fpkm.bed;
   sort -k1,1V -k2,2n -k3,3n "$i".gene.fpkm.bed >"$i".gene.fpkm.sort.bed;
   bedtools merge -i "$i".a.bed |awk '$0=$0"\tA"'> "$i".a.merge.bed;####merge consecutive regions in the file
   bedtools merge -i "$i".b.bed |awk '$0=$0"\tB"'> "$i".b.merge.bed;
   cat "$i".a.merge.bed "$i".b.merge.bed | sort -k1,1V -k2,2n -k3,3n > "$i".ab.sort.bed;####merge A and B compartments into one file and sort
   bedtools intersect -nonamecheck -a "$i".gene.fpkm.sort.bed -b "$i".ab.sort.bed -wa > 1.txt;
   sort 1.txt | uniq -d >2.txt; ###save duplicate regions to a new file, check which genes exist in both A and B compartments
   less 2.txt|awk '{print $1";"$2";"$3}'|while read id;do echo "sed -i '/$id/d' "$i".gene.fpkm.sort.bed">>tp.sh;done;###delete genes that exist in both A and B compartments
   bash tp.sh;
   less "$i".gene.fpkm.sort.bed|sed 's/;/\t/g' > "$i".gene.fpkm.sort.uniq.bed;
   bedtools intersect -nonamecheck -a "$i".gene.fpkm.sort.uniq.bed -b "$i".ab.sort.bed -wa -wb|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > "$i".ab.gene.fpkm.csv;
   less "$i".ab.gene.fpkm.csv| awk '$0=$0"\t'$i'"' >"$i".ab.gene.fpkm.name.csv;
done

