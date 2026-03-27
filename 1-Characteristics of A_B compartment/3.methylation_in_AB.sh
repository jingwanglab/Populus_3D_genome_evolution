#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do 
   cd /w/00/g/g00/wjl/hic/compartment/methylation/"$i";
   #########CHG
   bedtools intersect -nonamecheck -a CHG.bed -b "$i".a.bed -wa -wb |awk '{print $4"\t"$10}' |awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A"}}'  > "$i".a.chg.csv;
   bedtools intersect -nonamecheck -a CHG.bed -b "$i".b.bed -wa -wb |awk '{print $4"\t"$10}' |awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B"}}'  > "$i".b.chg.csv;
   cat "$i".a.chg.csv "$i".b.chg.csv > "$i".ab.chg.csv;
   #########CHH
   bedtools intersect -nonamecheck -a CHH.bed -b "$i".a.bed -wa -wb |awk '{print $4"\t"$10}' |awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A"}}'  > "$i".a.chh.csv;
   bedtools intersect -nonamecheck -a CHH.bed -b "$i".b.bed -wa -wb |awk '{print $4"\t"$10}' |awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B"}}'  > "$i".b.chh.csv;
   cat "$i".a.chh.csv "$i".b.chh.csv > "$i".ab.chh.csv;
   #########CpG
   bedtools intersect -nonamecheck -a CpG.bed -b "$i".a.bed -wa -wb |awk '{print $4"\t"$10}' |awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A"}}'  > "$i".a.cpg.csv;
   bedtools intersect -nonamecheck -a CpG.bed -b "$i".b.bed -wa -wb |awk '{print $4"\t"$10}' |awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B"}}'  > "$i".b.cpg.csv;
   cat "$i".a.cpg.csv "$i".b.cpg.csv > "$i".ab.cpg.csv;
done 



























