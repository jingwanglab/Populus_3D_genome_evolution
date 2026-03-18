#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   cd /w/00/g/g00/hic/tad.analysis/"$i"/gene.TE.expression;
   cp /w/00/g/g00/hic/compartment/gene.expression/"$i"/"$i".gene.fpkm.sort.bed ./;
   bedtools intersect -nonamecheck -a "$i".gene.fpkm.sort.bed -b tad.interior.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $4}'|awk '$0=$0"\ttad"'|awk '$0=$0"\t'$i'"' > tad.interior.fpkm.csv;
   bedtools intersect -nonamecheck -a "$i".gene.fpkm.sort.bed -b boundary.20kb.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $4}'|awk '$0=$0"\tboundary"'|awk '$0=$0"\t'$i'"'>boundary.20kb.fpkm.csv;
   cat tad.interior.fpkm.csv boundary.20kb.fpkm.csv > tad.boundary.20kb.fpkm.csv;
   
   cp /w/00/g/g00/hic/compartment/te.expression/"$i"/"$i".te.fpkm.sort.bed ./;
   bedtools intersect -nonamecheck -a "$i".te.fpkm.sort.bed -b tad.interior.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $4}'|awk '$0=$0"\ttad"'|awk '$0=$0"\t'$i'"' > tad.interior.te.fpkm.csv;
   bedtools intersect -nonamecheck -a "$i".te.fpkm.sort.bed -b boundary.20kb.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $4}'|awk '$0=$0"\tboundary"'|awk '$0=$0"\t'$i'"'>boundary.20kb.te.fpkm.csv;
   cat tad.interior.te.fpkm.csv boundary.20kb.te.fpkm.csv > tad.boundary.20kb.te.fpkm.csv;
done


