#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   cd /w/00/g/g00/hic/tad.analysis/"$i"/gene.TE.expression;
   ln -s /w/00/g/g00/hic/compartment/gene.expression/"$i"/"$i".gene.tpm.sort.bed ./;
   bedtools intersect -nonamecheck -a "$i".gene.tpm.sort.bed -b tad.interior.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $4}'|awk '$0=$0"\ttad"'|awk '$0=$0"\t'$i'"' > tad.interior.tpm.csv;
   bedtools intersect -nonamecheck -a "$i".gene.tpm.sort.bed -b boundary.20kb.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $4}'|awk '$0=$0"\tboundary"'|awk '$0=$0"\t'$i'"'>boundary.20kb.tpm.csv;
   cat tad.interior.tpm.csv boundary.20kb.tpm.csv > tad.boundary.20kb.tpm.csv;
   
   cp /w/00/g/g00/hic/compartment/te.expression/"$i"/"$i".te.tpm.sort.bed ./;
   bedtools intersect -nonamecheck -a "$i".te.tpm.sort.bed -b tad.interior.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $4}'|awk '$0=$0"\ttad"'|awk '$0=$0"\t'$i'"' > tad.interior.te.tpm.csv;
   bedtools intersect -nonamecheck -a "$i".te.tpm.sort.bed -b boundary.20kb.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $4}'|awk '$0=$0"\tboundary"'|awk '$0=$0"\t'$i'"'>boundary.20kb.te.tpm.csv;
   cat tad.interior.te.tpm.csv boundary.20kb.te.tpm.csv > tad.boundary.20kb.te.tpm.csv;
done


