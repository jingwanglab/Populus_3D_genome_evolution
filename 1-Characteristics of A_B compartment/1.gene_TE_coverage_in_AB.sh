#!/bin/bash

####gene and TE coverage in A/B
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   awk '{if ($4>0) print $1"\t"$2"\t"$3}' /w/00/g/g00/hic/juicer/"$i"/compartment/eigen.txt > "$i"/A.bed;####extract A compartment from juicer results
   awk '{if ($4<0) print $1"\t"$2"\t"$3}' /w/00/g/g00/hic/juicer/"$i"/compartment/eigen.txt > "$i"/B.bed;####extract B compartment from juicer results
   bedtools coverage -nonamecheck -a "$i"/A.bed -b "$i"/"$i".gene.bed > "$i"/A_gene.txt;#####calculate gene coverage in A compartment
   awk '$0=$0"\tA"' "$i"/A_gene.txt >"$i"/A1_gene.txt;
   bedtools coverage -nonamecheck -a "$i"/B.bed -b "$i"/"$i".gene.bed > "$i"/B_gene.txt;####calculate gene coverage in B compartment
   awk '$0=$0"\tB"' "$i"/B_gene.txt >"$i"/B1_gene.txt;
   bedtools coverage -nonamecheck -a "$i"/A.bed -b "$i"/"$i".TE.bed > "$i"/A_TE.txt;#####calculate TE coverage in A compartment
   awk '$0=$0"\tA"' "$i"/A_TE.txt >"$i"/A1_TE.txt;
   bedtools coverage -nonamecheck -a "$i"/B.bed -b "$i"/"$i".TE.bed > "$i"/B_TE.txt;####calculate TE coverage in B compartment
   awk '$0=$0"\tB"' "$i"/B_TE.txt >"$i"/B1_TE.txt;
   cat "$i"/A1_gene.txt "$i"/B1_gene.txt > "$i"/"$i".AB_gene.csv;
   cat "$i"/A1_TE.txt "$i"/B1_TE.txt > "$i"/"$i".AB_TE.csv;
done

