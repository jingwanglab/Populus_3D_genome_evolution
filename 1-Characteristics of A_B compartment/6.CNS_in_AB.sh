#!/usr/bin/bash
less Pkor.a.bed | awk '{print $0"\t"NR}'  > ./Pkor.a.bed;
less Pkor.b.bed | awk '{print $0"\t"NR}'  > ./Pkor.b.bed;
bedtools intersect -nonamecheck -a "$i".CNS_final.bed -b Pkor.a.bed  -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$8}' |awk '{a[$4]+=1}END{for(i in a){print a[i]"\t"i"\t""A"}}'| sort -k2,2n > a_cns.bed;
bedtools intersect -nonamecheck -a "$i".CNS_final.bed -b Pkor.b.bed  -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$8}' |awk '{a[$4]+=1}END{for(i in a){print a[i]"\t"i"\t""B"}}'| sort -k2,2n > b_cns.bed;
cat a_cns.bed b_cns.bed > ab.cns.csv;


