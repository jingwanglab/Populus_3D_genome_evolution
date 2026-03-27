#!/usr/bin/bash

ln -s /w/00/g/g00/hic/tad.analysis/Pkor/methylation.tad.20bin/Pkor.tadup.50kb.20bin.bed ./;
ln -s /w/00/g/g00/hic/tad.analysis/Pkor/methylation.tad.20bin/Pkor.tad.20bin.bed ./;
ln -s /w/00/g/g00/hic/tad.analysis/Pkor/methylation.tad.20bin/Pkor.taddown.50kb.20bin.bed ./;

bedtools shuffle -i CNS_final.bed -g Pkor.genome.size -chrom > CNS_random.bed;
### Calculate observed CNS enrichment in TADs and TAD boundaries (count CNS numbers)
bedtools intersect -nonamecheck -a CNS_final.bed -b Pkor.tadup.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > tadup.50kb.20bin.cns.csv;
bedtools intersect -nonamecheck -a CNS_final.bed -b Pkor.tad.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}'| awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > tad.20bin.cns.csv;
bedtools intersect -nonamecheck -a CNS_final.bed -b Pkor.taddown.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > taddown.50kb.20bin.cns.csv;
cat tadup.50kb.20bin.cns.csv tad.20bin.cns.csv taddown.50kb.20bin.cns.csv > tadup.tad.taddown.50kb.20bin.cns.csv;
### Calculate random CNS enrichment in TADs and TAD boundaries
bedtools intersect -nonamecheck -a CNS_random.bed -b Pkor.tadup.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > tadup.50kb.20bin.random.cns.csv;
bedtools intersect -nonamecheck -a CNS_random.bed -b Pkor.tad.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}'| awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > tad.20bin.random.cns.csv;
bedtools intersect -nonamecheck -a CNS_random.bed -b Pkor.taddown.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > taddown.50kb.20bin.random.cns.csv;
cat tadup.50kb.20bin.random.cns.csv tad.20bin.random.cns.csv taddown.50kb.20bin.random.cns.csv > tadup.tad.taddown.50kb.20bin.random.cns.csv;
paste tadup.tad.taddown.50kb.20bin.cns.csv tadup.tad.taddown.50kb.20bin.random.cns.csv | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3}' > tadup.tad.taddown.50kb.20bin.observed.random.cns.csv;
