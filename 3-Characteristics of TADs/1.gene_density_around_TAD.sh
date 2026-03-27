#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim
do
    cd /w/00/g/g00/hic/tad.analysis/"$i"/
    mkdir gene.density
    cd gene.density
    ln -s /w/00/g/g00/hic/compartment/gene.te.class/"$i"/gene.bed ./
    ####slide window into 10 bins for TAD up, body, down regions respectively
    bedtools makewindows -b ../tad.bed/"$i".tad.raw.bed -n 10 -i winnum > tad.10bin.bed
    bedtools makewindows -b ../tad.tad/"$i".tadup.50kb.bed -n 10 -i winnum > tadup.50kb.10bin.bed
    bedtools makewindows -b ../tad.tad/"$i".taddown.50kb.bed -n 10 -i winnum > taddown.50kb.10bin.bed
    bedtools intersect -nonamecheck -a gene.bed -b tadup.50kb.10bin.bed -wa -wb | awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.gene.density.bed
    bedtools intersect -nonamecheck -a gene.bed -b tad.10bin.bed -wa -wb | awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.gene.density.bed
    bedtools intersect -nonamecheck -a gene.bed -b taddown.50kb.10bin.bed -wa -wb | awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.gene.density.bed
    for file in tadup.gene.density.bed tad.gene.density.bed taddown.gene.density.bed
    do
        for bin in {1..10}
        do
            awk -v val="$bin" '$5 == val' "$file" | wc -l >> results.txt
        done
        wc -l < "$file" >> results.txt
    done
done