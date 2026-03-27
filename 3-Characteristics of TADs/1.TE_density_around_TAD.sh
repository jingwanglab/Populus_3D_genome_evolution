#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim
do
    cd /w/00/g/g00/hic/tad.analysis/"$i"/
    mkdir te.density
    cd te.density
    ln -s /w/00/g/g00/hic/compartment/gene.te.class/"$i"/TE.bed ./
    bedtools intersect -nonamecheck -a TE.bed -b ../gene.density/tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.te.density.bed
    bedtools intersect -nonamecheck -a TE.bed -b ../gene.density/tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.te.density.bed
    bedtools intersect -nonamecheck -a TE.bed -b ../gene.density/taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.te.density.bed
    for file in tadup.te.density.bed tad.te.density.bed taddown.te.density.bed
    do
        for bin in {1..10}
        do
            awk -v val="$bin" '$5 == val' "$file" | wc -l >> results.txt
        done
        wc -l < "$file" >> results.txt
    done
done