#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pdav Pade Pkor Psim;
do
  cd /usr_storage/wjl/atac/atac.changed_tad/"$i";
  ln -s /usr_storage/atac/atac.ab/"$i"/peak2gene_dis.bed ./
  #######tad
  bedtools intersect -a "$i".high.con.tad.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""H-conserved""\t""tad"}' > "$i".high.con.tad.peak2gene_dis.txt;
  bedtools intersect -a "$i".con.tad.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""Conserved""\t""tad"}' > "$i".con.tad.peak2gene_dis.txt;
  bedtools intersect -a "$i".spe.tad.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""Diverged""\t""tad"}' > "$i".spe.tad.peak2gene_dis.txt;
  cat "$i".high.con.tad.peak2gene_dis.txt "$i".con.tad.peak2gene_dis.txt "$i".spe.tad.peak2gene_dis.txt >  "$i".high.con.spe.tad.peak2gene_dis.txt;
  #######boundary
  bedtools intersect -a "$i".high.con.boundary.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""H-conserved""\t""boundary"}' > "$i".high.con.boundary.peak2gene_dis.txt;
  bedtools intersect -a "$i".con.boundary.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""Conserved""\t""boundary"}' > "$i".con.boundary.peak2gene_dis.txt;
  bedtools intersect -a "$i".spe.boundary.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""Diverged""\t""boundary"}' > "$i".spe.boundary.peak2gene_dis.txt;
  cat "$i".high.con.boundary.peak2gene_dis.txt "$i".con.boundary.peak2gene_dis.txt "$i".spe.boundary.peak2gene_dis.txt >  "$i".high.con.spe.boundary.peak2gene_dis.txt;
  #######interior
  bedtools intersect -a "$i".high.con.interior.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""H-conserved""\t""interior"}' > "$i".high.con.interior.peak2gene_dis.txt;
  bedtools intersect -a "$i".con.interior.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""Conserved""\t""interior"}' > "$i".con.interior.peak2gene_dis.txt;
  bedtools intersect -a "$i".spe.interior.bed -b peak2gene_dis.bed -nonamecheck -wa -wb| awk '{print $0"\t""Diverged""\t""interior"}' > "$i".spe.interior.peak2gene_dis.txt;
  cat "$i".high.con.interior.peak2gene_dis.txt "$i".con.interior.peak2gene_dis.txt "$i".spe.interior.peak2gene_dis.txt >  "$i".high.con.spe.interior.peak2gene_dis.txt;
  cat "$i".high.con.spe.tad.peak2gene_dis.txt "$i".high.con.spe.boundary.peak2gene_dis.txt "$i".high.con.spe.interior.peak2gene_dis.txt > "$i".changed_tad.atac.csv;
done