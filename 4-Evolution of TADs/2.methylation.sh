#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim;
do
   cd /w/00/g/g00/hic/tad.analysis/changed_tad.methylation/;
   mkdir "$i"; cd "$i";
   ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/"$i".high.con.tad.bed;
   ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/"$i".con.tad.bed;
   ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/"$i".spe.tad.bed;
   ## Create soft links to methylation data
   ln -s /w/00/g/g00/hic/methylation/"$i"/CHG.results ./;
   ln -s /w/00/g/g00/hic/methylation/"$i"/CHH.results ./;
   ln -s /w/00/g/g00/hic/methylation/"$i"/CpG.results ./;
   ## Obtain upstream and downstream 50kb positional information for highly conserved, conserved, and diverged TADs
   less "$i".high.con.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > "$i".high.con.tadup.50kb.bed;
   less "$i".high.con.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > "$i".high.con.taddown.50kb.bed;
   less "$i".con.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > "$i".con.tadup.50kb.bed;
   less "$i".con.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > "$i".con.taddown.50kb.bed;
   less "$i".spe.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > "$i".spe.tadup.50kb.bed;
   less "$i".spe.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > "$i".spe.taddown.50kb.bed;
   ## Divide highly conserved, conserved, and diverged TADs and their upstream/downstream 50kb regions into 20 bins
   bedtools makewindows -b "$i".high.con.tad.bed -n 20 -i winnum > "$i".high.con.tad.20bin.bed;
   bedtools makewindows -b "$i".high.con.tadup.50kb.bed -n 20 -i winnum > "$i".high.con.tadup.50kb.20bin.bed;
   bedtools makewindows -b "$i".high.con.taddown.50kb.bed -n 20 -i winnum > "$i".high.con.taddown.50kb.20bin.bed;
   bedtools makewindows -b "$i".con.tad.bed -n 20 -i winnum > "$i".con.tad.20bin.bed;
   bedtools makewindows -b "$i".con.tadup.50kb.bed -n 20 -i winnum > "$i".con.tadup.50kb.20bin.bed;
   bedtools makewindows -b "$i".con.taddown.50kb.bed -n 20 -i winnum > "$i".con.taddown.50kb.20bin.bed;
   bedtools makewindows -b "$i".spe.tad.bed -n 20 -i winnum > "$i".spe.tad.20bin.bed;
   bedtools makewindows -b "$i".spe.tadup.50kb.bed -n 20 -i winnum > "$i".spe.tadup.50kb.20bin.bed;
   bedtools makewindows -b "$i".spe.taddown.50kb.bed -n 20 -i winnum > "$i".spe.taddown.50kb.20bin.bed;
   ## Obtain methylation levels for each bin in highly conserved, conserved, and diverged TADs and their upstream/downstream 50kb regions, then merge
   ####CHG
   #######high conserved
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".high.con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".high.con.tad.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".high.con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".high.con.tadup.50kb.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".high.con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}'| sort -k2,2n > "$i".high.con.taddown.50kb.20bin.chg.bed;
   cat "$i".high.con.tadup.50kb.20bin.chg.bed "$i".high.con.tad.20bin.chg.bed  "$i".high.con.taddown.50kb.20bin.chg.bed | awk '{print $1"\t"NR"\t"$3"\t""H-conserved"}'> "$i".high.con.uptaddown.20bin.chg.csv;
   #######conserved
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".con.tad.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".con.tadup.50kb.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".con.taddown.50kb.20bin.chg.bed;
   cat "$i".con.tadup.50kb.20bin.chg.bed "$i".con.tad.20bin.chg.bed  "$i".con.taddown.50kb.20bin.chg.bed | awk '{print $1"\t"NR"\t"$3"\t""Conserved"}'> "$i".con.uptaddown.20bin.chg.csv;
   #######diverged
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".spe.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".spe.tad.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".spe.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".spe.tadup.50kb.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.results -b "$i".spe.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".spe.taddown.50kb.20bin.chg.bed;
   cat "$i".spe.tadup.50kb.20bin.chg.bed "$i".spe.tad.20bin.chg.bed  "$i".spe.taddown.50kb.20bin.chg.bed | awk '{print $1"\t"NR"\t"$3"\t""Diverged"}'> "$i".spe.uptaddown.20bin.chg.csv;
   cat "$i".high.con.uptaddown.20bin.chg.csv "$i".con.uptaddown.20bin.chg.csv "$i".spe.uptaddown.20bin.chg.csv > "$i".high.con.spe.uptaddown.20bin.chg.csv;
   ####CHH
   #######high conserved
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".high.con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".high.con.tad.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".high.con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".high.con.tadup.50kb.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".high.con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".high.con.taddown.50kb.20bin.chh.bed;
   cat "$i".high.con.tadup.50kb.20bin.chh.bed "$i".high.con.tad.20bin.chh.bed  "$i".high.con.taddown.50kb.20bin.chh.bed | awk '{print $1"\t"NR"\t"$3"\t""H-conserved"}'> "$i".high.con.uptaddown.20bin.chh.csv;
   #######conserved
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".con.tad.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".con.tadup.50kb.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".con.taddown.50kb.20bin.chh.bed;
   cat "$i".con.tadup.50kb.20bin.chh.bed "$i".con.tad.20bin.chh.bed  "$i".con.taddown.50kb.20bin.chh.bed | awk '{print $1"\t"NR"\t"$3"\t""Conserved"}'> "$i".con.uptaddown.20bin.chh.csv;
   #######diverged
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".spe.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".spe.tad.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".spe.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".spe.tadup.50kb.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.results -b "$i".spe.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}'| sort -k2,2n > "$i".spe.taddown.50kb.20bin.chh.bed;
   cat "$i".spe.tadup.50kb.20bin.chh.bed "$i".spe.tad.20bin.chh.bed  "$i".spe.taddown.50kb.20bin.chh.bed | awk '{print $1"\t"NR"\t"$3"\t""Diverged"}'> "$i".spe.uptaddown.20bin.chh.csv;
   cat "$i".high.con.uptaddown.20bin.chh.csv "$i".con.uptaddown.20bin.chh.csv "$i".spe.uptaddown.20bin.chh.csv > "$i".high.con.spe.uptaddown.20bin.chh.csv;
   ####CpG
   #######high conserved
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".high.con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".high.con.tad.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".high.con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".high.con.tadup.50kb.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".high.con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".high.con.taddown.50kb.20bin.cpg.bed;
   cat "$i".high.con.tadup.50kb.20bin.cpg.bed "$i".high.con.tad.20bin.cpg.bed  "$i".high.con.taddown.50kb.20bin.cpg.bed | awk '{print $1"\t"NR"\t"$3"\t""H-conserved"}'> "$i".high.con.uptaddown.20bin.cpg.csv;
   #######conserved
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".con.tad.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".con.tadup.50kb.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".con.taddown.50kb.20bin.cpg.bed;
   cat "$i".con.tadup.50kb.20bin.cpg.bed "$i".con.tad.20bin.cpg.bed  "$i".con.taddown.50kb.20bin.cpg.bed | awk '{print $1"\t"NR"\t"$3"\t""Conserved"}'> "$i".con.uptaddown.20bin.cpg.csv;
   #######diverged
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".spe.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".spe.tad.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".spe.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".spe.tadup.50kb.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.results -b "$i".spe.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}'| sort -k2,2n > "$i".spe.taddown.50kb.20bin.cpg.bed;
   cat "$i".spe.tadup.50kb.20bin.cpg.bed "$i".spe.tad.20bin.cpg.bed  "$i".spe.taddown.50kb.20bin.cpg.bed | awk '{print $1"\t"NR"\t"$3"\t""Diverged"}'> "$i".spe.uptaddown.20bin.cpg.csv;
   cat "$i".high.con.uptaddown.20bin.cpg.csv "$i".con.uptaddown.20bin.cpg.csv "$i".spe.uptaddown.20bin.cpg.csv > "$i".high.con.spe.uptaddown.20bin.cpg.csv;
done
