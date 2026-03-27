#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   cd /w/00/g/g00/wjl/tad.analysis/"$i"/
   mkdir pangene.density; cd pangene.density
   ln -s /w/00/g/g00/compartment/pangene/"$i"/"$i"_core_gene.bed ./
   ln -s /w/00/g/g00/compartment/pangene/"$i"/"$i"_dis_gene.bed ./
   ln -s /w/00/g/g00/compartment/pangene/"$i"/"$i"_pri_gene.bed ./
   
   bedtools intersect -nonamecheck -a "$i"_core_gene.bed -b ../gene.density/tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.core.gene.density.bed
   bedtools intersect -nonamecheck -a "$i"_core_gene.bed -b ../gene.density/tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.core.gene.density.bed
   bedtools intersect -nonamecheck -a "$i"_core_gene.bed -b ../gene.density/taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.core.gene.density.bed
   bedtools intersect -nonamecheck -a "$i"_dis_gene.bed -b ../gene.density/tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.dis.gene.density.bed
   bedtools intersect -nonamecheck -a "$i"_dis_gene.bed -b ../gene.density/tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.dis.gene.density.bed
   bedtools intersect -nonamecheck -a "$i"_dis_gene.bed -b ../gene.density/taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.dis.gene.density.bed
   bedtools intersect -nonamecheck -a "$i"_pri_gene.bed -b ../gene.density/tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.pri.gene.density.bed
   bedtools intersect -nonamecheck -a "$i"_pri_gene.bed -b ../gene.density/tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.pri.gene.density.bed
   bedtools intersect -nonamecheck -a "$i"_pri_gene.bed -b ../gene.density/taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.pri.gene.density.bed
   for bin in $(seq 1 10)
   do
      less tadup.core.gene.density.bed| awk '{if($5=="'$bin'") print $0}'|wc -l>>up.core.txt
      less tad.core.gene.density.bed |awk '{if($5=="'$bin'") print $0}'|wc -l>>tad.core.txt
      less taddown.core.gene.density.bed|awk '{if($5=="'$bin'") print $0}'|wc -l>>down.core.txt
      less tadup.dis.gene.density.bed| awk '{if($5=="'$bin'") print $0}'|wc -l>>up.dis.txt
      less tad.dis.gene.density.bed |awk '{if($5=="'$bin'") print $0}'|wc -l>>tad.dis.txt
      less taddown.dis.gene.density.bed|awk '{if($5=="'$bin'") print $0}'|wc -l>>down.dis.txt
      less tadup.pri.gene.density.bed| awk '{if($5=="'$bin'") print $0}'|wc -l>>up.pri.txt
      less tad.pri.gene.density.bed |awk '{if($5=="'$bin'") print $0}'|wc -l>>tad.pri.txt
      less taddown.pri.gene.density.bed|awk '{if($5=="'$bin'") print $0}'|wc -l>>down.pri.txt
   done
   cat up.core.txt tad.core.txt down.core.txt > core.txt
   cat up.dis.txt tad.dis.txt down.dis.txt > dis.txt
   cat up.pri.txt tad.pri.txt down.pri.txt > pri.txt
done

