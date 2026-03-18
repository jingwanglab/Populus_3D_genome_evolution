####11个杨属物种的core、dispensable、private基因在tad和tad边界的gene.density
1.
export PATH="/data/apps/bedtools/bedtools-2.26.0/bin:$PATH"
#!/usr/bin/bash
for i in changxuy changyey chuany dayey diany qingxiy qiongdaoy shany xiangy xiangyey xiaoyey;
do
   cd /w/00/g/g00/wjl/hic/tad.analysis/"$i"/;
   mkdir core.dis.pri.gene.density; cd core.dis.pri.gene.density;
   ####将tad和tad上下游与core、dispensable、private基因的.bed文件取交集，只保留完全包含在区域中的gene
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_core_gene.bed -b ../gene.density/tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.core.gene.density.bed;
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_core_gene.bed -b ../gene.density/tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.core.gene.density.bed;
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_core_gene.bed -b ../gene.density/taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.core.gene.density.bed;
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_dis_gene.bed -b ../gene.density/tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.dis.gene.density.bed;
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_dis_gene.bed -b ../gene.density/tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.dis.gene.density.bed;
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_dis_gene.bed -b ../gene.density/taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.dis.gene.density.bed;
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_pri_gene.bed -b ../gene.density/tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tadup.pri.gene.density.bed;
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_pri_gene.bed -b ../gene.density/tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > tad.pri.gene.density.bed;
   bedtools intersect -nonamecheck -a ../core.dis.pri.gene/"$i"_pri_gene.bed -b ../gene.density/taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > taddown.pri.gene.density.bed;
done

2.
vi 10bin.sh
#!/usr/bin/bash
for i in $(seq 1 10);
do
   less tadup.core.gene.density.bed| awk '{if($5=="'$i'") print $0}'|wc -l>>up.core.txt;
   less tad.core.gene.density.bed |awk '{if($5=="'$i'") print $0}'|wc -l>>tad.core.txt;
   less taddown.core.gene.density.bed|awk '{if($5=="'$i'") print $0}'|wc -l>>down.core.txt;
   less tadup.dis.gene.density.bed| awk '{if($5=="'$i'") print $0}'|wc -l>>up.dis.txt;
   less tad.dis.gene.density.bed |awk '{if($5=="'$i'") print $0}'|wc -l>>tad.dis.txt;
   less taddown.dis.gene.density.bed|awk '{if($5=="'$i'") print $0}'|wc -l>>down.dis.txt;
   less tadup.pri.gene.density.bed| awk '{if($5=="'$i'") print $0}'|wc -l>>up.pri.txt;
   less tad.pri.gene.density.bed |awk '{if($5=="'$i'") print $0}'|wc -l>>tad.pri.txt;
   less taddown.pri.gene.density.bed|awk '{if($5=="'$i'") print $0}'|wc -l>>down.pri.txt;
done

vi 1.sh
wc -l tadup.core.gene.density.bed >> 1.txt;
wc -l tad.core.gene.density.bed>> 1.txt;
wc -l taddown.core.gene.density.bed>> 1.txt;
wc -l tadup.dis.gene.density.bed>> 1.txt;
wc -l tad.dis.gene.density.bed>> 1.txt;
wc -l taddown.dis.gene.density.bed>> 1.txt;
wc -l tadup.pri.gene.density.bed>> 1.txt;
wc -l tad.pri.gene.density.bed>> 1.txt;
wc -l taddown.pri.gene.density.bed>> 1.txt;

#!/usr/bin/bash
for i in changyey chuany dayey diany qingxiy qiongdaoy shany xiangy xiangyey xiaoyey;
do
  cd /w/00/g/g00/wjl/hic/tad.analysis/"$i"/core.dis.pri.gene.density; 
  cp ../../changxuy/core.dis.pri.gene.density/10bin.sh ./;
  cp ../../changxuy/core.dis.pri.gene.density/1.sh ./;
  bash 10bin.sh;
  bash 1.sh;
  cat up.core.txt tad.core.txt down.core.txt > core.txt;
  cat up.dis.txt tad.dis.txt down.dis.txt > dis.txt;
  cat up.pri.txt tad.pri.txt down.pri.txt > pri.txt;
done



  


