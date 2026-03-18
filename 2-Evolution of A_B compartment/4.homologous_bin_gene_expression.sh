以香杨为参考，发生AB区室转换（A-B,B-A,A-A,B-B)的基因表达情况--homologous.bin.40kb.minmatch0.7
#!/usr/bin/bash
export PATH="/data/apps/bedtools/bedtools-2.26.0/bin:$PATH"
for i in changyey dayey diany qingxiy qiongdaoy shany xiangyey xiaoyey;
do
   cd /w/00/g/g00/wjl/hic/compartment/ab.change.homologous.bin.gene.expression/
   mkdir xiangy."$i"/;
   cd xiangy."$i";
   mkdir 40kb.minm0.7;
   cd 40kb.minm0.7;
   ##以香杨为参考，发生AB区室转换（A-B,B-A,A-A,B-B)的基因表达情况
   ######得到发生AB区室转换的相应.bed文件
   cp /w/00/g/g00/wjl/hic/compartment/ab.change.homologous.bin/xiangy."$i"/40kb.minm0.7/xiangy."$i".homologous.bin.rmAB.AB.csv ./;
   awk '{if($4=="A" && $8=="B" && $1==$5) print $1"\t"$2"\t"$3"\t""xiangy""\t""A-B"}' xiangy."$i".homologous.bin.rmAB.AB.csv|sort -k1,1V -k2,2n -k3,3n > A_B_Axiangy.homologous.bin;
   awk '{if($4=="A" && $8=="B" && $1==$5) print $5"\t"$6"\t"$7"\t""'$i'""\t""A-B"}' xiangy."$i".homologous.bin.rmAB.AB.csv|sort -k1,1V -k2,2n -k3,3n > A_B_B"$i".homologous.bin;
   awk '{if($4=="B" && $8=="A" && $1==$5) print $1"\t"$2"\t"$3"\t""xiangy""\t""B-A"}' xiangy."$i".homologous.bin.rmAB.AB.csv|sort -k1,1V -k2,2n -k3,3n > B_A_Bxiangy.homologous.bin;
   awk '{if($4=="B" && $8=="A" && $1==$5) print $5"\t"$6"\t"$7"\t""'$i'""\t""B-A"}' xiangy."$i".homologous.bin.rmAB.AB.csv|sort -k1,1V -k2,2n -k3,3n > B_A_A"$i".homologous.bin;
   awk '{if($4=="A" && $8=="A" && $1==$5) print $1"\t"$2"\t"$3"\t""xiangy""\t""A-A"}' xiangy."$i".homologous.bin.rmAB.AB.csv|sort -k1,1V -k2,2n -k3,3n > A_A_Axiangy.homologous.bin;
   awk '{if($4=="A" && $8=="A" && $1==$5) print $5"\t"$6"\t"$7"\t""'$i'""\t""A-A"}' xiangy."$i".homologous.bin.rmAB.AB.csv|sort -k1,1V -k2,2n -k3,3n > A_A_A"$i".homologous.bin;
   awk '{if($4=="B" && $8=="B" && $1==$5) print $1"\t"$2"\t"$3"\t""xiangy""\t""B-B"}' xiangy."$i".homologous.bin.rmAB.AB.csv|sort -k1,1V -k2,2n -k3,3n > B_B_Bxiangy.homologous.bin;
   awk '{if($4=="B" && $8=="B" && $1==$5) print $5"\t"$6"\t"$7"\t""'$i'""\t""B-B"}' xiangy."$i".homologous.bin.rmAB.AB.csv|sort -k1,1V -k2,2n -k3,3n > B_B_B"$i".homologous.bin;
   bedtools merge -i A_B_Axiangy.homologous.bin |awk '{print $0"\t""xiangy""\t""A-B"}' > A_B_Axiangy.homologous.bin.merge;
   bedtools merge -i B_A_Bxiangy.homologous.bin |awk '{print $0"\t""xiangy""\t""B-A"}' > B_A_Bxiangy.homologous.bin.merge;
   bedtools merge -i A_A_Axiangy.homologous.bin |awk '{print $0"\t""xiangy""\t""A-A"}' > A_A_Axiangy.homologous.bin.merge;
   bedtools merge -i B_B_Bxiangy.homologous.bin |awk '{print $0"\t""xiangy""\t""B-B"}' > B_B_Bxiangy.homologous.bin.merge;
   bedtools merge -i A_B_B"$i".homologous.bin |awk '{print $0"\t""'$i'""\t""A-B"}' > A_B_B"$i".homologous.bin.merge;
   bedtools merge -i B_A_A"$i".homologous.bin |awk '{print $0"\t""'$i'""\t""B-A"}' > B_A_A"$i".homologous.bin.merge;
   bedtools merge -i A_A_A"$i".homologous.bin |awk '{print $0"\t""'$i'""\t""A-A"}' > A_A_A"$i".homologous.bin.merge;
   bedtools merge -i B_B_B"$i".homologous.bin |awk '{print $0"\t""'$i'""\t""B-B"}' > B_B_B"$i".homologous.bin.merge;
   #####gene expression
   cp /w/00/g/g00/wjl/hic/tad.analysis/"$i"/gene.expression/"$i".gene.fpkm.m.sort.bed ../;
   cp /w/00/g/g00/wjl/hic/tad.analysis/xiangy/gene.expression/xiangy.gene.fpkm.m.sort.bed ../;
   bedtools intersect -nonamecheck -a ../xiangy.gene.fpkm.m.sort.bed -b A_B_Axiangy.homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > A_B_Axiangy.homologous.bin.merge.fpkm;
   bedtools intersect -nonamecheck -a ../xiangy.gene.fpkm.m.sort.bed -b B_A_Bxiangy.homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > B_A_Bxiangy.homologous.bin.merge.fpkm;
   bedtools intersect -nonamecheck -a ../xiangy.gene.fpkm.m.sort.bed -b A_A_Axiangy.homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > A_A_Axiangy.homologous.bin.merge.fpkm;
   bedtools intersect -nonamecheck -a ../xiangy.gene.fpkm.m.sort.bed -b B_B_Bxiangy.homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > B_B_Bxiangy.homologous.bin.merge.fpkm;
   cat A_B_Axiangy.homologous.bin.merge.fpkm B_A_Bxiangy.homologous.bin.merge.fpkm A_A_Axiangy.homologous.bin.merge.fpkm B_B_Bxiangy.homologous.bin.merge.fpkm > xiangy.homologous.bin.fpkm.csv;
   bedtools intersect -nonamecheck -a ../"$i".gene.fpkm.m.sort.bed -b A_B_B"$i".homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > A_B_B"$i".homologous.bin.merge.fpkm;
   bedtools intersect -nonamecheck -a ../"$i".gene.fpkm.m.sort.bed -b B_A_A"$i".homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > B_A_A"$i".homologous.bin.merge.fpkm;
   bedtools intersect -nonamecheck -a ../"$i".gene.fpkm.m.sort.bed -b A_A_A"$i".homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > A_A_A"$i".homologous.bin.merge.fpkm;
   bedtools intersect -nonamecheck -a ../"$i".gene.fpkm.m.sort.bed -b B_B_B"$i".homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > B_B_B"$i".homologous.bin.merge.fpkm;
   cat A_B_B"$i".homologous.bin.merge.fpkm B_A_A"$i".homologous.bin.merge.fpkm A_A_A"$i".homologous.bin.merge.fpkm B_B_B"$i".homologous.bin.merge.fpkm > "$i".homologous.bin.fpkm.csv;
   cat xiangy.homologous.bin.fpkm.csv "$i".homologous.bin.fpkm.csv > xiangy."$i".homologous.bin.fpkm.csv;
done