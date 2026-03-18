#####11个杨属物种的core、dispensable、private基因在tad和tad边界的分布情况
1.
export PATH="/data/apps/bedtools/bedtools-2.26.0/bin:$PATH"
cd /w/00/g/g00/wjl/hic/tad.analysis/
vi core.dis.pri.gene.boundarytadothers.20kb.sh
#!/usr/bin/bash
for i in changxuy changyey chuany dayey diany qingxiy qiongdaoy shany xiangy xiangyey xiaoyey;
do
   cd /w/00/g/g00/wjl/hic/tad.analysis/"$i"/;
   mkdir core.dis.pri.gene ;cd core.dis.pri.gene;
   #########得到core、dispensable、private基因的.bed文件(位置信息)
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($4 in a)print $0}'  /w/00/g/g00/wjl/hic/compartment/core.dis.pri.gene/core.gene/"$i"_core_gene.list /w/00/g/g00/wjl/hic/compartment/gene.class/"$i"/geneid.bed > "$i"_core_gene.bed;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($4 in a)print $0}'  /w/00/g/g00/wjl/hic/compartment/core.dis.pri.gene/dispensable.gene/"$i"_dis_gene.list /w/00/g/g00/wjl/hic/compartment/gene.class/"$i"/geneid.bed > "$i"_dis_gene.bed;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($4 in a)print $0}'  /w/00/g/g00/wjl/hic/compartment/core.dis.pri.gene/private.gene/"$i"_pri_gene.list /w/00/g/g00/wjl/hic/compartment/gene.class/"$i"/geneid.bed > "$i"_pri_gene.bed;
   #########将tad和tad边界以及其他区域与core、dispensable、private基因的.bed文件取交集，只保留完全包含在区域中的gene
   cp /w/00/g/g00/wjl/hic/tad.analysis/"$i"/gene.class/boundary.tad.others.20kb.bed ./;
   bedtools intersect -nonamecheck -a "$i"_core_gene.bed -b boundary.tad.others.20kb.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > boundary.tad.others.20kb.core.bed;
   bedtools intersect -nonamecheck -a "$i"_dis_gene.bed -b boundary.tad.others.20kb.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > boundary.tad.others.20kb.dis.bed;
   bedtools intersect -nonamecheck -a "$i"_pri_gene.bed -b boundary.tad.others.20kb.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > boundary.tad.others.20kb.pri.bed;
   ####统计在tad、tad边界以及其他区域中各个基因类型的个数（core、dispensable、private）
   less boundary.tad.others.20kb.core.bed |awk '{if($5=="boundary") print $0}'|wc -l >>1.txt;
   less boundary.tad.others.20kb.core.bed |awk '{if($5=="tad") print $0}'|wc -l >>1.txt;
   less boundary.tad.others.20kb.core.bed |awk '{if($5=="others") print $0}'|wc -l >>1.txt;
   less boundary.tad.others.20kb.dis.bed |awk '{if($5=="boundary") print $0}'|wc -l >>1.txt;
   less boundary.tad.others.20kb.dis.bed |awk '{if($5=="tad") print $0}'|wc -l >>1.txt;
   less boundary.tad.others.20kb.dis.bed |awk '{if($5=="others") print $0}'|wc -l >>1.txt;
   less boundary.tad.others.20kb.pri.bed |awk '{if($5=="boundary") print $0}'|wc -l >>1.txt;
   less boundary.tad.others.20kb.pri.bed |awk '{if($5=="tad") print $0}'|wc -l >>1.txt;
   less boundary.tad.others.20kb.pri.bed |awk '{if($5=="others") print $0}'|wc -l >>1.txt;
done