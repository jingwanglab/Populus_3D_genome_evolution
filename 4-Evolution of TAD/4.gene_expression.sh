######以香杨为参考保守TAD的基因表达水平比较
1.将TAD分为高度保守的（mcscanx+liftoverB）、一般保守的（mcscanx）和特异的
export PATH="/data/apps/bedtools/bedtools-2.26.0/bin:$PATH"
vi xiangy.con.tad.gene.expression.sh
#!/usr/bin/bash
for i in changxuy changyey chuany dayey diany qingxiy qiongdaoy shany xiangyey xiaoyey;
do
   cd /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.gene.expression/xiangy."$i".chr;
   ln -s /w/00/g/g00/wjl/hic/tad.analysis/"$i"/gene.expression/tad.raw.bed ./"$i".tad.raw.bed;
   ln -s /w/00/g/g00/wjl/hic/compartment/gene.class/"$i"/geneid.m.bed ./"$i".geneid.m.bed;
   ln -s /w/00/g/g00/wjl/hic/compartment/ab.change.gene.expression/xiangy."$i".chr/"$i".geneid.fpkm ./;
   ################
   less "$i"_xiangy.boundary.pass | awk '{print $1"\t"$2"\t"$3}' > "$i".high.con.tad.bed;#############得到高度保守的tad
   less /w/00/g/g00/wjl/hic/tad.analysis/conservation/xiangy."$i".chr/"$i".xiangy.tad.con.70.finall.csv |awk '{print $1"\t"$2"\t"$3}' > "$i".mcscanx.con.tad.bed;######得到mcscanx结果的保守tad
   sort "$i".mcscanx.con.tad.bed "$i".high.con.tad.bed "$i".high.con.tad.bed |uniq -u > "$i".con.tad.bed;######得到一般保守的tad
   sort "$i".tad.raw.bed "$i".mcscanx.con.tad.bed "$i".mcscanx.con.tad.bed |uniq -u  > "$i".spe.tad.bed;#######得到不保守的tad
   #####保守和不保守的tad.bed文件分别与基因表达文件取交集
   bedtools intersect -a "$i".geneid.m.bed -b "$i".high.con.tad.bed -wa -wb -nonamecheck | awk '{if($2>$6 && $3<$7) print $4"\t"$5"\t"$6"\t"$7}'  > "$i".high.con.tad.geneid;
   bedtools intersect -a "$i".geneid.m.bed -b "$i".con.tad.bed -wa -wb -nonamecheck | awk '{if($2>$6 && $3<$7) print $4"\t"$5"\t"$6"\t"$7}'  > "$i".con.tad.geneid;
   bedtools intersect -a "$i".geneid.m.bed -b "$i".spe.tad.bed -wa -wb -nonamecheck | awk '{if($2>$6 && $3<$7) print $4"\t"$5"\t"$6"\t"$7}'  > "$i".spe.tad.geneid;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0"\t""H-conserved"}' "$i".high.con.tad.geneid "$i".geneid.fpkm  > "$i".high.con.tad.fpkm;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0"\t""Conserved"}' "$i".con.tad.geneid "$i".geneid.fpkm  > "$i".con.tad.fpkm;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0"\t""Specific"}' "$i".spe.tad.geneid "$i".geneid.fpkm > "$i".spe.tad.fpkm;
   cat "$i".high.con.tad.fpkm "$i".con.tad.fpkm "$i".spe.tad.fpkm > "$i".high.con.spe.tad.fpkm.csv;
done
