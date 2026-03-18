#####conserved no-coding sequence(cns)在TAD和TAD边界的分布情况-香杨为参考
vi cns.enrichment.tad.boundary.sh

#!/usr/bin/bash
#conda activate hicexplorer
export PATH="/data/apps/bedtools/bedtools-2.26.0/bin:$PATH"
for i in xiangy;
do 
  cd /w/00/g/g00/wjl/hic/tad.analysis/"$i"/cns.enrichment;
  ###准备数据文件
  python chr.lg.bed_position_reverse.py;
  ln -s ../gene.density/tadup.50kb.10bin.bed ./;
  ln -s ../gene.density/tad.10bin.bed ./;
  ln -s ../gene.density/taddown.50kb.10bin.bed ./;
  bedtools shuffle -i "$i".CNS_final.lg.bed -g /w/00/g/g00/wjl/hic/hicpro/"$i"/new_out1_LG/"$i".size -chrom > "$i".CNS_random.lg.bed;
  ###计算observed cns在tad和tad边界的富集情况（计算cns的个数）
  bedtools intersect -nonamecheck -a "$i".CNS_final.lg.bed -b tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > tadup.50kb.10bin.cns.csv;
  bedtools intersect -nonamecheck -a "$i".CNS_final.lg.bed -b tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}'| awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > tad.10bin.cns.csv;
  bedtools intersect -nonamecheck -a "$i".CNS_final.lg.bed -b taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > taddown.50kb.10bin.cns.csv;
  cat tadup.50kb.10bin.cns.csv tad.10bin.cns.csv taddown.50kb.10bin.cns.csv > tadup.tad.taddown.50kb.10bin.cns.csv;
  ####计算random cns在tad和tad边界的富集情况
  bedtools intersect -nonamecheck -a "$i".CNS_random.lg.bed -b tadup.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > tadup.50kb.10bin.random.cns.csv;
  bedtools intersect -nonamecheck -a "$i".CNS_random.lg.bed -b tad.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}'| awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > tad.10bin.random.cns.csv;
  bedtools intersect -nonamecheck -a "$i".CNS_random.lg.bed -b taddown.50kb.10bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > taddown.50kb.10bin.random.cns.csv;
  cat tadup.50kb.10bin.random.cns.csv tad.10bin.random.cns.csv taddown.50kb.10bin.random.cns.csv > tadup.tad.taddown.50kb.10bin.random.cns.csv;
  paste tadup.tad.taddown.50kb.10bin.cns.csv tadup.tad.taddown.50kb.10bin.random.cns.csv | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3}' > tadup.tad.taddown.50kb.10bin.observed.random.cns.csv;
done

####20bin
vi cns.enrichment.tad.boundary.20bin.sh
#!/usr/bin/bash
#conda activate hicexplorer
export PATH="/data/apps/bedtools/bedtools-2.26.0/bin:$PATH"
for i in xiangy;
do 
  cd /w/00/g/g00/wjl/hic/tad.analysis/"$i"/cns.enrichment;
  ###准备数据文件
  ln -s ../methylation.uptaddown.20bin/"$i".tadup.50kb.20bin.bed ./tadup.50kb.20bin.bed;
  ln -s ../methylation.uptaddown.20bin/"$i".tad.20bin.bed ./tad.20bin.bed;
  ln -s ../methylation.uptaddown.20bin/"$i".taddown.50kb.20bin.bed ./taddown.50kb.20bin.bed;
  ###计算observed cns在tad和tad边界的富集情况（计算cns的个数）
  bedtools intersect -nonamecheck -a "$i".CNS_final.lg.bed -b tadup.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > tadup.50kb.20bin.cns.csv;
  bedtools intersect -nonamecheck -a "$i".CNS_final.lg.bed -b tad.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}'| awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > tad.20bin.cns.csv;
  bedtools intersect -nonamecheck -a "$i".CNS_final.lg.bed -b taddown.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > taddown.50kb.20bin.cns.csv;
  cat tadup.50kb.20bin.cns.csv tad.20bin.cns.csv taddown.50kb.20bin.cns.csv > tadup.tad.taddown.50kb.20bin.cns.csv;
  ####计算random cns在tad和tad边界的富集情况
  bedtools intersect -nonamecheck -a "$i".CNS_random.lg.bed -b tadup.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > tadup.50kb.20bin.random.cns.csv;
  bedtools intersect -nonamecheck -a "$i".CNS_random.lg.bed -b tad.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}'| awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > tad.20bin.random.cns.csv;
  bedtools intersect -nonamecheck -a "$i".CNS_random.lg.bed -b taddown.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > taddown.50kb.20bin.random.cns.csv;
  cat tadup.50kb.20bin.random.cns.csv tad.20bin.random.cns.csv taddown.50kb.20bin.random.cns.csv > tadup.tad.taddown.50kb.20bin.random.cns.csv;
  paste tadup.tad.taddown.50kb.20bin.cns.csv tadup.tad.taddown.50kb.20bin.random.cns.csv | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3}' > tadup.tad.taddown.50kb.20bin.observed.random.cns.csv;
done



