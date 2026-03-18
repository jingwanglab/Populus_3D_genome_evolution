######以香杨为参考保守TAD的甲基化水平比较
1.将TAD分为高度保守的（mcscanx+liftoverB）、一般保守的（mcscanx）和特异的
conda activate hicexplorer#####调用该环境里面的python库
export PATH="/data/apps/bedtools/bedtools-2.26.0/bin:$PATH"

1.1 
vi conserved.tad.methylation.sh
#!/usr/bin/bash
for i in changxuy changyey chuany dayey diany qingxiy qiongdaoy shany xiangyey xiaoyey;
do
   cd /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.methylation/xiangy."$i".chr;
   ln -s /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.gene.expression/xiangy."$i".chr/"$i".high.con.tad.bed;
   ln -s /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.gene.expression/xiangy."$i".chr/"$i".con.tad.bed;
   ln -s /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.gene.expression/xiangy."$i".chr/"$i".spe.tad.bed;
   ######调整high conserved、conserved、specific tad的位置，由Chr变为LG
   python chr.lg.high.con.bed_position_reverse.py;
   python chr.lg.con.bed_position_reverse.py;
   python chr.lg.spe.bed_position_reverse.py;
   ######软连接甲基化数据
   ln -s /w/00/g/g00/wjl/hic/tad.analysis/"$i"/methylation.boundary.50bin/CHG.bed ./CHG.bed;
   ln -s /w/00/g/g00/wjl/hic/tad.analysis/"$i"/methylation.boundary.50bin/CHH.bed ./CHH.bed;
   ln -s /w/00/g/g00/wjl/hic/tad.analysis/"$i"/methylation.boundary.50bin/CpG.bed ./CpG.bed;
   #######获得high conserved、conserved、specific tad的上下游50kb的位置信息
   less "$i".high.con.tad.m.bed | awk '{print $1"\t"$2-50000"\t"$2}' > "$i".high.con.tadup.50kb.m.bed;
   less "$i".high.con.tad.m.bed | awk '{print $1"\t"$3"\t"$3+50000}' > "$i".high.con.taddown.50kb.m.bed;
   less "$i".con.tad.m.bed | awk '{print $1"\t"$2-50000"\t"$2}' > "$i".con.tadup.50kb.m.bed;
   less "$i".con.tad.m.bed | awk '{print $1"\t"$3"\t"$3+50000}' > "$i".con.taddown.50kb.m.bed;
   less "$i".spe.tad.m.bed | awk '{print $1"\t"$2-50000"\t"$2}' > "$i".spe.tadup.50kb.m.bed;
   less "$i".spe.tad.m.bed | awk '{print $1"\t"$3"\t"$3+50000}' > "$i".spe.taddown.50kb.m.bed;
   #######将high conserved、conserved、specific tad以及其上下游50kb的位置划分为20个bin
   bedtools makewindows -b "$i".high.con.tad.m.bed -n 20 -i winnum > "$i".high.con.tad.20bin.bed;
   bedtools makewindows -b "$i".high.con.tadup.50kb.m.bed -n 20 -i winnum > "$i".high.con.tadup.50kb.20bin.bed;
   bedtools makewindows -b "$i".high.con.taddown.50kb.m.bed -n 20 -i winnum > "$i".high.con.taddown.50kb.20bin.bed;
   bedtools makewindows -b "$i".con.tad.m.bed -n 20 -i winnum > "$i".con.tad.20bin.bed;
   bedtools makewindows -b "$i".con.tadup.50kb.m.bed -n 20 -i winnum > "$i".con.tadup.50kb.20bin.bed;
   bedtools makewindows -b "$i".con.taddown.50kb.m.bed -n 20 -i winnum > "$i".con.taddown.50kb.20bin.bed;
   bedtools makewindows -b "$i".spe.tad.m.bed -n 20 -i winnum > "$i".spe.tad.20bin.bed;
   bedtools makewindows -b "$i".spe.tadup.50kb.m.bed -n 20 -i winnum > "$i".spe.tadup.50kb.20bin.bed;
   bedtools makewindows -b "$i".spe.taddown.50kb.m.bed -n 20 -i winnum > "$i".spe.taddown.50kb.20bin.bed;
   #######得到high conserved、conserved、specific tad以及其上下游50kb中每个bin的甲基化水平再合并
   ####CHG
   #######high conserved
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".high.con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".high.con.tad.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".high.con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".high.con.tadup.50kb.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".high.con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}'| sort -k2,2n > "$i".high.con.taddown.50kb.20bin.chg.bed;
   cat "$i".high.con.tadup.50kb.20bin.chg.bed "$i".high.con.tad.20bin.chg.bed  "$i".high.con.taddown.50kb.20bin.chg.bed | awk '{print $1"\t"NR"\t"$3"\t""H-conserved"}'> "$i".high.con.uptaddown.20bin.chg.csv;
   #######conserved
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".con.tad.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".con.tadup.50kb.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".con.taddown.50kb.20bin.chg.bed;
   cat "$i".con.tadup.50kb.20bin.chg.bed "$i".con.tad.20bin.chg.bed  "$i".con.taddown.50kb.20bin.chg.bed | awk '{print $1"\t"NR"\t"$3"\t""Conserved"}'> "$i".con.uptaddown.20bin.chg.csv;
   #######specific
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".spe.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".spe.tad.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".spe.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".spe.tadup.50kb.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ./CHG.bed -b "$i".spe.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".spe.taddown.50kb.20bin.chg.bed;
   cat "$i".spe.tadup.50kb.20bin.chg.bed "$i".spe.tad.20bin.chg.bed  "$i".spe.taddown.50kb.20bin.chg.bed | awk '{print $1"\t"NR"\t"$3"\t""Specific"}'> "$i".spe.uptaddown.20bin.chg.csv;
   cat "$i".high.con.uptaddown.20bin.chg.csv "$i".con.uptaddown.20bin.chg.csv "$i".spe.uptaddown.20bin.chg.csv > "$i".high.con.spe.uptaddown.20bin.chg.csv;
   ####CHH
   #######high conserved
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".high.con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".high.con.tad.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".high.con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".high.con.tadup.50kb.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".high.con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".high.con.taddown.50kb.20bin.chh.bed;
   cat "$i".high.con.tadup.50kb.20bin.chh.bed "$i".high.con.tad.20bin.chh.bed  "$i".high.con.taddown.50kb.20bin.chh.bed | awk '{print $1"\t"NR"\t"$3"\t""H-conserved"}'> "$i".high.con.uptaddown.20bin.chh.csv;
   #######conserved
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".con.tad.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".con.tadup.50kb.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".con.taddown.50kb.20bin.chh.bed;
   cat "$i".con.tadup.50kb.20bin.chh.bed "$i".con.tad.20bin.chh.bed  "$i".con.taddown.50kb.20bin.chh.bed | awk '{print $1"\t"NR"\t"$3"\t""Conserved"}'> "$i".con.uptaddown.20bin.chh.csv;
   #######specific
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".spe.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".spe.tad.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".spe.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".spe.tadup.50kb.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ./CHH.bed -b "$i".spe.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}'| sort -k2,2n > "$i".spe.taddown.50kb.20bin.chh.bed;
   cat "$i".spe.tadup.50kb.20bin.chh.bed "$i".spe.tad.20bin.chh.bed  "$i".spe.taddown.50kb.20bin.chh.bed | awk '{print $1"\t"NR"\t"$3"\t""Specific"}'> "$i".spe.uptaddown.20bin.chh.csv;
   cat "$i".high.con.uptaddown.20bin.chh.csv "$i".con.uptaddown.20bin.chh.csv "$i".spe.uptaddown.20bin.chh.csv > "$i".high.con.spe.uptaddown.20bin.chh.csv;
   ####CpG
   #######high conserved
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".high.con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".high.con.tad.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".high.con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".high.con.tadup.50kb.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".high.con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".high.con.taddown.50kb.20bin.cpg.bed;
   cat "$i".high.con.tadup.50kb.20bin.cpg.bed "$i".high.con.tad.20bin.cpg.bed  "$i".high.con.taddown.50kb.20bin.cpg.bed | awk '{print $1"\t"NR"\t"$3"\t""H-conserved"}'> "$i".high.con.uptaddown.20bin.cpg.csv;
   #######conserved
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".con.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".con.tad.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".con.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".con.tadup.50kb.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".con.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}' | sort -k2,2n > "$i".con.taddown.50kb.20bin.cpg.bed;
   cat "$i".con.tadup.50kb.20bin.cpg.bed "$i".con.tad.20bin.cpg.bed  "$i".con.taddown.50kb.20bin.cpg.bed | awk '{print $1"\t"NR"\t"$3"\t""Conserved"}'> "$i".con.uptaddown.20bin.cpg.csv;
   #######specific
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".spe.tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' | awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""tad"}}'| sort -k2,2n > "$i".spe.tad.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".spe.tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""up"}}' | sort -k2,2n > "$i".spe.tadup.50kb.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ./CpG.bed -b "$i".spe.taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}'| awk '{a[$2]+=$1;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""down"}}'| sort -k2,2n > "$i".spe.taddown.50kb.20bin.cpg.bed;
   cat "$i".spe.tadup.50kb.20bin.cpg.bed "$i".spe.tad.20bin.cpg.bed  "$i".spe.taddown.50kb.20bin.cpg.bed | awk '{print $1"\t"NR"\t"$3"\t""Specific"}'> "$i".spe.uptaddown.20bin.cpg.csv;
   cat "$i".high.con.uptaddown.20bin.cpg.csv "$i".con.uptaddown.20bin.cpg.csv "$i".spe.uptaddown.20bin.cpg.csv > "$i".high.con.spe.uptaddown.20bin.cpg.csv;
done
