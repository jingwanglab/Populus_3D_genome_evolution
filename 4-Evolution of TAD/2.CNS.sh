#####conserved no-coding sequence(cns)在高度、一般保守和特异TAD和TAD边界的分布情况-香杨为参考

vi high.con.spe.tad.boundary.cns.enrichment.sh
#!/usr/bin/bash
export PATH="/data/apps/bedtools/bedtools-2.26.0/bin:$PATH"
for i in changxuy changyey chuany dayey diany qingxiy qiongdaoy shany xiangyey xiaoyey;
do
  cd /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.cns.enrichment/;
  mkdir xiangy."$i"; cd xiangy."$i";
  ###准备数据
  ##bedtools shuffle -i ../xiangy.CNS_final.bed -g /w/00/g/g00/wjl/hic/hicpro/xiangy/new_out2/xiangy.size -chrom |sort -k1,1V -k2,2n -k3,3n > ../xiangy.CNS_random.bed;
  ln -s /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.gene.expression/xiangy."$i".chr/xiangy.high.con.tad.bed;
  ln -s /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.gene.expression/xiangy."$i".chr/xiangy.con.tad.bed;
  ln -s /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.gene.expression/xiangy."$i".chr/xiangy.spe.tad.bed;
  #######获得high conserved、conserved、specific tad的上下游50kb的位置信息
  less xiangy.high.con.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > xiangy.high.con.tadup.50kb.bed;
  less xiangy.high.con.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > xiangy.high.con.taddown.50kb.bed;
  less xiangy.con.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > xiangy.con.tadup.50kb.bed;
  less xiangy.con.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > xiangy.con.taddown.50kb.bed;
  less xiangy.spe.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > xiangy.spe.tadup.50kb.bed;
  less xiangy.spe.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > xiangy.spe.taddown.50kb.bed;   
  #######将high conserved、conserved、specific tad以及其上下游50kb的位置划分为20个bin
  bedtools makewindows -b xiangy.high.con.tad.bed -n 20 -i winnum > xiangy.high.con.tad.20bin.bed;
  bedtools makewindows -b xiangy.high.con.tadup.50kb.bed -n 20 -i winnum > xiangy.high.con.tadup.50kb.20bin.bed;
  bedtools makewindows -b xiangy.high.con.taddown.50kb.bed -n 20 -i winnum > xiangy.high.con.taddown.50kb.20bin.bed;
  bedtools makewindows -b xiangy.con.tad.bed -n 20 -i winnum > xiangy.con.tad.20bin.bed;
  bedtools makewindows -b xiangy.con.tadup.50kb.bed -n 20 -i winnum > xiangy.con.tadup.50kb.20bin.bed;
  bedtools makewindows -b xiangy.con.taddown.50kb.bed -n 20 -i winnum > xiangy.con.taddown.50kb.20bin.bed;
  bedtools makewindows -b xiangy.spe.tad.bed -n 20 -i winnum > xiangy.spe.tad.20bin.bed;
  bedtools makewindows -b xiangy.spe.tadup.50kb.bed -n 20 -i winnum > xiangy.spe.tadup.50kb.20bin.bed;
  bedtools makewindows -b xiangy.spe.taddown.50kb.bed -n 20 -i winnum > xiangy.spe.taddown.50kb.20bin.bed;
  for j in high.con con spe;
  do
    ######observed
    bedtools intersect -nonamecheck -a ../xiangy.CNS_final.bed -b xiangy."$j".tadup.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > xiangy."$j".tadup.50kb.20bin.cns;
    bedtools intersect -nonamecheck -a ../xiangy.CNS_final.bed -b xiangy."$j".tad.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > xiangy."$j".tad.20bin.cns;
    bedtools intersect -nonamecheck -a ../xiangy.CNS_final.bed -b xiangy."$j".taddown.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > xiangy."$j".taddown.50kb.20bin.cns;
    cat xiangy."$j".tadup.50kb.20bin.cns xiangy."$j".tad.20bin.cns xiangy."$j".taddown.50kb.20bin.cns > xiangy."$j".uptaddown.50kb.20bin.cns;
    ######random
    bedtools intersect -nonamecheck -a ../xiangy.CNS_random.bed -b xiangy."$j".tadup.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > xiangy."$j".tadup.50kb.20bin.random.cns;
    bedtools intersect -nonamecheck -a ../xiangy.CNS_random.bed -b xiangy."$j".tad.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > xiangy."$j".tad.20bin.random.cns;
    bedtools intersect -nonamecheck -a ../xiangy.CNS_random.bed -b xiangy."$j".taddown.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > xiangy."$j".taddown.50kb.20bin.random.cns;
    cat xiangy."$j".tadup.50kb.20bin.random.cns xiangy."$j".tad.20bin.random.cns xiangy."$j".taddown.50kb.20bin.random.cns > xiangy."$j".uptaddown.50kb.20bin.random.cns;
  done
  ###observed/random
  paste xiangy.high.con.uptaddown.50kb.20bin.cns xiangy.high.con.uptaddown.50kb.20bin.random.cns | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3"\t""H-conserved"}'  > xiangy.high.con.uptaddown.50kb.20bin.observed.random.cns.csv;
  paste xiangy.con.uptaddown.50kb.20bin.cns xiangy.con.uptaddown.50kb.20bin.random.cns | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3"\t""Conserved"}'  > xiangy.con.uptaddown.50kb.20bin.observed.random.cns.csv;
  paste xiangy.spe.uptaddown.50kb.20bin.cns xiangy.spe.uptaddown.50kb.20bin.random.cns | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3"\t""Specific"}'  > xiangy.spe.uptaddown.50kb.20bin.observed.random.cns.csv;
  cat xiangy.high.con.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.con.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.spe.uptaddown.50kb.20bin.observed.random.cns.csv > xiangy."$i".high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv;
done


####计算整体的cns在高度、一般保守和特异TAD和TAD边界的分布情况
vi all.high.con.spe.tad.boundary.average.cns.enrichment.sh

#!/usr/bin/bash
cat xiangy.changxuy/xiangy.changxuy.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.changyey/xiangy.changyey.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.chuany/xiangy.chuany.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.dayey/xiangy.dayey.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.diany/xiangy.diany.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.qingxiy/xiangy.qingxiy.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.qiongdaoy/xiangy.qiongdaoy.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.shany/xiangy.shany.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.xiangyey/xiangy.xiangyey.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv xiangy.xiaoyey/xiangy.xiaoyey.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv > all.xiangy.others.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv;
less all.xiangy.others.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv | awk '{if($7=="H-conserved") print $0}' | awk '{a[$5]+=$3;b[$5]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t"$7}}' | sort -k2,2n > all.xiangy.others.high.con.uptaddown.50kb.20bin.observed.random.cns.csv;
less all.xiangy.others.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv | awk '{if($7=="Conserved") print $0}' | awk '{a[$5]+=$3;b[$5]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t"$7}}' | sort -k2,2n > all.xiangy.others.con.uptaddown.50kb.20bin.observed.random.cns.csv;
less all.xiangy.others.high.con.spe.uptaddown.50kb.20bin.observed.random.cns.csv | awk '{if($7=="Specific") print $0}' | awk '{a[$5]+=$3;b[$5]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t"$7}}' | sort -k2,2n > all.xiangy.others.spe.uptaddown.50kb.20bin.observed.random.cns.csv;
cat all.xiangy.others.high.con.uptaddown.50kb.20bin.observed.random.cns.csv all.xiangy.others.con.uptaddown.50kb.20bin.observed.random.cns.csv all.xiangy.others.spe.uptaddown.50kb.20bin.observed.random.cns.csv > all.xiangy.others.high.con.spe.uptaddown.50kb.20bin.cns.enrichment.csv
  
  
