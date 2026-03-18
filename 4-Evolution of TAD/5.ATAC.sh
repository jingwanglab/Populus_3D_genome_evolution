####在high.con.spe.boundary.interior.tad中的ATAC peak与最近基因的距离比较
vi atac.peak2gene_dis.high.con.spe.boundary.interior.tad.sh
#!/usr/bin/bash
export PATH="/usr_storage/software/bedtools2/bin:$PATH"
##conda activate hicexplorer
for i in changxuy changyey chuany dayey diany qingxiy shany xiaoyey;
do
  cd /usr_storage/wjl/atac/atac.tad/"$i";
  python bed_position_reverse.py;
  ##################计算boundary、interior和TAD中的ATAC peak与最近基因的距离
  #######tad
  bedtools intersect -a "$i".high.con.tad.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""H-conserved""\t""tad"}' > "$i".high.con.tad.peak2gene_dis.txt;
  bedtools intersect -a "$i".con.tad.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""Conserved""\t""tad"}' > "$i".con.tad.peak2gene_dis.txt;
  bedtools intersect -a "$i".spe.tad.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""Specific""\t""tad"}' > "$i".spe.tad.peak2gene_dis.txt;
  cat "$i".high.con.tad.peak2gene_dis.txt "$i".con.tad.peak2gene_dis.txt "$i".spe.tad.peak2gene_dis.txt >  "$i".high.con.spe.tad.peak2gene_dis.txt;
  #######boundary
  bedtools intersect -a "$i".high.con.boundary.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""H-conserved""\t""boundary"}' > "$i".high.con.boundary.peak2gene_dis.txt;
  bedtools intersect -a "$i".con.boundary.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""Conserved""\t""boundary"}' > "$i".con.boundary.peak2gene_dis.txt;
  bedtools intersect -a "$i".spe.boundary.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""Specific""\t""boundary"}' > "$i".spe.boundary.peak2gene_dis.txt;
  cat "$i".high.con.boundary.peak2gene_dis.txt "$i".con.boundary.peak2gene_dis.txt "$i".spe.boundary.peak2gene_dis.txt >  "$i".high.con.spe.boundary.peak2gene_dis.txt;
  #######interior
  bedtools intersect -a "$i".high.con.interior.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""H-conserved""\t""interior"}' > "$i".high.con.interior.peak2gene_dis.txt;
  bedtools intersect -a "$i".con.interior.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""Conserved""\t""interior"}' > "$i".con.interior.peak2gene_dis.txt;
  bedtools intersect -a "$i".spe.interior.bed -b peak2gene_dis.m.txt -nonamecheck -wa -wb| awk '{print $0"\t""Specific""\t""interior"}' > "$i".spe.interior.peak2gene_dis.txt;
  cat "$i".high.con.interior.peak2gene_dis.txt "$i".con.interior.peak2gene_dis.txt "$i".spe.interior.peak2gene_dis.txt >  "$i".high.con.spe.interior.peak2gene_dis.txt;
  cat "$i".high.con.spe.tad.peak2gene_dis.txt "$i".high.con.spe.boundary.peak2gene_dis.txt "$i".high.con.spe.interior.peak2gene_dis.txt > "$i".high.con.spe.tad.boundary.interior.peak2gene_dis.csv;
done