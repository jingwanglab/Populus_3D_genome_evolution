########以香杨为参考，保守边界、保守interior和保守TAD的ka、ks、ka/ks比较
########
1.计算kaks-wgdi（/usr_storage/wjl/ka.ks）
1.1准备好cds和pep文件，共线性分析输出文件(支持MCScanX的共线性分析结果)
1.2修改配置文件
vi total.conf
[ks]
cds_file = pkor_ppse_cds_longest.fa
pep_file = pkor_ppse_pep_longest.fa
align_software = muscle
pairs_file = changxuy_xiangy.chr.uniq.collinearity.txt
ks_file = xiangy.changxuy.ks
1.3 计算kaks
wgdi -ks total.conf
1.4 将输出的结果进行汇总（cat），用于下一步分析（这里只用了ng86的结果）
结果文件：all.xiangy.ka.ks.csv
2.保守边界、保守interior和保守TAD的ka、ks、ka/ks比较
vi xiangy.high.con.spe.boundary.interior.tad.ka.ks.sh
#!/usr/bin/bash
for i in boundary interior tad;
do
   cd /w/00/g/g00/wjl/hic/tad.analysis/conserved.tad.ka.ks;
   #####分别与ka.ks文件取交集
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($2 in a) print $1"\t"$2"\t""H-conserved""\t""'$i'""\t"$3"\t"$4"\t"$5}' all.xiangy.high.con."$i".geneid all.xiangy.ka.ks.csv  > all.xiangy.high.con."$i".ka.ks.csv;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($2 in a) print $1"\t"$2"\t""Conserved""\t""'$i'""\t"$3"\t"$4"\t"$5}' all.xiangy.con."$i".geneid all.xiangy.ka.ks.csv  > all.xiangy.con."$i".ka.ks.csv;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($2 in a) print $1"\t"$2"\t""Specific""\t""'$i'""\t"$3"\t"$4"\t"$5}' all.xiangy.spe."$i".geneid all.xiangy.ka.ks.csv > all.xiangy.spe."$i".ka.ks.csv;
   cat all.xiangy.high.con."$i".ka.ks.csv all.xiangy.con."$i".ka.ks.csv all.xiangy.spe."$i".ka.ks.csv > all.xiangy.high.con.spe."$i".ka.ks.csv;
done
3.将boundary interior tad进行汇总用于R画图
cat all.xiangy.high.con.spe.boundary.ka.ks.csv all.xiangy.high.con.spe.interior.ka.ks.csv all.xiangy.high.con.spe.tad.ka.ks.csv > all.xiangy.high.con.spe.boundary.interior.tad.ka.ks.csv













