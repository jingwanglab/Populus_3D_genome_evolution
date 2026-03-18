#!/usr/bin/bash
# 1. Identification of homologous bins
# 1.1 Genome alignment using LAST
ref_genome=$1
query_genome=$2
basename1=`basename $genome1`
ref=${basename1%%.*}
basename2=`basename $genome2`
query=${basename2%%.*}

lastdb -P20 -uNEAR -R01 $ref $ref_genome
last-train -P20 --revsym --matsym --gapsym -E0.05 -C2 $ref $query_genome  > $ref-$query.mat
lastal -P20 -m50 -E0.05 -C2 -p $ref-$query.mat  $ref  $query_genome |last-split -m1 > $ref-$query.maf
maf-swap $ref-$query.maf | awk -v idx="$ref" -v idx2="$query" '/^s/ {$2 = (++s % 2 ? idx"." : idx2".") $2}1' | last-split -m1 | maf-swap > $ref-$query-2.maf


# 1.2 Convert maf file to chain file
mafToPsl $ref $query $ref-$query-2.maf output.psl
axtChain -linearGap=medium -psl output.psl -faQ $query -faT $ref $(ref)_$(query).chain

# 1.3 This step generates 40kb solution bed file
bedtools makewindows -g genomesize -w 40000 -i srcwinnum > $ref.bed

# 1.4 This step generates mapping file at 40kb resolution
liftOver $ref.bed -s $(ref)_$(query).chain map.loc unmap.bed -minMatch=0.7

# 1.5 Convert mapping bed file from bin-bed to bin-bin. Input file is mapping bed file, Output file is bin2bin file
python runbedtobin.py -i map.loc -s 40000 -o outputfile

# 1.6 Obtain positional information of homologous bins
awk  -F'\t' 'NR==FNR{a[$4]=$1"\t"$2"\t"$3;next}NR>FNR{if($4 in a) print $1"\t"$2"\t"$3"\t"a[$4]}' spe1.40000.bed map.40000.minm0.7_loc | awk '{if($1==$4) print $0}' > spe2.spe1.40kb.minm0.7_homologous.bin

# 2. Conversion status of A/B compartments for homologous bins
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim;
do
   cd /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/40kb.minm0.7;
   mv "$i".Pkor.40kb.minm0.7_homologous.bin "$i".Pkor.homologous.bin;
   less "$i".Pkor.homologous.bin | awk '{print $1"\t"$2"\t"$3}' > "$i".homologous.bin;
   less "$i".Pkor.homologous.bin | awk '{print $4"\t"$5"\t"$6}' > Pkor.homologous.bin;
   ## Identify homologous bins present in both A/B compartments and remove those with conflicting A/B assignments
   cp /w/00/g/g00/hic/compartment/gene.expression/Pkor/Pkor.ab.sort.bed ../Pkor.ab.merge.bed;
   cp /w/00/g/g00/hic/compartment/gene.expression/"$i"/"$i".ab.sort.bed ../"$i".ab.merge.bed;
   bedtools intersect -nonamecheck -a Pkor.homologous.bin -b ../Pkor.ab.merge.bed -wa -wb |awk '{if($2>=$5 && $3<=$6) print $1","$2","$3}' > 1.txt;
   bedtools intersect -nonamecheck -a "$i".homologous.bin -b ../"$i".ab.merge.bed -wa -wb |awk '{if($2>=$5 && $3<=$6) print $1","$2","$3}' > 2.txt;
   paste 2.txt 1.txt > 3."$i".Pkor.txt;
   less "$i".Pkor.homologous.bin | awk '{print $1","$2","$3"\t"$4","$5","$6}'  > "$i".Pkor.homologous.bin.1
   awk  -F'\t' 'NR==FNR{a[$1]=$1;b[$2]=$2;next}NR>FNR{if($1 in a && $2 in b) print $0}' 3."$i".Pkor.txt "$i".Pkor.homologous.bin.1 > "$i".Pkor.homologous.bin.rmAB
   sed -i 's/,/\t/g' "$i".Pkor.homologous.bin.rmAB
   ## Analyze whether homologous bins are in A or B compartments
   awk '{print $1"\t"$2"\t"$3}' "$i".Pkor.homologous.bin.rmAB> "$i".homologous.bin.rmAB;
   awk '{print $4"\t"$5"\t"$6}' "$i".Pkor.homologous.bin.rmAB> Pkor.homologous.bin.rmAB;
   bedtools intersect -nonamecheck -a Pkor.homologous.bin.rmAB -b ../Pkor.ab.merge.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$7}'> Pkor.homologous.bin.rmAB.AB;
   bedtools intersect -nonamecheck -a "$i".homologous.bin.rmAB -b ../"$i".ab.merge.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$7}'> "$i".homologous.bin.rmAB.AB;
   paste Pkor.homologous.bin.rmAB.AB "$i".homologous.bin.rmAB.AB > Pkor."$i".homologous.bin.AB.csv; 
   ## Count the number of A/B compartment conversions
   awk '{if($4=="A" && $8=="B" && $1==$5) print $0}' Pkor."$i".homologous.bin.AB.csv |wc -l >>one.bin.result.txt;# Count number of A to B 
   awk '{if($4=="B" && $8=="A" && $1==$5) print $0}' Pkor."$i".homologous.bin.AB.csv |wc -l >>one.bin.result.txt;# Count number of B to A
   awk '{if($4=="A" && $8=="A" && $1==$5) print $0}' Pkor."$i".homologous.bin.AB.csv |wc -l >>one.bin.result.txt;# Count number of A compartments without changed
   awk '{if($4=="B" && $8=="B" && $1==$5) print $0}' Pkor."$i".homologous.bin.AB.csv |wc -l >>one.bin.result.txt;# Count number of B compartments without changed

done

