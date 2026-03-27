#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim;
do
   cd /w/00/g/g00/hic/compartment/ab.change.homologous.bin.gene.expression/;
   cd Pkor."$i";
   ln -s /w/00/g/g00/hic/compartment/gene.expression/"$i"/"$i".gene.tpm.sort.bed ./;
   ln -s /w/00/g/g00/hic/compartment/gene.expression/Pkor/Pkor.gene.tpm.sort.bed ./;
   bedtools intersect -nonamecheck -a Pkor.gene.tpm.sort.bed -b A2B_APkor.homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > A2B_APkor.homologous.bin.merge.tpm;
   bedtools intersect -nonamecheck -a Pkor.gene.tpm.sort.bed -b B2A_BPkor.homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > B2A_BPkor.homologous.bin.merge.tpm;
   bedtools intersect -nonamecheck -a Pkor.gene.tpm.sort.bed -b A2A_APkor.homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > A2A_APkor.homologous.bin.merge.tpm;
   bedtools intersect -nonamecheck -a Pkor.gene.tpm.sort.bed -b B2B_BPkor.homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > B2B_BPkor.homologous.bin.merge.tpm;
   cat A2B_APkor.homologous.bin.merge.tpm B2A_BPkor.homologous.bin.merge.tpm A2A_APkor.homologous.bin.merge.tpm B2B_BPkor.homologous.bin.merge.tpm > Pkor.homologous.bin.tpm.csv;
   bedtools intersect -nonamecheck -a "$i".gene.tpm.sort.bed -b A2B_B"$i".homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > A2B_B"$i".homologous.bin.merge.tpm;
   bedtools intersect -nonamecheck -a "$i".gene.tpm.sort.bed -b B2A_A"$i".homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > B2A_A"$i".homologous.bin.merge.tpm;
   bedtools intersect -nonamecheck -a "$i".gene.tpm.sort.bed -b A2A_A"$i".homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > A2A_A"$i".homologous.bin.merge.tpm;
   bedtools intersect -nonamecheck -a "$i".gene.tpm.sort.bed -b B2B_B"$i".homologous.bin.merge -wa -wb |awk '{if($2>=$6 && $3<=$7) print $4"\t"$8"\t"$9}' > B2B_B"$i".homologous.bin.merge.tpm;
   cat A2B_B"$i".homologous.bin.merge.tpm B2A_A"$i".homologous.bin.merge.tpm A2A_A"$i".homologous.bin.merge.tpm B2B_B"$i".homologous.bin.merge.tpm > "$i".homologous.bin.tpm.csv;
   cat Pkor.homologous.bin.tpm.csv "$i".homologous.bin.tpm.csv > Pkor."$i".homologous.bin.tpm.csv;
done