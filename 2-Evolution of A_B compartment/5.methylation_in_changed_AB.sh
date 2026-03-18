#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim;
do
   cd /w/00/g/g00/hic/compartment/ab.change.homologous.bin.methylation/;
   mkdir Pkor."$i";
   cd Pkor."$i";
   ##########soft link homologous.bin files with A/B compartment changes (A-B, B-A, A-A, B-B)
   ln -s /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/A2B_A.Pkor.homologous.bin.merge ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/B2A_B.Pkor.homologous.bin.merge ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/A2B_B."$i".homologous.bin.merge ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/B2A_A."$i".homologous.bin.merge ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/A2A_A.Pkor.homologous.bin.merge ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/B2B_B.Pkor.homologous.bin.merge ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/B2B_B."$i".homologous.bin.merge ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.homologous.bin/Pkor."$i"/A2A_A."$i".homologous.bin.merge ./;
   ###########get geneids from homologous.bins with A-B, B-A, A-A, B-B transitions
   ln -s /w/00/g/g00/hic/compartment/gene.class/Pkor/geneid.bed ./Pkor.geneid.bed;
   ln -s /w/00/g/g00/hic/compartment/gene.class/"$i"/geneid.bed ./"$i".geneid.bed;
   bedtools intersect -nonamecheck -a A2B_A.Pkor.homologous.bin.merge -b Pkor.geneid.bed -wa -wb |awk '{if($7>=$2 && $8<=$3) print $9"\t"$4"\t"$5}' > A2B_A.Pkor.homologous.bin.geneid;
   bedtools intersect -nonamecheck -a B2A_B.Pkor.homologous.bin.merge -b Pkor.geneid.bed -wa -wb |awk '{if($7>=$2 && $8<=$3) print $9"\t"$4"\t"$5}' > B2A_B.Pkor.homologous.bin.geneid;
   bedtools intersect -nonamecheck -a A2A_A.Pkor.homologous.bin.merge -b Pkor.geneid.bed -wa -wb |awk '{if($7>=$2 && $8<=$3) print $9"\t"$4"\t"$5}' > A2A_A.Pkor.homologous.bin.geneid;
   bedtools intersect -nonamecheck -a B2B_B.Pkor.homologous.bin.merge -b Pkor.geneid.bed -wa -wb |awk '{if($7>=$2 && $8<=$3) print $9"\t"$4"\t"$5}' > B2B_B.Pkor.homologous.bin.geneid;
   bedtools intersect -nonamecheck -a A2B_B."$i".homologous.bin.merge -b "$i".geneid.bed -wa -wb |awk '{if($7>=$2 && $8<=$3) print $9"\t"$4"\t"$5}' > A2B_B."$i".homologous.bin.geneid;
   bedtools intersect -nonamecheck -a B2A_A."$i".homologous.bin.merge -b "$i".geneid.bed -wa -wb |awk '{if($7>=$2 && $8<=$3) print $9"\t"$4"\t"$5}' > B2A_A."$i".homologous.bin.geneid;
   bedtools intersect -nonamecheck -a A2A_A."$i".homologous.bin.merge -b "$i".geneid.bed -wa -wb |awk '{if($7>=$2 && $8<=$3) print $9"\t"$4"\t"$5}' > A2A_A."$i".homologous.bin.geneid;
   bedtools intersect -nonamecheck -a B2B_B."$i".homologous.bin.merge -b "$i".geneid.bed -wa -wb |awk '{if($7>=$2 && $8<=$3) print $9"\t"$4"\t"$5}' > B2B_B."$i".homologous.bin.geneid;
   ##########intersect genes with methylation results files and calculate average methylation level per bin (A-B, B-A)
   ln -s /w/00/g/g00/hic/compartment/ab.change.geneupdown.methylation/Pkor/Pkor.CHG.results ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.geneupdown.methylation/Pkor/Pkor.CHH.results ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.geneupdown.methylation/Pkor/Pkor.CpG.results ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.geneupdown.methylation/"$i"/"$i".CHG.results ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.geneupdown.methylation/"$i"/"$i".CHH.results ./;
   ln -s /w/00/g/g00/hic/compartment/ab.change.geneupdown.methylation/"$i"/"$i".CpG.results ./;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2B_A.Pkor.homologous.bin.geneid ./Pkor.CHG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""A-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2B_A.Pkor.chg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2B_A.Pkor.homologous.bin.geneid ./Pkor.CHH.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""A-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2B_A.Pkor.chh;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2B_A.Pkor.homologous.bin.geneid ./Pkor.CpG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""A-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2B_A.Pkor.cpg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2A_B.Pkor.homologous.bin.geneid ./Pkor.CHG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""B-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2A_B.Pkor.chg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2A_B.Pkor.homologous.bin.geneid ./Pkor.CHH.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""B-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2A_B.Pkor.chh;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2A_B.Pkor.homologous.bin.geneid ./Pkor.CpG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""B-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2A_B.Pkor.cpg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2B_B."$i".homologous.bin.geneid ./"$i".CHG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2B_B."$i".chg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2B_B."$i".homologous.bin.geneid ./"$i".CHH.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2B_B."$i".chh;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2B_B."$i".homologous.bin.geneid ./"$i".CpG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2B_B."$i".cpg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2A_A."$i".homologous.bin.geneid ./"$i".CHG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2A_A."$i".chg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2A_A."$i".homologous.bin.geneid ./"$i".CHH.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2A_A."$i".chh;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2A_A."$i".homologous.bin.geneid ./"$i".CpG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2A_A."$i".cpg;
   ############intersect genes with methylation results files and calculate average methylation level per bin (A-A, B-B)
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2A_A.Pkor.homologous.bin.geneid ./Pkor.CHG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""A-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2A_A.Pkor.chg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2A_A.Pkor.homologous.bin.geneid ./Pkor.CHH.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""A-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2A_A.Pkor.chh;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2A_A.Pkor.homologous.bin.geneid ./Pkor.CpG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""A-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2A_A.Pkor.cpg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2B_B.Pkor.homologous.bin.geneid ./Pkor.CHG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""B-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2B_B.Pkor.chg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2B_B.Pkor.homologous.bin.geneid ./Pkor.CHH.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""B-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2B_B.Pkor.chh;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2B_B.Pkor.homologous.bin.geneid ./Pkor.CpG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""Pkor""\t""B-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2B_B.Pkor.cpg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2B_B."$i".homologous.bin.geneid ./"$i".CHG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2B_B."$i".chg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2B_B."$i".homologous.bin.geneid ./"$i".CHH.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2B_B."$i".chh;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  B2B_B."$i".homologous.bin.geneid ./"$i".CpG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""B-B"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > B2B_B."$i".cpg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2A_A."$i".homologous.bin.geneid ./"$i".CHG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2A_A."$i".chg;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2A_A."$i".homologous.bin.geneid ./"$i".CHH.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2A_A."$i".chh;
   awk  -F'\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a) print $0}'  A2A_A."$i".homologous.bin.geneid ./"$i".CpG.results | awk '{a[$2]+=$3;b[$2]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t""'$i'""\t""A-A"}}' | sed 's/_/\t/g' |sed 's/up/aup/g' | sed 's/body/bbody/g' | sed 's/down/cdown/g'|sort -k2,2 -k3n |awk '{print $1"\t"NR"\t"$4"\t"$5}' > A2A_A."$i".cpg;
   ##########merge files
   cat A2A_A.Pkor.chg A2B_A.Pkor.chg B2B_B.Pkor.chg B2A_B.Pkor.chg  > Pkor.all.chg.csv;
   cat A2A_A.Pkor.chh A2B_A.Pkor.chh B2B_B.Pkor.chh B2A_B.Pkor.chh  > Pkor.all.chh.csv;
   cat A2A_A.Pkor.cpg A2B_A.Pkor.cpg B2B_B.Pkor.cpg B2A_B.Pkor.cpg  > Pkor.all.cpg.csv;
   cat A2A_A."$i".chg B2A_A."$i".chg B2B_B."$i".chg A2B_B."$i".chg  > "$i".all.chg.csv;
   cat A2A_A."$i".chh B2A_A."$i".chh B2B_B."$i".chh A2B_B."$i".chh  > "$i".all.chh.csv;
   cat A2A_A."$i".cpg B2A_A."$i".cpg B2B_B."$i".cpg A2B_B."$i".cpg  > "$i".all.cpg.csv;
done