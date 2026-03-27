#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim;
do
  cd /w/00/g/g00/compartment/ab.change.homologous.bin.core.dis.pri;
  mkdir Pkor."$i";cd Pkor."$i";
  ## Create soft links to geneid files with A/B compartment conversions
  ln -s ../../ab.change.homologous.bin.methylation/xiangy."$i"/40kb.minm0.7/xiangy.homologous.bin.abchange.geneid ./;
  ln -s ../../ab.change.homologous.bin.methylation/xiangy."$i"/40kb.minm0.7/"$i".homologous.bin.abchange.geneid ./;
  #ln -s /w/00/g/g00/hic/compartment/ab.change.core.dis.pri/all.core.dis.pri.geneid ../;
  ## Classify genes with A/B compartment conversions into core, dis, and pri categories
  for j in A-A A-B B-B B-A;
  do
     awk  -F'\t' 'NR==FNR{a[$1]=$3;next}NR>FNR{if($1 in a) print $0"\t"a[$1]}' Pkor.homologous.bin.abchange.geneid ../all.core.dis.pri.geneid |awk '{if($3=="'$j'") print $0}' | awk '{a[$2]+=1}END{for(i in a){print a[i]"\t"i"\t""'$j'"}}' | awk '{print $1"\t"$2"\t"$3"\t""Pkor_'$i'"}'> Pkor."$j".core.dis.pri.csv;
     awk  -F'\t' 'NR==FNR{a[$1]=$3;next}NR>FNR{if($1 in a) print $0"\t"a[$1]}' "$i".homologous.bin.abchange.geneid ../all.core.dis.pri.geneid |awk '{if($3=="'$j'") print $0}' | awk '{a[$2]+=1}END{for(i in a){print a[i]"\t"i"\t""'$j'"}}' | awk '{print $1"\t"$2"\t"$3"\t""'$i'_Pkor"}'> "$i"."$j".core.dis.pri.csv;
  done
  cat Pkor.A-A.core.dis.pri.csv Pkor.A-B.core.dis.pri.csv Pkor.B-B.core.dis.pri.csv Pkor.B-A.core.dis.pri.csv > Pkor.abchange.pangene.csv;
  cat "$i".A-A.core.dis.pri.csv "$i".A-B.core.dis.pri.csv "$i".B-B.core.dis.pri.csv "$i".B-A.core.dis.pri.csv > "$i".abchange.pangene.csv;
done

