#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim;
do
  cd /w/00/g/g00/hic/compartment/ab.change.homologous.bin/"$i";
  awk '{if($4!=$8 && $1==$5 ) print $5"\t"$6"\t"$7"\t""'$i'""\t""changed"}' Pkor."$i".homologous.bin.AB.csv|sort -k1,1V -k2,2n -k3,3n > "$i".changed.homologous.bin;
  awk '{if($4==$8 && $1==$5 ) print $5"\t"$6"\t"$7"\t""'$i'""\t""conserved"}' Pkor."$i".homologous.bin.AB.csv|sort -k1,1V -k2,2n -k3,3n > "$i".conserved.homologous.bin;
  bedtools coverage -a "$i".changed.homologous.bin -b TE.bed -nonamecheck | awk '{print $9"\t""Changed"}' > "$i".changed.homologous.bin.te.coverage;
  bedtools coverage -a "$i".conserved.homologous.bin -b TE.bed -nonamecheck | awk '{print $9"\t""Conserved"}' > "$i".conserved.homologous.bin.te.coverage;
  cat "$i".changed.homologous.bin.te.coverage "$i".conserved.homologous.bin.te.coverage > "$i".changed.conserved.homologous.bin.te.coverage.csv;
  bedtools coverage -a "$i".changed.homologous.bin -b sv.bed -nonamecheck | awk '{print $9"\t""Changed"}' > "$i".changed.homologous.bin.sv.coverage;
  bedtools coverage -a "$i".conserved.homologous.bin -b sv.bed -nonamecheck | awk '{print $9"\t""Conserved"}' > "$i".conserved.homologous.bin.sv.coverage;
  cat "$i".changed.homologous.bin.sv.coverage "$i".conserved.homologous.bin.sv.coverage > "$i".changed.conserved.homologous.bin.sv.coverage.csv;
done
