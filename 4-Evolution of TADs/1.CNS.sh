#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim;
do
  cd /w/00/g/g00/hic/tad.analysis/changed_tad.CNS_enrichment/;
  mkdir Pkor."$i"; cd Pkor."$i";
  ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/Pkor.high.con.tad.bed;
  ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/Pkor.con.tad.bed;
  ln -s /w/00/g/g00/hic/tad.analysis/changed_tad/"$i"/Pkor.spe.tad.bed;
  ln -s /w/00/g/g00/hic/tad.analysis/tad.CNS_enrichment/CNS_final.bed ./;
  ln -s /w/00/g/g00/hic/tad.analysis/tad.CNS_enrichment/CNS_random.bed ./;
  ## Obtain upstream and downstream 50kb positional information for highly conserved, conserved, and diverged TADs
  less Pkor.high.con.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > Pkor.high.con.tadup.50kb.bed;
  less Pkor.high.con.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > Pkor.high.con.taddown.50kb.bed;
  less Pkor.con.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > Pkor.con.tadup.50kb.bed;
  less Pkor.con.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > Pkor.con.taddown.50kb.bed;
  less Pkor.spe.tad.bed | awk '{print $1"\t"$2-50000"\t"$2}' > Pkor.spe.tadup.50kb.bed;
  less Pkor.spe.tad.bed | awk '{print $1"\t"$3"\t"$3+50000}' > Pkor.spe.taddown.50kb.bed;   
  ## Divide highly conserved, conserved, and diverged TADs and their upstream/downstream 50kb regions into 20 bins
  bedtools makewindows -b Pkor.high.con.tad.bed -n 20 -i winnum > Pkor.high.con.tad.20bin.bed;
  bedtools makewindows -b Pkor.high.con.tadup.50kb.bed -n 20 -i winnum > Pkor.high.con.tadup.50kb.20bin.bed;
  bedtools makewindows -b Pkor.high.con.taddown.50kb.bed -n 20 -i winnum > Pkor.high.con.taddown.50kb.20bin.bed;
  bedtools makewindows -b Pkor.con.tad.bed -n 20 -i winnum > Pkor.con.tad.20bin.bed;
  bedtools makewindows -b Pkor.con.tadup.50kb.bed -n 20 -i winnum > Pkor.con.tadup.50kb.20bin.bed;
  bedtools makewindows -b Pkor.con.taddown.50kb.bed -n 20 -i winnum > Pkor.con.taddown.50kb.20bin.bed;
  bedtools makewindows -b Pkor.spe.tad.bed -n 20 -i winnum > Pkor.spe.tad.20bin.bed;
  bedtools makewindows -b Pkor.spe.tadup.50kb.bed -n 20 -i winnum > Pkor.spe.tadup.50kb.20bin.bed;
  bedtools makewindows -b Pkor.spe.taddown.50kb.bed -n 20 -i winnum > Pkor.spe.taddown.50kb.20bin.bed;
  for j in high.con con spe;
  do
    ######observed CNS
    bedtools intersect -nonamecheck -a CNS_final.bed -b Pkor."$j".tadup.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > Pkor."$j".tadup.50kb.20bin.cns;
    bedtools intersect -nonamecheck -a CNS_final.bed -b Pkor."$j".tad.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > Pkor."$j".tad.20bin.cns;
    bedtools intersect -nonamecheck -a CNS_final.bed -b Pkor."$j".taddown.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > Pkor."$j".taddown.50kb.20bin.cns;
    cat Pkor."$j".tadup.50kb.20bin.cns Pkor."$j".tad.20bin.cns Pkor."$j".taddown.50kb.20bin.cns > Pkor."$j".uptaddown.50kb.20bin.cns;
    ######random CNS
    bedtools intersect -nonamecheck -a CNS_random.bed -b Pkor."$j".tadup.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""up"}}'| sort -k2,2n > Pkor."$j".tadup.50kb.20bin.random.cns;
    bedtools intersect -nonamecheck -a CNS_random.bed -b Pkor."$j".tad.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""tad"}}'| sort -k2,2n > Pkor."$j".tad.20bin.random.cns;
    bedtools intersect -nonamecheck -a CNS_random.bed -b Pkor."$j".taddown.50kb.20bin.bed -wa -wb |awk '{if($2>$6 && $3<$7) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | awk '{a[$5]+=1}END{for(i in a){print a[i]"\t"i"\t""down"}}'| sort -k2,2n > Pkor."$j".taddown.50kb.20bin.random.cns;
    cat Pkor."$j".tadup.50kb.20bin.random.cns Pkor."$j".tad.20bin.random.cns Pkor."$j".taddown.50kb.20bin.random.cns > Pkor."$j".uptaddown.50kb.20bin.random.cns;
  done
  ## Calculate observed/random enrichment ratio
  paste Pkor.high.con.uptaddown.50kb.20bin.cns Pkor.high.con.uptaddown.50kb.20bin.random.cns | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3"\t""H-conserved"}'  > Pkor.high.con.uptaddown.50kb.20bin.observed.random.cns.csv;
  paste Pkor.con.uptaddown.50kb.20bin.cns Pkor.con.uptaddown.50kb.20bin.random.cns | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3"\t""Conserved"}'  > Pkor.con.uptaddown.50kb.20bin.observed.random.cns.csv;
  paste Pkor.spe.uptaddown.50kb.20bin.cns Pkor.spe.uptaddown.50kb.20bin.random.cns | awk '{print $1"\t"$4"\t"$1/$4"\t"$2"\t"NR"\t"$3"\t""Diverged"}'  > Pkor.spe.uptaddown.50kb.20bin.observed.random.cns.csv;
  cat Pkor.high.con.uptaddown.50kb.20bin.observed.random.cns.csv Pkor.con.uptaddown.50kb.20bin.observed.random.cns.csv Pkor.spe.uptaddown.50kb.20bin.observed.random.cns.csv > Pkor."$i".changed_tad.cns;
done

## Calculate overall CNS distribution patterns around highly conserved, conserved, and diverged TADs.
cd ..;
cat Pkor.Ppse/Pkor.Ppse.changed_tad.cns Pkor.Pwua/Pkor.Pwua.changed_tad.cns Pkor.Psze/Pkor.Psze.changed_tad.cns Pkor.Plas/Pkor.Plas.changed_tad.cns Pkor.Pyun/Pkor.Pyun.changed_tad.cns Pkor.Prot/Pkor.Prot.changed_tad.cns Pkor.Pqio/Pkor.Pqio.changed_tad.cns Pkor.Pdav/Pkor.Pdav.changed_tad.cns Pkor.Pade/Pkor.Pade.changed_tad.cns Pkor.Psim/Pkor.Psim.changed_tad.cns > all.changed_tad.flank50kb.20bin.cns;
less all.changed_tad.flank50kb.20bin.cns | awk '{if($7=="H-conserved") print $0}' | awk '{a[$5]+=$3;b[$5]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t"$7}}' | sort -k2,2n > all.hcon.flank50kb.20bin.cns;
less all.changed_tad.flank50kb.20bin.cns | awk '{if($7=="Conserved") print $0}' | awk '{a[$5]+=$3;b[$5]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t"$7}}' | sort -k2,2n > all.con.flank50kb.20bin.cns;
less all.changed_tad.flank50kb.20bin.cns | awk '{if($7=="Diverged") print $0}' | awk '{a[$5]+=$3;b[$5]+=1}END{for(i in a){print a[i]/b[i]"\t"i"\t"$7}}' | sort -k2,2n > all.spe.flank50kb.20bin.cns;
cat all.hcon.flank50kb.20bin.cns all.con.flank50kb.20bin.cns all.spe.flank50kb.20bin.cns > all.changed_tad.flank50kb.20bin.cns.enrichment.csv
  
  
