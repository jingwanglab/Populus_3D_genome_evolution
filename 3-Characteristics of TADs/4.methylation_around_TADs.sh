#!/usr/bin/bash
for i in Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pkor Pade Psim;
do
   cd /w/00/g/g00/hic/tad.analysis/"$i"/;
   mkdir methylation.tad.20bin; cd methylation.tad.20bin;
   # Create soft links to TAD bed files
   ln -s /w/00/g/g00/hic/tad.analysis/"$i"/tad.bed/"$i".tad.raw.bed ./;
   ln -s /w/00/g/g00/hic/tad.analysis/"$i"/tad.bed/"$i".tadup.50kb.bed ./;
   ln -s /w/00/g/g00/hic/tad.analysis/"$i"/tad.bed/"$i".taddown.50kb.bed ./;
   # Generate 20-bin windows for TAD and flanking regions
   bedtools makewindows -b ../methylation/"$i".tad.raw.bed -n 20 -i winnum > "$i".tad.20bin.bed;
   bedtools makewindows -b "$i".tadup.50kb.bed -n 20 -i winnum > "$i".tadup.50kb.20bin.bed;
   bedtools makewindows -b "$i".taddown.50kb.bed -n 20 -i winnum > "$i".taddown.50kb.20bin.bed;
   # Extract methylation values for each methylation type (CHG, CHH, CpG)
   ####CHG
   bedtools intersect -nonamecheck -a ../methylation/CHG.bed -b "$i".tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".tad.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ../methylation/CHG.bed -b "$i".tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".tadup.50kb.20bin.chg.bed;
   bedtools intersect -nonamecheck -a ../methylation/CHG.bed -b "$i".taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".taddown.50kb.20bin.chg.bed;
   ####CHH
   bedtools intersect -nonamecheck -a ../methylation/CHH.bed -b "$i".tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".tad.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ../methylation/CHH.bed -b "$i".tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".tadup.50kb.20bin.chh.bed;
   bedtools intersect -nonamecheck -a ../methylation/CHH.bed -b "$i".taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".taddown.50kb.20bin.chh.bed;
   ####CpG
   bedtools intersect -nonamecheck -a ../methylation/CpG.bed -b "$i".tad.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".tad.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ../methylation/CpG.bed -b "$i".tadup.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".tadup.50kb.20bin.cpg.bed;
   bedtools intersect -nonamecheck -a ../methylation/CpG.bed -b "$i".taddown.50kb.20bin.bed -wa -wb |awk '{print $4"\t"$10}' > "$i".taddown.50kb.20bin.cpg.bed;

   # Calculate average methylation for each sliding window
   for bin in $(seq 1 20)
   do
	   ###CHG
	   less "$i".tadup.50kb.20bin.chg.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.up.chg.txt;
	   less "$i".tadup.50kb.20bin.chg.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.up.chg.txt;
	   less "$i".tad.20bin.chg.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.tad.chg.txt;
	   less "$i".tad.20bin.chg.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.tad.chg.txt;
	   less "$i".taddown.50kb.20bin.chg.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.down.chg.txt;
	   less "$i".taddown.50kb.20bin.chg.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.down.chg.txt;
	   ####CHH
	   less  "$i".tadup.50kb.20bin.chh.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.up.chh.txt;
	   less  "$i".tadup.50kb.20bin.chh.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.up.chh.txt;
	   less  "$i".tad.20bin.chh.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.tad.chh.txt;
	   less  "$i".tad.20bin.chh.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.tad.chh.txt;
	   less  "$i".taddown.50kb.20bin.chh.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.down.chh.txt;
	   less  "$i".taddown.50kb.20bin.chh.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.down.chh.txt;
	   ####CpG
	   less  "$i".tadup.50kb.20bin.cpg.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.up.cpg.txt;
	   less  "$i".tadup.50kb.20bin.cpg.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.up.cpg.txt;
	   less  "$i".tad.20bin.cpg.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.tad.cpg.txt;
	   less  "$i".tad.20bin.cpg.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.tad.cpg.txt;
	   less  "$i".taddown.50kb.20bin.cpg.bed|awk '{if($2=="'$bin'") sum+=$1}END{print '$i'"\t"sum}' >>sum.down.cpg.txt;
	   less  "$i".taddown.50kb.20bin.cpg.bed|awk '{if($2=="'$bin'") print $2}' |wc -l >> wc.down.cpg.txt;
	done
	
   # Calculate final average methylation for each methylation type
   ### CHG final averages
   paste sum.up.chg.txt wc.up.chg.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".up.20bin.chg.final.csv;
   paste sum.tad.chg.txt wc.tad.chg.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".tad.20bin.chg.final.csv;
   paste sum.down.chg.txt wc.down.chg.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".down.20bin.chg.final.csv;
   cat "$i".up.20bin.chg.final.csv "$i".tad.20bin.chg.final.csv "$i".down.20bin.chg.final.csv | awk '{print $2"\t"$3"\t"NR}' > "$i".uptaddown.20bin.chg.final.csv;
   
   ### CHH final averages
   paste sum.up.chh.txt wc.up.chh.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".up.20bin.chh.final.csv;
   paste sum.tad.chh.txt wc.tad.chh.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".tad.20bin.chh.final.csv;
   paste sum.down.chh.txt wc.down.chh.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".down.20bin.chh.final.csv;
   cat "$i".up.20bin.chh.final.csv "$i".tad.20bin.chh.final.csv "$i".down.20bin.chh.final.csv | awk '{print $2"\t"$3"\t"NR}' > "$i".uptaddown.20bin.chh.final.csv;
   
   ### CpG final averages
   paste sum.up.cpg.txt wc.up.cpg.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".up.20bin.cpg.final.csv;
   paste sum.tad.cpg.txt wc.tad.cpg.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".tad.20bin.cpg.final.csv;
   paste sum.down.cpg.txt wc.down.cpg.txt | awk -v sp="$i" '{print $1"\t"$2/$3"\t"sp}' > "$i".down.20bin.cpg.final.csv;
   cat "$i".up.20bin.cpg.final.csv "$i".tad.20bin.cpg.final.csv "$i".down.20bin.cpg.final.csv | awk '{print $2"\t"$3"\t"NR}' > "$i".uptaddown.20bin.cpg.final.csv;
   rm sum.*.txt wc.*.txt
done
