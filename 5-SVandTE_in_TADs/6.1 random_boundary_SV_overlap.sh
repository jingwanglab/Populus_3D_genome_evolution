#!/bin/bash
species_list="Ppse Pwua Psze Plas Pyun Prot Pqio Pdav Pade Psim"
for species in $species_list; do
    cd "$species"
	## Step 1: Generate random boundary regions (1000 random shuffles)
    mkdir -p bound_random
    count=1    
    while ((count < 1001)); do
        bedtools shuffle \
            -i "boundary.20kb.bed" \
            -g "${species}_genome.bed" \
            -chrom \
            > "./bound_random/random_${count}.bed"       
        ((count++))
    done
    ## Step 2: Create output directories for SV overlap results
	mkdir -p bound_random_sv/indel
	mkdir -p bound_random_sv/inv
	mkdir -p bound_random_sv/trans
	mkdir -p bound_random_sv/dup
	cd ./bound_random
	for file in ./*.bed
	do
		bedtools intersect -a $file -b ../indel.bed -wao > ../bound_random_sv/indel/${file%.*}".overlap"
		bedtools intersect -a $file -b ../inv.bed -wao > ../bound_random_sv/inv/${file%.*}".overlap"
		bedtools intersect -a $file -b ../trans.bed -wao > ../bound_random_sv/trans/${file%.*}".overlap"
		bedtools intersect -a $file -b ../dup.bed -wao > ../bound_random_sv/dup/${file%.*}".overlap"
	done
	cd ../..
done

