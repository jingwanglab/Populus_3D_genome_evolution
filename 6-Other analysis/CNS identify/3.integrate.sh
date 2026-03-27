#!/bin/bash

#integrate maf file

for ((i=1; i<=9; i++))
do
	body="Chr0"$i"_body.maf"
	whole="Chr0"$i".maf"
	cat head.maf $body > $whole
done

for ((i=10; i<=19; i++))
do 
	body="Chr"$i"_body.maf"
	whole="Chr"$i".maf"
	cat head.maf $body > $whole
done
