#!/bin/bash

#extract chromosome by chromosome
gawk 'BEGIN{RS=""} /Pkor.Chr01/ {print $0 "\n" > "Chr01_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr02/ {print $0 "\n" > "Chr02_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr03/ {print $0 "\n" > "Chr03_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr04/ {print $0 "\n" > "Chr04_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr05/ {print $0 "\n" > "Chr05_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr06/ {print $0 "\n" > "Chr06_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr07/ {print $0 "\n" > "Chr07_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr08/ {print $0 "\n" > "Chr08_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr09/ {print $0 "\n" > "Chr09_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr10/ {print $0 "\n" > "Chr10_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr11/ {print $0 "\n" > "Chr11_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr12/ {print $0 "\n" > "Chr12_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr13/ {print $0 "\n" > "Chr13_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr14/ {print $0 "\n" > "Chr14_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr15/ {print $0 "\n" > "Chr15_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr16/ {print $0 "\n" > "Chr16_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr17/ {print $0 "\n" > "Chr17_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr18/ {print $0 "\n" > "Chr18_body.maf"}' Pkor_as_ref.maf
gawk 'BEGIN{RS=""} /Pkor.Chr19/ {print $0 "\n" > "Chr19_body.maf"}' Pkor_as_ref.maf

