#!/bin/bash

#SBATCH -p short
#SBATCH -J chip
#SBATCH -e zjob-chip_%j.error
#SBATCH -o zjob-chip_%j.out
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --time=2:00:00

mm9=/primary/vari/genomicdata/genomes/mm9
db=~/szabo/brandon_estrem/db/mm9
raw=/home/brandon.estrem/szabo/brandon_estrem/03.zbtb16_mouse/analysis/raw_data/chip_Lovelace2016
out=/home/brandon.estrem/szabo/brandon_estrem/03.zbtb16_mouse/analysis/output/zbtb16_ChIP

## Analyze ChIP data from Lovelace 2016

## Format ChIP data
conda run -n sci python3 formatChIP.py $raw/zbtb16_ChIP.raw.csv $out/zbtb16_ChIP.txt

## Get ChIP sites that overlap with G4 (1424; 46% of ZBTB16 binding sites)
bedtools intersect -wa -a $out/zbtb16_ChIP.txt -b $db/nonB/kouzine_G4.mm9.bed | uniq > $out/zbtb16_ChIP_G4.bed
bedtools getfasta -fi $mm9/mm9.fa -bed $out/zbtb16_ChIP_G4.bed > $out/zbtb16_ChIP_G4.fa

## Expand G4 regions to 2kb windows (2025; 66% of ZBTB16 binding sites)
awk -v OFS='\t' -v adj=1000 '{mid = int(($3+$2)/2)} {print $1, mid - adj, mid + adj, $4, $5}' $db/nonB/kouzine_G4.mm9.bed > $db/nonB/kouzine_G4.mm9.2kb.tmp
awk -v OFS='\t' '{if ($2 < 0) $2 = 0; if ($2 < $3) print}' $db/nonB/kouzine_G4.mm9.2kb.tmp > $db/nonB/kouzine_G4.mm9.2kb.tmp2
bedtools slop -i $db/nonB/kouzine_G4.mm9.2kb.tmp2 -g "$mm9/mm9.fa.fai" -b 0 > $db/nonB/kouzine_G4.mm9.2kb.bed
rm $db/nonB/kouzine_G4.mm9.2kb.tmp $db/nonB/kouzine_G4.mm9.2kb.tmp2

bedtools intersect -wa -a $out/zbtb16_ChIP.txt -b $db/nonB/kouzine_G4.mm9.2kb.bed | uniq > $out/zbtb16_ChIP_G4.2kb.bed

## Get G4 sequences that ZBTB16 binds
bedtools intersect -wa -a $db/nonB/kouzine_G4.mm9.bed -b $out/zbtb16_ChIP.txt > $out/G4_zbtb16_ChIP.bed
bedtools getfasta -fi $mm9/mm9.fa -bed $out/G4_zbtb16_ChIP.bed > $out/G4_zbtb16_ChIP.fa

## Get G4 sequences from top 20 ZBTB16 G4 binding sites
bedtools intersect -wa -a $db/nonB/kouzine_G4.mm9.bed -b $out/zbtb16_top20.bed > $out/G4_zbtb16.top20.bed
bedtools getfasta -fi $mm9/mm9.fa -bed $out/G4_zbtb16.top20.bed > $out/G4_zbtb16.top20.fa

## Shuffle sequences to see what G4s would overlap by random chance (485; 15%)
awk -v OFS='\t' '{print $1, $2, $3}' $out/zbtb16_ChIP.txt > $out/zbtb16_ChIP.bed
bedtools shuffle -i $out/zbtb16_ChIP.bed -g $mm9/mm9.fa.fai -excl $db/mm9-blacklist.bed -seed 123 -noOverlapping > $out/zbtb16_ChIP.shuf.bed
bedtools intersect -wa -a $out/zbtb16_ChIP.shuf.bed -b $db/nonB/kouzine_G4.mm9.bed | uniq > $out/zbtb16_ChIP_G4.shuf.bed
bedtools intersect -wa -a $db/nonB/kouzine_G4.mm9.bed -b $out/zbtb16_ChIP.shuf.bed > $out/G4_zbtb16_ChIP.shuf.bed



## Get SALL4 sites that overlap with G4 (65%)
conda run -n sci python3 formatChIP.py $out/sall4 sall4_ChIP.raw.csv sall4_ChIP.txt
bedtools intersect -wa -a $out/sall4/sall4_ChIP.txt -b $db/nonB/kouzine_G4.mm9.bed | uniq > $out/sall4/sall4_ChIP_G4.txt
awk -v OFS='\t' '{if (NF==0) next} {print $0}' $out/sall4/sall4_ChIP_G4.txt > $out/sall4/sall4_ChIP_G4.bed
bedtools intersect -wa -a $db/nonB/kouzine_G4.mm9.bed -b $out/sall4/sall4_ChIP.txt > $out/sall4/G4_sall4_ChIP.bed


# MEME motif
module load bbc2/meme/meme-5.5.4
meme $out/G4_zbtb16_ChIP.fa -oc $out/meme_g4
meme $out/zbtb16_ChIP_G4.fa -oc $out/meme_chip
meme $out/zbtb16_ChIP_G4.fa -oc $out/meme_chipRC -revcomp -dna
