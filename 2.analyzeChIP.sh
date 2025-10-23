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
out=/home/brandon.estrem/szabo/brandon_estrem/02.zbtb16/analysis/output/zbtb16_ChIP

## Analyze ChIP data from Lovelace 2016

## Format ChIP data
conda run -n sci python3 formatChIP.py $out zbtb16_ChIP.raw.csv zbtb16_ChIP.txt

# ## Get ChIP sites that overlap with G4 (46% of ZBTB16 binding sites)
bedtools intersect -wa -a $out/zbtb16_ChIP.txt -b $db/nonB/kouzine_G4.mm9.bed | uniq > $out/zbtb16_ChIP_G4.txt
awk -v OFS='\t' '{if (NF==0) next} {print $0}' $out/zbtb16_ChIP_G4.txt > $out/zbtb16_ChIP_G4.bed
bedtools getfasta -fi $mm9/mm9.fa -bed $out/zbtb16_ChIP_G4.bed > $out/zbtb16_ChIP_G4.fa

## Get G4 sequences that ZBTB16 binds
bedtools intersect -wa -a $db/nonB/kouzine_G4.mm9.bed -b $out/zbtb16_ChIP.txt > $out/G4_zbtb16_ChIP.bed
bedtools getfasta -fi $mm9/mm9.fa -bed $out/G4_zbtb16_ChIP.bed > $out/G4_zbtb16_ChIP.fa

## Get G4 sequences from top 20 ZBTB16 G4 binding sites
bedtools intersect -wa -a $db/nonB/kouzine_G4.mm9.bed -b $out/zbtb16_top20.bed > $out/G4_zbtb16.top20.bed
bedtools getfasta -fi $mm9/mm9.fa -bed $out/G4_zbtb16.top20.bed > $out/G4_zbtb16.top20.fa

## Shuffle sequences to see what G4s would overlap by random chance
awk -v OFS='\t' '{print $1, $2, $3}' $out/zbtb16_ChIP.txt > $out/zbtb16_ChIP.bed
bedtools shuffle -i $out/zbtb16_ChIP.bed -g $mm9/mm9.fa.fai -excl $db/mm9-blacklist.bed -seed 123 -noOverlapping > $out/zbtb16_ChIP.shuf.bed
bedtools intersect -wa -a $out/zbtb16_ChIP.shuf.bed -b $db/nonB/kouzine_G4.mm9.bed | uniq > $out/zbtb16_ChIP_G4.shuf.txt
bedtools intersect -wa -a $db/nonB/kouzine_G4.mm9.bed -b $out/zbtb16_ChIP.shuf.bed > $out/G4_zbtb16_ChIP.shuf.bed




## Get SALL4 sites that overlap with G4 
conda run -n sci python3 formatChIP.py $out/sall4 sall4_ChIP.raw.csv sall4_ChIP.txt
bedtools intersect -wa -a $out/sall4/sall4_ChIP.txt -b $db/nonB/kouzine_G4.mm9.bed | uniq > $out/sall4/sall4_ChIP_G4.txt
awk -v OFS='\t' '{if (NF==0) next} {print $0}' $out/sall4/sall4_ChIP_G4.txt > $out/sall4/sall4_ChIP_G4.bed
bedtools intersect -wa -a $db/nonB/kouzine_G4.mm9.bed -b $out/sall4/sall4_ChIP.txt > $out/sall4/G4_sall4_ChIP.bed


# MEME motif
module load bbc2/meme/meme-5.5.4
meme $out/G4_zbtb16_ChIP.fa -oc $out/meme_g4
meme $out/zbtb16_ChIP_G4.fa -oc $out/meme_chip
meme $out/zbtb16_ChIP_G4.fa -oc $out/meme_chipRC -revcomp -dna
