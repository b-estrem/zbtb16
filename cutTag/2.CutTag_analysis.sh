#!/bin/bash

#SBATCH -p short
#SBATCH -J CutTagPeaks
#SBATCH -e zjob-CutTagPeaks_%j.error
#SBATCH -o zjob-CutTagPeaks_%j.out
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --time=12:00:00
t=$SLURM_CPUS_PER_TASK

## Call and analyze peaks from mapped CUT&Tag data
module load bbc2/macs2/macs2-2.2.7.1

base=/home/brandon.estrem/szabo/brandon_estrem/03.zbtb16_mouse/analysis/output/cut_tag
db=/home/brandon.estrem/szabo/brandon_estrem/db/mm10/nonB

## Tried narrow peak, broad seems better

##H3K4me3
# for file in $base/bed/*H3K4me3_070126.bed
# do
#     name="$(basename -- $file .bed)"
#     ctrl_file="$(basename -- $file H3K4me3_070126.bed)"
#     macs2 callpeak -t $file -g mm --keep-dup auto --nomodel --extsize 147 --outdir $base/macs --broad -n $name ## Broad peak
#     macs2 callpeak -t $file -g mm --keep-dup auto --nomodel --extsize 147 --outdir $base/macs -c $base/bed/$ctrl_file\IgG_070126.bed --broad -n $name.ctrl ## Broad peak with IgG ctrl
# done

# ## G4
# for file in $base/bed/*G4_070126.bed
# do
#     name="$(basename -- $file .bed)"
#     ctrl_file="$(basename -- $file G4_070126.bed)"
#     macs2 callpeak -t $file -g mm --keep-dup auto --nomodel --extsize 147 --outdir $base/macs --broad -n $name ## Broad peak
#     macs2 callpeak -t $file -g mm --keep-dup auto --nomodel --extsize 147 --outdir $base/macs -c $base/bed/$ctrl_file\IgG_070126.bed --broad -n $name.ctrl ## Broad peak with IgG ctrl

#     bedtools intersect -a $base/macs/$name\_peaks.broadPeak -b $db/G4.bed -u > $base/macs/intersect/$name.G4.bed
#     bedtools intersect -a $base/macs/$name.ctrl_peaks.broadPeak -b $db/G4.bed -u > $base/macs/intersect/$name.G4_ctrl.bed
# done

# ## IgG
# for file in $base/bed/*IgG_070126.bed
# do
#     name="$(basename -- $file .bed)"
#     macs2 callpeak -t $file -g mm --keep-dup auto --nomodel --extsize 147 --outdir $base/macs --broad -n $name ## Broad peak

#     bedtools intersect -a $base/macs/$name\_peaks.broadPeak -b $db/G4.bed -u > $base/macs/intersect/$name.G4.bed
# done


## Intended to compare reads/G4 but decided too crude, MACS is better. Did not do analysis
# for file in $base/bigwig/*G4*.bigwig $base/bigwig/*IgG*.bigwig
# do
#     name="$(basename -- $file .bigwig)"
#     bigWigToBedGraph $file $name.bg

#     bedtools coverage -a $db/G4.slop10.bed -b $base/bigwig/$name.bg > $base/bigwig/G4cov/$name.cov
# done