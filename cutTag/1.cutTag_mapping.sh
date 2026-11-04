#!/bin/bash

#SBATCH -p short
#SBATCH -J CutTag
#SBATCH -e zjob-mapCutTag_%j.error
#SBATCH -o zjob-mapCutTag_%j.out
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --time=2:00:00
t=$SLURM_CPUS_PER_TASK

tmp=~/hpctmp/brandon.estrem/cutTag
raw=/home/brandon.estrem/szabo/brandon_estrem/03.zbtb16_mouse/analysis/raw_data/cuttag_20260810
dm6=/home/brandon.estrem/szabo/brandon_estrem/db/dm6
mm9=/primary/vari/genomicdata/genomes/mm9/bowtie2
db=~/szabo/brandon_estrem/db
nonB=/home/brandon.estrem/szabo/brandon_estrem/db/mm9/nonB

name=BE_mm_CUTTag_WT_A_G4_070126

#################################################
##### Trim and QC of data
module load bbc2/fastqc/fastqc-0.12.1 bbc2/Trimmomatic/trimmomatic-0.39

## Check each step output and see how many reads lost, determine rm or keep dups
for name in BE_mm_CUTTag_WT_A_G4_070126 BE_mm_CUTTag_WT_A_H3K4me3_070126
do
    # fastqc $raw/$name.R1.fq.gz

    # java -jar $TRIMM PE -phred33 -threads $t $raw/$name.R1.fq.gz $raw/$name.R2.fq.gz -baseout $tmp/$name.trimmed.fq.gz ILLUMINACLIP:$db/NexteraPE-PE.fa:3:30:10:3:True MINLEN:20
    # fastqc $tmp/$name.trimmed_1P.fq.gz

    ## Map to Drosophila Genome for normalization
    ## test # of reads that map to both genomes - need to estimate library size?
    bowtie2 -p $t -x $dm6/dm6 -1 $tmp/$name.trimmed_1P.fq.gz -2 $tmp/$name.trimmed_2P.fq.gz -S $out/$name.dm6.sam --end-to-end --very-sensitive --nounal --no-mixed --no-discordant -I 10 -X 700

    bowtie2 -p $t -x $mm9/mm9 -1 $tmp/$name.trimmed_1P.fq.gz -2 $tmp/$name.trimmed_2P.fq.gz -S $out/$name.sam --end-to-end --very-sensitive --nounal --no-mixed --no-discordant -I 10 -X 700
    
    samtools view -@ $t -bf 3 -q 30 $out/$name.sam | samtools sort -@ $t - -o $out/$name.bam

    module load bbc2/picard/picard-3.0.0
    java -Xms8g -Xmx16g -jar $PICARD MarkDuplicates \
        -I $out/$name.bam \
        -O $out/$name.noDup.bam \
        -M $out/$name.duplicates.txt \
        --REMOVE_DUPLICATES true \
        --QUIET true \
        --VERBOSITY WARNING
    samtools index -@ $t $out/$name.bam

    ## Normalize by reads

    ## Normalize by spike in

    ## Make bg/bw file?

    ## MACS/SEACR, overlap with G4

    ## FRIP? other analysis?
done
module unload bbc2/fastqc/fastqc-0.12.1 bbc2/Trimmomatic/trimmomatic-0.39

# ## Map reads to drosophila genome for normalization

# ## Map reads to mm10
# bowtie2 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700