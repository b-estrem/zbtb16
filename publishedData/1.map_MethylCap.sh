#!/bin/bash

#SBATCH -p short
#SBATCH -J methylCap
#SBATCH -e zjob-meCap_%j.error
#SBATCH -o zjob-meCap_%j.out
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --time=48:00:00
t=$SLURM_CPUS_PER_TASK

mm10=/primary/vari/genomicdata/genomes/mm10
db=~/szabo/brandon_estrem/db/mm10
raw=~/szabo/brandon_estrem/02.zbtb16/analysis/raw_data/mCap_maezawa2025
out=~/szabo/brandon_estrem/02.zbtb16/analysis/output

# module load bbc2/picard/picard-3.0.0

# for file in $raw/*.fq.gz
# do
#     name="$(basename -- $file .fq.gz)"
#     >&2 echo $name 

#     bowtie2 -p $t --local -x $mm10/mm10 -U $file -S $out/bam/$name.sam
#     samtools view -@ $t -b -q 30 $out/bam/$name.sam | samtools sort -@ $t - -o $out/bam/$name.noMask.bam
#     bedtools intersect -v -a $out/bam/$name.noMask.bam -b $db/mm10-blacklist.v2.bed > $out/bam/$name.mask.bam
#     java -Xms8g -Xmx16g -jar $PICARD MarkDuplicates \
#         -I $out/bam/$name.mask.bam \
#         -O $out/bam/$name.bam \
#         -M $out/bam/$name.duplicates.txt \
#         --REMOVE_DUPLICATES true \
#         --QUIET true \
#         --VERBOSITY WARNING
#     rm $out/bam/$name.sam $out/bam/$name.noMask.bam $out/bam/$name.mask.bam

#     ## Calculate normalization factor
#     reads="$(samtools view -c $out/bam/$name.bam)"
#     norm="$(echo "scale=3; 100000000/$reads" | bc)"

#     ## Create files for igv
#     bedtools genomecov -scale $norm -ibam $out/bam/$name.bam -bg | sort -k1,1 -k2,2n > $out/coverage/$name.norm.bg
#     bedGraphToBigWig $out/coverage/$name.norm.bg $mm10/mm10.fa.fai $out/coverage/$name.norm.bw
# done

# ## picard interferes with running R script
# module unload bbc2/picard/picard-3.0.0

## conda create -n mCap_env bioconda::bioconductor-medips bioconda::bioconductor-bsgenome r-tidyverse
source /home/brandon.estrem/miniconda3/etc/profile.d/conda.sh
conda activate mCap_env
Rscript medips.R