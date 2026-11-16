#!/bin/bash

#SBATCH -p long
#SBATCH -J CutTag
#SBATCH -e zjob-mapCutTag_%j.error
#SBATCH -o zjob-mapCutTag_%j.out
#SBATCH --nodes=1
#SBATCH -c 64
#SBATCH --time=10:00:00
set -euo pipefail
t=$SLURM_CPUS_PER_TASK

raw=/home/brandon.estrem/szabo/brandon_estrem/03.zbtb16_mouse/analysis/raw_data/cuttag_20260810
tmp=~/hpctmp/brandon.estrem/cutTag
db=~/szabo/brandon_estrem/db
mm10=/primary/vari/genomicdata/genomes/mm10
hybrid=$db/hybrid_mm10_dm6
out=/home/brandon.estrem/szabo/brandon_estrem/03.zbtb16_mouse/analysis/output/cutTag

# #################################################
# ##### Things that only need to be run once
# ## Make concatenated genome - better than mapping to each genome individually
# cp /primary/vari/genomicdata/genomes/dm6/dm6.fa $db
# ## prefix spike-in chromosomes so they can be separated later
# sed 's/^>/>dm_/' /primary/vari/genomicdata/genomes/dm6/dm6.fa > $hybrid/dm6_prefixed.fa
# cat $mm10/mm10.fa $hybrid/dm6_prefixed.fa > $hybrid/hybrid_mm10_dm6.fa
# bowtie2-build --threads 16 $hybrid/hybrid_mm10_dm6.fa $hybrid/hybrid_mm10_dm6

# ## Make libStats file with header
# echo -e "name\tfastq\tmm10_reads\tdm6_reads" > libStats.tsv

#################################################
##### Trim and QC of data
module load bbc2/fastqc/fastqc-0.12.1 bbc2/Trimmomatic/trimmomatic-0.39 bbc2/cutadapt/cutadapt-4.3

for file in $raw/*.R1.fq.gz
do
    name="$(basename -- $file .R1.fq.gz)"
    echo "Processing $name" >&2

    # java -jar $TRIMM PE -phred33 -threads $t $raw/$name.R1.fq.gz $raw/$name.R2.fq.gz -baseout $tmp/$name.trimmed.fq.gz ILLUMINACLIP:$db/NexteraPE-PE.fa:2:30:10:2:True MINLEN:20
    echo "Trimming adapters with Cutadapt..." >&2
    cutadapt --report=minimal -j 0 -m 20 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o $tmp/$name.trimmed.R1.fq.gz -p $tmp/$name.trimmed.R2.fq.gz $raw/$name.R1.fq.gz $raw/$name.R2.fq.gz
    fastqc -q $tmp/$name.trimmed.R1.fq.gz

    ## Map to hybrid genome for normalization
    echo "Mapping to hybrid genome with bowtie2..." >&2
    echo "$name" >&2
    bowtie2 -p $t -x $hybrid/hybrid_mm10_dm6 -1 $tmp/$name.trimmed.R1.fq.gz -2 $tmp/$name.trimmed.R2.fq.gz -S $out/bam/$name.hybrid.sam --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 1000
    
    samtools view -@ $t -bf 3 -q 30 $out/bam/$name.hybrid.sam | samtools sort -@ $t - -o $out/bam/$name.hybrid.bam
    samtools index -@ $t $out/bam/$name.hybrid.bam
    rm $out/bam/$name.hybrid.sam
    
    ## Remove dm6 reads, leaving behind mm10 mapped reads
    echo "Making an mm10 only bam file..." >&2
    chroms=$(samtools idxstats $out/bam/$name.hybrid.bam | cut -f1 | grep -v -e '^dm_' -e '^\*$')
    samtools view -@ $t -b $out/bam/$name.hybrid.bam $chroms > $out/bam/$name.mm10.bam
    samtools index $out/bam/$name.mm10.bam

    ## Calculate normalization by spike in and 10M reads
    echo "Calculating normalization factors..." >&2
    dm_reads=$(samtools idxstats $out/bam/$name.hybrid.bam | awk '$1 ~ /^dm_/ {n+=$3} END{print n}')
    mm10_reads=$(samtools idxstats $out/bam/$name.hybrid.bam | awk '$1 !~ /^dm_/ && $1!="*" {n+=$3} END{print n}')
    norm_dm6="$(echo "scale=6; 10000/$dm_reads" | bc)"
    norm_mm10="$(echo "scale=6; 10000000/$mm10_reads" | bc)"

    # ## Add read counts to file
    # lines=$(zcat "$raw/$name.R1.fq.gz" | wc -l | tr -d ' ') # Lines in a Fq file
    # fq=$(( lines / 4 * 2 )) # Find num of reads, multiply by two for num paired reads
    # echo -e "${name}\t${fq}\t${reads}\t${dm_reads}" >> libStats.tsv
    
    ## Convert to bedpe
    echo "Convert bam to bed file of fragments..." >&2
    samtools sort $out/bam/$name.mm10.bam -n -@ $t -o $out/bam/$name.nsort.bam 
    bedtools bamtobed -bedpe -i $out/bam/$name.nsort.bam |
    awk -v OFS='\t' '$1==$4 && $1!="." {
        s = ($2<$5) ? $2 : $5;
        e = ($3>$6) ? $3 : $6;
        print $1, s, e }' |
    sort -k1,1 -k2,2n -k3,3n > $out/coverage/$name.frag.bed
    rm $out/bam/$name.nsort.bam 

    ## Make bg/bw file
    echo "Convert to normalized bw..." >&2
    if [[ $dm_reads > 100 ]]; then 
        bedtools genomecov -i $out/coverage/$name.frag.bed -g $mm10/mm10.fa.fai -bg -scale $norm_dm6 | sort -k1,1 -k2,2n > $out/coverage/$name.norm_dm6.bg
        bedGraphToBigWig $out/coverage/$name.norm_dm6.bg $mm10/mm10.fa.fai $out/coverage/$name.norm_dm6.bw
    fi 

    bedtools genomecov -i $out/coverage/$name.frag.bed -g $mm10/mm10.fa.fai -bg -scale $norm_mm10 | sort -k1,1 -k2,2n > $out/coverage/$name.norm10M.bg
    bedGraphToBigWig $out/coverage/$name.norm10M.bg $mm10/mm10.fa.fai $out/coverage/$name.norm10M.bw
done

module unload bbc2/fastqc/fastqc-0.12.1 bbc2/Trimmomatic/trimmomatic-0.39
