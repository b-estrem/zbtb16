#!/bin/bash

#SBATCH -p long
#SBATCH -J CutTag
#SBATCH -e zjob-mapCutTag_%j.error
#SBATCH -o zjob-mapCutTag_%j.out
#SBATCH --nodes=1
#SBATCH -c 64
#SBATCH --time=10:00:00
t=$SLURM_CPUS_PER_TASK

raw=/home/brandon.estrem/szabo/brandon_estrem/03.zbtb16_mouse/analysis/raw_data/cuttag_20260810
tmp=~/hpctmp/brandon.estrem/cutTag
db=~/szabo/brandon_estrem/db
mm10=/primary/vari/genomicdata/genomes/mm10
hybrid=$db/hybrid_mm10_dm6
out=/home/brandon.estrem/szabo/brandon_estrem/03.zbtb16_mouse/analysis/output/cutTag

#################################################
##### Things that only need to be run once
## Make concatenated genome - better than mapping to each genome individually
cp /primary/vari/genomicdata/genomes/dm6/dm6.fa $db
## prefix spike-in chromosomes so they can be separated later
sed 's/^>/>dm_/' /primary/vari/genomicdata/genomes/dm6/dm6.fa > $hybrid/dm6_prefixed.fa
cat $mm10/mm10.fa $hybrid/dm6_prefixed.fa > $hybrid/hybrid_mm10_dm6.fa
bowtie2-build --threads 16 $hybrid/hybrid_mm10_dm6.fa $hybrid/hybrid_mm10_dm6

## Make libStats file with header
echo -e "name\tfastq\tmm10_reads\tdm6_reads" > libStats.tsv

#################################################
##### Trim and QC of data
module load bbc2/fastqc/fastqc-0.12.1 bbc2/Trimmomatic/trimmomatic-0.39

for name in $raw/$name.R1.fq.gz
do
    fastqc $raw/$name.R1.fq.gz

    java -jar $TRIMM PE -phred33 -threads $t $raw/$name.R1.fq.gz $raw/$name.R2.fq.gz -baseout $tmp/$name.trimmed.fq.gz ILLUMINACLIP:$db/NexteraPE-PE.fa:3:30:10:3:True MINLEN:20
    fastqc $tmp/$name.trimmed_1P.fq.gz

    ## Map to hybrid genome for normalization
    echo "$name" >&2
    bowtie2 -p $t -x $hybrid/hybrid_mm10_dm6 -1 $tmp/$name.trimmed_1P.fq.gz -2 $tmp/$name.trimmed_2P.fq.gz -S $out/bam/$name.hybrid.sam --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 700
    
    samtools view -@ $t -bf 3 -q 30 $out/bam/$name.hybrid.sam | samtools sort -@ $t - -o $out/bam/$name.hybrid.bam
    samtools index -@ $t $out/bam/$name.hybrid.bam
    rm $out/bam/$name.hybrid.sam
    
    ## Remove dm6 reads, leaving behind mm10 mapped reads
    chroms=$(samtools idxstats $out/bam/$name.hybrid.bam | cut -f1 | grep -v -e '^dm_' -e '^\*$')
    samtools view -@ $t -b $out/bam/$name.hybrid.bam $chroms > $out/bam/$name.mm10.bam
    samtools index $out/bam/$name.mm10.bam

    ## Calculate normalization by spike in
    dm_reads=$(samtools idxstats $out/bam/$name.hybrid.bam | awk '$1 ~ /^dm_/ {n+=$3} END{print n}')
    mm10_reads=$(samtools idxstats $out/bam/$name.hybrid.bam | awk '$1 !~ /^dm_/ && $1!="*" {n+=$3} END{print n}')
    norm_dm6="$(echo "scale=3; 10000/$reads" | bc)"
    
    ## Calculate normalization to 10M reads 
    reads="$(samtools view -c $out/bam/$name.mm10.bam)"
    norm="$(echo "scale=3; 10000000/$reads" | bc)"

    ## Add read counts to file
    lines=$(zcat "$raw/$name.R1.fq.gz" | wc -l | tr -d ' ') # Lines in a Fq file
    fq=$(( lines / 4 * 2 )) # Find num of reads, multiply by two for num paired reads
    echo -e "${name}\t${fq}\t${reads}\t${dm_reads}" >> libStats.tsv
    
    ## Make bg/bw file
    bedtools genomecov -ibam ${out}/bam/${name}.mm10.bam -bg -scale ${norm_dm6} | sort -k1,1 -k2,2n > $out/coverage/$name.norm_dm6.bg
    bedGraphToBigWig $out/coverage/$name.norm_dm6.bg $mm10/mm10.fa.fai $out/coverage/$name.norm_dm6.bw
    
    bedtools genomecov -ibam ${out}/bam/${name}.mm10.bam -bg -scale ${norm} | sort -k1,1 -k2,2n > $out/coverage/$name.norm10M.bg
    bedGraphToBigWig $out/coverage/$name.norm10M.bg $mm10/mm10.fa.fai $out/coverage/$name.norm10M.bw
done

module unload bbc2/fastqc/fastqc-0.12.1 bbc2/Trimmomatic/trimmomatic-0.39
