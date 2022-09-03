#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l scratch=50G
#$ -l h_rt=24:0:0
#$ -pe smp 4
#$ -t 1-12 

## RNA-seq alignment with Hisat and FeatureCounts


DIR=RNASeq_2022-02/
READDIR=$DIR/reads/ #location of reads in fastq format

FILE="ArraySubmissionFile.txt" #a tab delimited file in which filenames are specified
CountMethod="FeatureCounts" #FeatureCounts or FADU

SAMPLE=$(sed "${SGE_TASK_ID}q;d" $FILE | cut -f1)
#GTF="GCF_000024265.1_ASM2426v1_genomic.gtf"
INDEX="Elenta2243_hisat2"
GFF="GCF_000024265.1_ASM2426v1_genomic.gff"

# hisat2-build GCF_000024265.1_ASM2426v1_genomic.fna Elenta2243_hisat2

echo $(date) $SAMPLE $INDEX

RUNTMP=$( mktemp -d -t -p /scratch RNASeq_XXXXXX ) #<-----------Here
echo $RUNTMP

cp $READDIR/${SAMPLE}* $RUNTMP

mkdir -p quality
fastp \
  --in1 ${RUNTMP}/${SAMPLE}*_R1_001.fastq.gz \
  --in2 ${RUNTMP}/${SAMPLE}*_R2_001.fastq.gz \
  --out1 ${RUNTMP}/${SAMPLE}*_R1_001.filt.fastq.gz \
  --out2 ${RUNTMP}/${SAMPLE}*_R2_001.filt.fastq.gz \
  --trim_poly_g \
  --cut_front \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 15 \
  --json quality/${SAMPLE}.json \
  --html quality/${SAMPLE}.html \
  --thread $NSLOTS
  
mkdir -p alignments4
conda activate hisat2_2.2.1

  #Skip multimaps or get them all (-a)
hisat2 -x $INDEX \
     -1 "$READDIR/${SAMPLE}_R1_001.fastq.gz" \
     -2 "$READDIR/${SAMPLE}_R2_001.fastq.gz" \
     -X 1000 -a --no-spliced-alignment \
     --very-sensitive\
     --threads $NSLOTS \
     --un-conc-gz alignments4/${SAMPLE}"_unaligned.fastq.gz" -S alignments4/${SAMPLE}.sam

samtools view -S -b alignments4/${SAMPLE}.sam | samtools sort > alignments4/${SAMPLE}.bam
rm alignments4/${SAMPLE}.sam
samtools index alignments4/${SAMPLE}.bam alignments4/${SAMPLE}.bai



conda deactivate

