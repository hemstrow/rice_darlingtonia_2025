#!/bin/bash -l

c1=$1
c2=$2
ref=$3

bn=$(basename "$c1" .gz)
bn=$(basename "$bn" .fq)
bn=$(basename "$bn" .fastq)
dir=$(dirname "$c1")
c1=$(basename "$c1")
c2=$(basename "$c2")

echo "Beginning filtering runs. Input files: ${dir}/${c1} and ${dir}/${c2}. Outputfile:${dir}/${bn}.sort.flt.bam"
bwa mem $ref ${dir}/${c1} ${dir}/${c2} | samtools view -Sb - | samtools sort - -n -o ${dir}/${bn}.sort.bam # align and sort by name
samtools flagstat ${dir}/${bn}.sort.bam > ${dir}/${bn}.sort.bam.flagstat
samtools fixmate -r -m ${dir}/${bn}.sort.bam ${dir}/${bn}.fixmate.bam # fixmate
rm ${dir}/${bn}.sort.bam
samtools sort -o ${dir}/${bn}.psort.bam ${dir}/${bn}.fixmate.bam # sort by position
rm ${dir}/${bn}.fixmate.bam
samtools markdup -r ${dir}/${bn}.psort.bam ${dir}/${bn}.markdup.bam # remove dups
rm ${dir}/${bn}.psort.bam
samtools view -q 5 -b ${dir}/${bn}.markdup.bam > ${dir}/${bn}.q1.bam # remove poorly mapped
rm ${dir}/${bn}.markdup.bam
samtools sort -n -o ${dir}/${bn}.namesort.bam ${dir}/${bn}.q1.bam # sort by name again
rm ${dir}/${bn}.q1.bam
samtools fixmate -m ${dir}/${bn}.namesort.bam ${dir}/${bn}.fixmate.bam # fix mates again
rm ${dir}/${bn}.namesort.bam
samtools sort -o ${dir}/${bn}.sort.flt.bam ${dir}/${bn}.fixmate.bam # sort
rm ${dir}/${bn}.fixmate.bam
samtools index ${dir}/${bn}.sort.flt.bam # index

# stats
samtools stats ${dir}/${bn}.sort.flt.bam > ${dir}/${bn}.sort.flt.bam.stats
samtools flagstat ${dir}/${bn}.sort.flt.bam > ${dir}/${bn}.sort.flt.bam.flagstat
