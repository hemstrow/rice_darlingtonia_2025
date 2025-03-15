#!/bin/bash -l
#SBATCH --mem=20G
#SBATCH -t 6-24:00:00
#SBATCH -J rad_split

read_dir="/home/hemstrow/cody/additional_reads"
r1=SOMM_514_NAGTGG_R1.fq
r2=SOMM_514_NAGTGG_R2.fq
barcodes="barcodes.txt"
output="/home/hemstrow/cody/additional_reads/SOMM_514_radsplit_"
leading_bases=2
leading_errors=1
barcode_errors=0
chunk_size=1000


Rscript run_rad_split.R $read_dir $r1 $r2 $barcodes $output $leading_bases $leading_errors $barcode_errors $chunk_size