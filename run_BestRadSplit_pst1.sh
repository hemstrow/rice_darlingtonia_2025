#!/bin/bash
#SBATCH --mem=20G
#SBATCH -t 6-12:00:00
#SBATCH -J p1_rad_split

f1=SOMM474_R1_CAACAT.fastq
f2=SOMM474_R3_CAACAT.fastq
out=plate_3

perl ~/scripts/BarcodeSplitListBestRadPairedEnd.pl $f1 $f2 GGACAAGCTATGCAG,GGAAACATCGTGCAG,GGACATTGGCTGCAG,GGACCACTGTTGCAG,GGAACGTGATTGCAG,GGCGCTGATCTGCAG,GGCAGATCTGTGCAG,GGATGCCTAATGCAG,GGAACGAACGTGCAG,GGAGTACAAGTGCAG,GGCATCAAGTTGCAG,GGAGTGGTCATGCAG,GGAACAACCATGCAG,GGAACCGAGATGCAG,GGAACGCTTATGCAG,GGAAGACGGATGCAG,GGAAGGTACATGCAG,GGACACAGAATGCAG,GGACAGCAGATGCAG,GGACCTCCAATGCAG,GGACGCTCGATGCAG,GGACGTATCATGCAG,GGACTATGCATGCAG,GGAGAGTCAATGCAG,GGAGATCGCATGCAG,GGAGCAGGAATGCAG,GGAGTCACTATGCAG,GGATCCTGTATGCAG,GGATTGAGGATGCAG,GGCAACCACATGCAG,GGCAAGACTATGCAG,GGCAATGGAATGCAG,GGCACTTCGATGCAG,GGCAGCGTTATGCAG,GGCATACCAATGCAG,GGCCAGTTCATGCAG,GGCCGAAGTATGCAG,GGCCGTGAGATGCAG,GGCCTCCTGATGCAG,GGCGAACTTATGCAG,GGCGACTGGATGCAG,GGCGCATACATGCAG,GGCTCAATGATGCAG,GGCTGAGCCATGCAG,GGCTGGCATATGCAG,GGGAATCTGATGCAG,GGGACTAGTATGCAG,GGGAGCTGAATGCAG,GGGATAGACATGCAG,GGGCCACATATGCAG,GGGCGAGTAATGCAG,GGGCTAACGATGCAG,GGGCTCGGTATGCAG,GGGGAGAACATGCAG,GGGGTGCGAATGCAG,GGGTACGCAATGCAG,GGGTCGTAGATGCAG,GGGTCTGTCATGCAG,GGGTGTTCTATGCAG,GGTAGGATGATGCAG,GGTATCAGCATGCAG,GGTCCGTCTATGCAG,GGTCTTCACATGCAG,GGTGAAGAGATGCAG,GGTGGAACAATGCAG,GGTGGCTTCATGCAG,GGTGGTGGTATGCAG,GGTTCACGCATGCAG,GGACACGAGATGCAG,GGAAGAGATCTGCAG,GGAAGGACACTGCAG,GGAATCCGTCTGCAG,GGAATGTTGCTGCAG,GGACACTGACTGCAG,GGACAGATTCTGCAG,GGAGATGTACTGCAG,GGAGCACCTCTGCAG,GGAGCCATGCTGCAG,GGAGGCTAACTGCAG,GGATAGCGACTGCAG,GGACGACAAGTGCAG,GGATTGGCTCTGCAG,GGCAAGGAGCTGCAG,GGCACCTTACTGCAG,GGCCATCCTCTGCAG,GGCCGACAACTGCAG,GGAGTCAAGCTGCAG,GGCCTCTATCTGCAG,GGCGACACACTGCAG,GGCGGATTGCTGCAG,GGCTAAGGTCTGCAG,GGGAACAGGCTGCAG,GGGACAGTGCTGCAG,GGGAGTTAGCTGCAG,GGGATGAATCTGCAG,GGGCCAAGACTGCAG $out

