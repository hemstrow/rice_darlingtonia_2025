#!/bin/bash -l
bamlist=bamlist.txt
mapQ=$1
minQ=$2
pval="${3}"
minInd=$4

# loop through each region file, make a beagle file, and then cat that onto a global file which is being produced
# do this for maf .05 and .01
for (( i=1; i<=20; i++ )); do
  trf=regionfiles/rf_${i}.txt
  echo "Region file: ${trf}."

  # maf .05
  angsd -bam $bamlist -out darlingtonia_maf05_${i} -doMajorMinor 1 -minMapQ $mapQ -minQ $minQ -SNP_pval $pval \
    -doMaf 1 -minMaf 0.05 -minInd $minInd -doGlf 2 -P 6 -GL 1 \
    -rf $trf

  zcat darlingtonia_maf05_${i}.beagle.gz | tail -q -n +2 >> darlingtonia_maf05.beagle

  # maf .01
  angsd -bam $bamlist -out darlingtonia_maf01_${i} -doMajorMinor 1 -minMapQ $mapQ -minQ $minQ -SNP_pval $pval \
    -doMaf 1 -minMaf 0.01 -minInd $minInd -doGlf 2 -P 6 -GL 1 \
    -rf $trf

  zcat darlingtonia_maf01_${i}.beagle.gz | tail -q -n +2 >> darlingtonia_maf01.beagle

done

# add headers to the combined files, clean, and zip
zcat darlingtonia_maf05_1.beagle.gz | head -n 1 > header.beagle
cat header.beagle darlingtonia_maf05.beagle > darlingtonia_maf05_combined.beagle
rm darlingtonia_maf05.beagle
gzip -c darlingtonia_maf05_combined.beagle > darlingtonia_maf05_combined.beagle.gz

zcat darlingtonia_maf01_1.beagle.gz | head -n 1 > header.beagle
cat header.beagle darlingtonia_maf01.beagle > darlingtonia_maf01_combined.beagle
rm darlingtonia_maf01.beagle
gzip -c darlingtonia_maf01_combined.beagle > darlingtonia_maf01_combined.beagle.gz

rm header.beagle

