#!/bin/bash

echo "Retreiving Canonical PSet URLs...\n"
# urlArr=(`Rscript -e 'library(PharmacoGx); fl <- availablePSets(); cat(paste0(fl$Download, " "))'`)
urlArr=(`cut -d ',' -f9 availablePSets.csv`)
echo "DONE\n"

echo "Beginning PSet downloads...\n"
for url in ${urlArr[@]:2}; do
    echo "...Downloading from $url...\n"
    wget $url
    echo "...Done\n\n"
done
echo "DONE\n"

echo "Trimming PSet file names to remove ?download=1..."
fileArr=(`ls | grep download=1`)
for file in ${fileArr[@]}; do
    fileName=`echo $file | sed s/rds.*/rds/g`
    mv $file $fileName
    echo "$file renamed $fileName\n"
done
echo ""

echo "downloadPSets.sh finished.\n"