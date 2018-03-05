#!/usr/bin/env bash

# Recursively get all files for each cancer.
find ./Cancers/ -name \*.gz -print > AllCancerFiles.txt
# Obtain only the individual patient files.
cat ./AllCancerFiles.txt | grep -v ".snvs.indels.maf.gz" > AllCancerFilesOut.txt
rm ./AllCancerFiles.txt

# For each file get the counts of mutation types.
while read line; do
  echo $line
  gzip -cd $line | awk '{ print $6 }' | awk '{ print $6 }' | sort | uniq -c | sort -nr
done <AllCancerFiles.txt
