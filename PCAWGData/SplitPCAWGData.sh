#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir "$DIR/Cancers"

while read line
do
  echo "${line//-/}"
  mkdir "$DIR/Cancers/${line//-/}"
  gzip -cd $DIR/final_consensus_12aug_passonly_snv_indel.maf.gz | grep "$line" | gzip > "$DIR/Cancers/${line//-/}/$line.snvs.indels.maf.gz"
done < CancerTypes.txt