#!/usr/bin/env bash

# This script only works if the sample names are the last part in the directory.
# Before execution, move to the same directory as desired for the output directories.
# E.g. my/directory/to/something/samplename.csv
#
# 1st argument = file with list of ids, DO NOT include tags like -RNA or _scrn. Those are added automatically
# 2nd argument = fastq directory
# 3rd argument = overall batch name. Give a name to this entire batch. The swarm file will be called
#                your_batch_name_RNA.swarm

SWARM_FILE=$3_RNA.swarm

# make sure there isn't already a swarm file
if test -f "$SWARM_FILE"; then
    rm $SWARM_FILE
fi

while IFS="" read -r p || [ -n "$p" ]
do
  printf '%s\n' "$p"
  echo "cellranger count --chemistry=ARC-v1 --id=$p-RNA --transcriptome=/fdb/cellranger/refdata-gex-GRCh38-2020-A
  --fastqs=$2/"$p"_scrn/ \
  --localcores=12 --localmem=64" >> $SWARM_FILE
done < $1

# run the swarm file
swarm -f $SWARM_FILE -g 66 -t 12 --module cellranger/7.0.0
