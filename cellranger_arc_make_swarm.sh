#!/usr/bin/env bash

# This script only works if the sample names are the last part in the directory. E.g. my/directory/to/something/samplename.csv
# Before execution, move to the same directory as desired for the output directories.
# 1st argument = file with list of ids
# 2nd argument = library directory
# 3rd argument = overall batch name. Give a name to this entire batch

SWARM_FILE=$3.swarm

# make sure there isn't already a swarm file
if test -f $SWARM_FILE; then
    rm $SWARM_FILE
fi

while IFS="" read -r p || [ -n "$p" ]
do
  printf '%s\n' "$p"
  echo "cellranger-arc count --id="$p"-ARC_test --reference=/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A\
  --libraries=$2/$p.csv \
  --localcores=24 --localmem=128" >> $SWARM_FILE
done < $1

# run the swarm file
swarm -f $SWARM_FILE -g 129 -t 24 --module cellranger-arc/2.0.2
