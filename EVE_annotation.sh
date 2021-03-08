#!/bin/bash

#-- inputs files 
db1=$1 # only viral protein from nr without retrovirus sequences
db2=$2      # all nr proteins
hostFasta=$3

host_dmd1=$(echo $3 | cut -d'.' -f1).dmd1
host_bed1=$(echo $3 | cut -d'.' -f1).bed1
host_dmdFasta1=$(echo $3 | cut -d'.' -f1).fasta1

host_dmd2=$(echo $3 | cut -d'.' -f1).dmd2
host_bed2=$(echo $3 | cut -d'.' -f1).bed2

#-- run tools
diamond blastx --more-sensitive --max-target-seqs 1000 --range-culling --min-score 40 --outfmt 6 -b 30 -c 1 -d $db1 -q $hostFasta -o $host_dmd1 # diamond 0.9.21 

python3 bestHitsToFasta.py $host_dmd1 $host_bed1 # custom python3 script 

bedtools getfasta -fi $hostFasta -bed $host_bed1 -fo $host_dmdFasta1 # cut genomic region which have a viral best hit

diamond blastx --more-sensitive --max-target-seqs 1000  --range-culling  --min-score 40 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle -b 30 -c 1 -d $db2 -q $host_dmdFasta1 -o $host_dmd2

# taxonomizR

Rscript --vanilla absolutCoord.R $host_dmd2 $host_bed2

# script clement
