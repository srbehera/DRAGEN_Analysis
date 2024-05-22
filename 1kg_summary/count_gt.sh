#!/bin/bash

shard=$1
infile=../data/msvcfs.hash/shard-$shard/dragen.vcf.gz
outfile=../data/var_info/shard-$shard/dragen.gtc.vcf.gz
mkdir -p $(dirname $outfile)
python3 count_gt.py $infile $cohort_def $outfile 
