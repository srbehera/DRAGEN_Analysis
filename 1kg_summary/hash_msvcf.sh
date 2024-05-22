#!/bin/bash

shard=$1
in_msvcf=../data/msvcfs/shard-$shard/dragen.vcf.gz
out_vcf=../data/msvcfs.hash/shard-$shard/dragen.vcf.gz
mkdir -p $(dirname $out_vcf)
python3 hash_msvcf.py $in_msvcf $out_vcf
