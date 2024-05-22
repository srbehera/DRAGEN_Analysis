#!/bin/bash

shard=$1
in_vcf=../data/msvcfs.hash/shard-$shard/dragen.vcf.gz
out_json=../data/var_stats/shard-$shard/dragen.count_var_sample.json
python3 count_var_sample.py $in_vcf $out_json false
    
