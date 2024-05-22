#!/bin/bash

shard=$1
in_vcf=../data/var_info/shard-$shard/dragen.gtc.vcf.gz
out_json=../data/var_stats/shard-$shard/dragen.count_var_cohort.json
python3 count_var_cohort.py $in_vcf $out_json 1kg3202 false
    
