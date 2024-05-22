#!/bin/bash

shard=$1
gtc_fname=../data/var_info/shard-$shard/dragen.gtc.vcf.gz
anno_fname=../data/anno/shard-$shard/dragen.anno.vcf.gz
out_fname=../data/var_info/shard-$shard/dragen.gtc-anno.vcf.gz
python3 join_gtc_anno.py $gtc_fname $anno_fname $out_fname
