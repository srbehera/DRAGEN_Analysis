#!/bin/bash

shard=$1
in_json=../data/anno/shard-$shard/dragen.anno.json.gz
out_vcf=../data/anno/shard-$shard/dragen.anno.vcf.gz
python3 convert_anno.py $in_json $out_vcf


