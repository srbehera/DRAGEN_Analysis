#!/bin/bash

shard=$1
in_vcf=../data/sites/shard-$shard/dragen.sites.vcf.gz
site_vcf=$json_folder/shard-$shard/dragen.sites.4col.vcf.gz
out_json=$json_folder/shard-$shard/dragen.anno
mkdir -p $(dirname $out_json)

echo "writing $site_vcf"
zcat $in_vcf | awk '{if (/^#/) print; else printf("%s\t%s\t.\t%s\t%s\t.\t.\t.\n", $1,$2,$4,$5)}' | bgzip -f > $site_vcf
tabix -f $site_vcf

echo "writing $out_json"
/usr/bin/time -v $dotnet \
    $nirvana_dll \
    -c $nirvana_cache \
    --sd $nirvana_suppl \
    -r $nirvana_ref \
    -i $site_vcf \
    -o $out_json \
    2>&1 | tee $out_json.log

