
ns_cut=$1
in_fname=../data/var_info/shard-$shard/dragen.gtc-anno.vcf.gz
out_vcf_fname=../data/var_stats/shard-$shard/dragen.freq-novel-func.nscut-$ns_cut.vcf.gz
out_json_fname=../data/var_stats/shard-$shard/dragen.freq-novel-func.nscut-$ns_cut.json
pass_only=false

python3 get_allele_freq_novel_func.py $in_fname $out_vcf_fname $out_json_fname $pass_only $ns_cut
