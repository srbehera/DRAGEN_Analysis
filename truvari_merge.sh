source activate Truvari-4.1.0-dev

input="list.txt"
while IFS= read -r line
do
  ${bcftool}bcftools norm -N -m-any -O z -o ${line}.norm.sv.vcf.gz ${SV}${line}.sv.vcf.gz
  ${tabix} -p vcf ${line}.norm.sv.vcf.gz
  ${bcftool}bcftools norm -N -m-any -O z -o ${line}.norm.cnv.vcf.gz ${CNV}${line}.cnv.vcf.gz
  ${tabix} -p vcf ${line}.norm.cnv.vcf.gz
  ${bcftool}bcftools norm -N -m-any -O z -o ${line}.norm.repeats.vcf.gz ${STR}${line}.repeats.vcf.gz
  ${tabix} -p vcf ${line}.norm.repeats.vcf.gz
  python dragen-sv_merge.py ${line}.norm.repeats.vcf.gz ${line}.norm.sv.vcf.gz ${line}.norm.cnv.vcf.gz | ${vcftool}vcf-sort | ${bgzip} > ${line}.merged.vcf.gz
  ${tabix} -p vcf ${line}.merged.vcf.gz
  rm ${line}.norm.sv.vcf.gz
  rm ${line}.norm.sv.vcf.gz.tbi
  rm ${line}.norm.cnv.vcf.gz
  rm ${line}.norm.cnv.vcf.gz.tbi
  rm ${line}.norm.repeats.vcf.gz
  rm ${line}.norm.repeats.vcf.gz.tbi
done < "$input"

ulimit -n 5000
#${bcftool}bcftools merge -l list_merged_norm_all.txt -O z -o merged_without_norm.vcf.gz
#${tabix} -p vcf merged_without_norm.vcf.gz

${bcftool}bcftools merge -l list_merged_norm_all.txt | ${bcftool}bcftools norm --check-ref s --fasta-ref ~/Sairam/Proj6_ICA/Ref_hg38/hg38.fa -N -m-any -O z -o  merged_norm_Truvari_bcftools_3202.vcf.gz
${tabix} -p vcf merged_norm_Truvari_bcftools_3202.vcf.gz

#${bcftool}bcftools merge -l list_100_merged.txt | ${bcftool}bcftools norm -N -m-any -O z -o merged.Truvari.bcftools.100.vcf.gz
#${tabix} -p vcf merged.Truvari.bcftools.100.vcf.gz
