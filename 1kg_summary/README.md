# DRAGEN 1KG reprocess summary stats scripts

This is a collection of scripts used for DRAGEN 1KG summary stats. 

Please contact Zhuoyi Huang <zhuang@illumina.com> for questions or access of end-to-end workflow.

Below is the order of execution.

## GT stats
- hash_msvcf.sh
- count_gt.sh
- count_var_cohort.sh
- count_var_sample.sh

## annotation and join with gt stats 
- annotate.sh
- convert_anno.sh
- join_gtc_anno.sh
- get_allele_freq_novel_func.sh

## cohort, sample level stats
- summary1/summarize_cohort_var_count.py
- summary1/summarize_cohort_var_count_print.py
- summary1/summarize_sample_var_count.py
- summary1/summarize_sample_var_count_print.py

## AF, annotation stats
- summary2/summary_freq_novel_func.py
- summary2/print_freq_novel.py
- summary2/print_freq_novel_func.py
