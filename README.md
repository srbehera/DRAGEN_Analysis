# DRAGEN_Analysis
BED files
1. CMRG_5030.bed - This file contains all 5,030 challenging medically relevant gene (CMRG) regions  
2. FixItFelix.collapsed.sorted.bed	- This contains the regions that are impacted due to collapsed errors in the GRCh38 reference
3. FixItFelix.duplicated.sorted.bed -  This contains the regions that are impacted due to duplicated errors in the GRCh38 reference
4. FixItFelix_12_CMRG.sorted.bed - This file contains the CMRG regions that are impacted due to either collapsed or duplicated errors

Counting basepairs
1. count_basepair.py - This script counts the total impacted basepairs of all variants at the sample level. It can be run with Tuvari-4.1-dev (https://github.com/ACEnglish/truvari)
2. count_bp_SNV_INDEL_chr21.json - Sample output file of the script
