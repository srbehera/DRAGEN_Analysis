import json
import pysam
import truvari
from collections import Counter
import sys

vcf = pysam.VariantFile(sys.argv[1])
cnt = Counter()
for entry in vcf:
    sz = max(1, truvari.entry_size(entry))
    for sample in entry.samples:
        if truvari.entry_is_present(entry, sample):
            cnt[sample] += sz
print(json.dumps(cnt, indent=4))
