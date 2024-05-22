import sys
import json
import os
import gzip

# count var by type: snp, ins, del, by pop
# number of samples = 3202 (for discovery not allele frequency)
# callset = dragen v4.0.3 v3.7.8 nygc2020

def get_var_type(ref, alt):
    if alt == '<NON_REF>' or alt == '*':
        raise Exception(f'can not determine variant type for alt={alt}')

    if len(ref) == 1 and len(alt) == 1:
        return 'snv'
    if len(ref) == len(alt) and ref[1:] == alt[1:]:
        return 'snv'

    if len(alt) > len(ref) and alt.endswith(ref[1:]) and ref[0] == alt[0]:
        return 'ins'
    if len(ref) > len(alt) and ref.endswith(alt[1:]) and ref[0] == alt[0]:
        return 'del'

    return 'mnv'


def count_var(in_vcf_path, out_json_path, cohort_name='1kg3202', pass_only=False):

    pop_list = [f'{cohort_name}.{x}'
                for x in ['all', 'afr', 'eur', 'eas', 'sas', 'amr']]
    vtype_list = ['snv', 'ins', 'del', 'mnv', 'loci']
    count_dict = {p:{v:0 for v in vtype_list} for p in pop_list}

    print(f'Reading {in_vcf_path}')
    for line in gzip.open(in_vcf_path):
        line = line.decode().strip()
        if line.startswith('#'):
            continue
        data = line.strip().split('\t')

        # filter non-pass
        if pass_only == True and data[6] != 'PASS':
            continue

        gtc = json.loads(data[7])
        alist = [data[3]] + data[4].split(',')

        has_mnv = False
        for pop in pop_list:
            cohort = f'{pop}.both'
            alt_index_set = set()
            for gt, count in gtc[cohort].items():
                gt = gt.replace('|', '/')
                if gt in {'X', '.', './.', '0', '0/0'}:
                    continue
                if count == 0:
                    continue
                for i in gt.split('/'):
                    i = int(i)
                    if i == 0 or alist[i] == '*' or alist[i] == '<NON_REF>':
                        continue
                    alt_index_set.add(i)
            for i in alt_index_set:
                vtype = get_var_type(alist[0], alist[i])
                count_dict[pop][vtype] += 1
                if vtype == 'mnv':
                    has_mnv = True

            if len(alt_index_set):
                count_dict[pop]['loci'] += 1
            else:
                if pop == '1kg3202.all':
                    print(f'*WARNING*: This line contains no variant in {pop}: ', line)
        if has_mnv:
            print(f'*WARNING*: This line contains MNV: ', line)

    print(f'Writing {out_json_path}')
    os.makedirs(os.path.dirname(os.path.realpath(out_json_path)), exist_ok=True)
    with open(out_json_path, 'w') as f:
        print(json.dumps(count_dict, indent=4), file=f)


in_vcf_path, out_json_path, cohort_name, pass_only = sys.argv[1:5]
pass_only = bool(pass_only.lower() == 'true')
count_var(in_vcf_path, out_json_path, cohort_name, pass_only)
