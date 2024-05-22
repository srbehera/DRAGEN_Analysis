import os
import sys
import json
import gzip


def load_gt_map(gtmap_str, nsample):

    gtmap_dict = dict()
    gt00_slist = {str(i) for i in range(nsample)}

    if gtmap_str != '':
        for x in gtmap_str.split(','):
            gt, slist = x.split('?')
            slist = slist.split(':')
            gtmap_dict[gt] = slist
            for s in slist:
                gt00_slist.remove(s)

    gtmap_dict['0/0'] = sorted(gt00_slist)
    return gtmap_dict


def subset_gtmap(sample_cohort_dict, gtmap_dict):

    cohort_gtmap_dict = dict()
    for gt, slist in gtmap_dict.items():
        for sample in slist:
            cohort_list = sample_cohort_dict[sample]
            for cohort in cohort_list:
                if cohort not in cohort_gtmap_dict:
                    cohort_gtmap_dict[cohort] = dict()
                if gt not in cohort_gtmap_dict[cohort]:
                    cohort_gtmap_dict[cohort][gt] = 0
                cohort_gtmap_dict[cohort][gt] += 1

    return cohort_gtmap_dict


def count_gt(input_hash_mscvf_path, cohort_def_path, out_vcf_path):

    print(f'Reading {cohort_def_path}')
    sample_cohort_dict = dict()
    with open(cohort_def_path) as f:
        for line in f:
            cohort_name, nsample, sample_ids = line.strip().split()
            for sample_id in sample_ids.split(','):
                if sample_id not in sample_cohort_dict:
                    sample_cohort_dict[sample_id] = []
                sample_cohort_dict[sample_id].append(cohort_name)

    if out_vcf_path.endswith('.vcf.gz'):
        out_vcf_path = out_vcf_path[:-3]

    print(f'Reading {input_hash_mscvf_path}')
    print(f'Writing {out_vcf_path}')
    nsample = 0
    os.makedirs(os.path.dirname(os.path.realpath(out_vcf_path)), exist_ok=True)
    fout = open(out_vcf_path, 'w')
    for line in gzip.open(input_hash_mscvf_path):
        line = line.decode().strip()
        if line.startswith('#'):
            if line.startswith('##SampleInfo'):
                sample_list = line.split('"')[-2].split(',')
                sample_index_dict = {s:str(i) for i,s in enumerate(sample_list)}
                nsample = len(sample_list)
                #     sample_cohort_dict = {sample_index_dict[s]:clist for s,clist in sample_cohort_dict.items()}

                sample_cohort_dict = {
                    sample_index_dict[s]:clist
                    for s,clist in sample_cohort_dict.items()
                    # only keep samples in the hashed msvcf
                    if s in sample_index_dict
                }

                if nsample == 2504:
                    # filter 1kg3202.* cohort
                    for s in sample_cohort_dict:
                        sample_cohort_dict[s] = [
                            c for c in sample_cohort_dict[s] if '2504' in c]

            if line.startswith('##fileformat=VCFv4.2') or line.startswith('#CHROM'):
                print(line, file=fout)
            continue

        chrom, pos, rsid, ref, alt, qual, filt, info = line.strip().split()[:8]

        gtmap_str = info.split(';')[-1].split('=')[-1]
        gtmap_dict = load_gt_map(gtmap_str, nsample)
        cohort_gtmap_dict = subset_gtmap(sample_cohort_dict, gtmap_dict)
        info = json.dumps(cohort_gtmap_dict).replace(' ','')
        out = '\t'.join([chrom, pos, rsid, ref, alt, qual, filt, info])
        print(out, file=fout)

    fout.close()

    print(f'Compressing and indexing output VCF: {out_vcf_path}')
    os.system(f'bgzip -f {out_vcf_path}')
    os.system(f'tabix -f {out_vcf_path}.gz')



input_hash_mscvf_path, cohort_def_path, out_vcf_path = sys.argv[1:4]
count_gt(input_hash_mscvf_path, cohort_def_path, out_vcf_path)
