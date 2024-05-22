import os
import sys
import json
import gzip


def load_gt_map(gtmap_str, index_sample_dict):

    gtmap_dict = dict()
    gt00_slist = {s for i,s in index_sample_dict.items()}

    if gtmap_str != '':
        for x in gtmap_str.split(','):
            gt, slist = x.split('?')
            slist = slist.split(':')
            gtmap_dict[gt] = sorted([index_sample_dict[i] for i in slist])
            for i in slist:
                gt00_slist.remove(index_sample_dict[i])

    gtmap_dict['0/0'] = sorted(gt00_slist)
    return gtmap_dict


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


def count_sample_var(input_hash_mscvf_path, out_json_path, pass_only=False):

    sample_var_stat = dict()
    print(f'Reading {input_hash_mscvf_path}')
    index_sample_dict = dict()
    nline = 0
    for line in gzip.open(input_hash_mscvf_path):
        # if nline == 1000:
        #     break

        line = line.decode().strip()
        if line.startswith('#'):
            if line.startswith('##SampleInfo'):
                sample_list = line.split('"')[-2].split(',')
                index_sample_dict = {str(i):s for i,s in enumerate(sample_list)}
            continue

        chrom, pos, rsid, ref, alt, qual, filt, info = line.strip().split()[:8]
        if pass_only == True and filt != 'PASS':
            continue

        gtmap_str = info.split(';')[-1].split('=')[-1]
        gtmap_dict = load_gt_map(gtmap_str, index_sample_dict)

        alist = [ref] + [a for a in alt.split(',')]

        for gt, slist in gtmap_dict.items():
            if len(slist) == 0 and gt != '0/0':
                print(f'Skip empty sample list with gt = {gt} at line: {line}')
                exit(0)
                continue
            gt = gt.replace('|', '/')
            if gt in {'X','.','./.','0','0/0'}:
                continue
            alt_index_set = set()
            for i in gt.split('/'):
                i = int(i)
                if i == 0 or  alist[i] in {'*', '<NON_REF>'}:
                    continue
                alt_index_set.add(alist[i])
            for alt in alt_index_set:
                # one of snv, ins, del, mnv
                vtype = get_var_type(ref, alt)
                for s in slist:
                    if s not in sample_var_stat:
                        sample_var_stat[s] = {
                            'snv':0, 'ins':0, 'del':0, 'mnv':0, 'alleles':0, 'loci':0}
                    sample_var_stat[s][vtype] += 1
                    sample_var_stat[s]['alleles'] += 1
            if len(alt_index_set):
                for s in slist:
                    sample_var_stat[s]['loci'] += 1

        nline += 1

    print(f'Writing {out_json_path}')
    os.makedirs(os.path.dirname(os.path.realpath(out_json_path)), exist_ok=True)
    fout = open(out_json_path, 'w')
    print(json.dumps(sample_var_stat, indent=4), file=fout)
    fout.close()


input_hash_mscvf_path, out_json_path, pass_only = sys.argv[1:4]
pass_only = bool(pass_only.lower() == 'true')
count_sample_var(input_hash_mscvf_path, out_json_path, pass_only)
