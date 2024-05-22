import json
import os
from os.path import dirname, realpath
import sys
import gzip


def get_vtype(ref, alt):

    if len(ref) == 1 and len(alt) == 1:
        return 'snv'
    elif len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:
        return 'del'
    elif len(ref) == 1 and len(alt) > 1 and ref[0] == alt[0]:
        return 'ins'
    else:
        return 'mnv'


def get_allele_freq_novel_func(in_fname, out_vcf_fname, out_json_fname, pass_only=False, ns_cut=0.95):

    # chr5.gtc-anno.vcf.gz
    if out_vcf_fname.endswith('.gz'):
        out_vcf_fname = out_vcf_fname[:-3]
    os.makedirs(dirname(realpath(out_vcf_fname)), exist_ok=True)

    print(f'Reading {in_fname}')
    print(f'Writing {out_vcf_fname}')

    fout = open(out_vcf_fname, 'w')
    print(f'##fileformat=VCFv4.2', file=fout)
    print(f'##INFO=<ID=VarType,Number=1,Type=String,Description="Variant Type: snv, del, ins, mnv">', file=fout)
    print(f'##INFO=<ID=Novelty,Number=1,Type=String,Description="Novelty compared to dbSNP v155: known, novel">', file=fout)
    print(f'##INFO=<ID=FreqBin,Number=1,Type=String,Description="Allele frequency bins: singleton, rare, intmed, common">', file=fout)
    print(f'##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">', file=fout)
    print(f'##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">', file=fout)
    print(f'##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">', file=fout)
    print('\t'.join(
        f'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO'.split()), file=fout)

    counter = {}
    counter2 = {}
    nline = 0
    for line in gzip.open(in_fname):
        line = line.decode().strip()
        if line.startswith('#'):
            continue
        data = line.split('\t')
        data[2] = '.'

        if pass_only and data[6] != 'PASS':
            new_info = {'FILTERED': 'SITE_NON_PASS'}
            data[7] = json.dumps(new_info).replace(' ', '')
            vcf_line = data[:8]
            print('\t'.join(vcf_line), file=fout)
            continue

        info = json.loads(data[-1])

        cohort = '1kg2504.all.both' if data[0] != 'chrY' else '1kg2504.all.male'
        if data[0] != 'chrY':
            n = sum(info['gtc'][cohort].values())
            if n != 2504:
                raise Exception(f'*ERROR* chr1-22,X NS: {n}!=2504: {line}')
        else:
            n = sum(info['gtc'][cohort].values())
            if n != 1233:
                raise Exception(f'*ERROR* chrY NS: {n}!=1233: {line}')

        ns = 0
        an = 0
        alt_list = data[4].split(',')
        nalt = len(alt_list)
        alt_count = [0 for i in range(nalt)]
        for gt, cnt in info['gtc'][cohort].items():
            gt = gt.replace('|', '/')
            if gt in {'./.', '.|.', '.', 'X'}:
                continue
            ns += cnt
            if '/' in gt:
                an += cnt*2
            else:
                an += cnt
            for i in gt.split('/'):
                i = int(i)
                if i == 0:
                    continue
                alt_count[i-1] += cnt

        if an == 0 or ns == 0:
            # no sample has GT
            new_info = {'FILTERED': 'SITE_ZERO_NS_AN_GT',
                        'NS_GT': ns, 'AN_GT': an}
            data[7] = json.dumps(new_info).replace(' ', '')
            vcf_line = data[:8]
            print('\t'.join(vcf_line), file=fout)
            continue

        # filter variants with GT/data < ns_cut (e.g. 0.50 or 0.95)
        if data[0] != 'chrY' and ns < 2504 * ns_cut:
            # print(f'Skip low ns line {ns}: {line}')
            new_info = {'FILTERED': 'SITE_LOW_GT_RATE', 'NS_GT': ns}
            data[7] = json.dumps(new_info).replace(' ', '')
            vcf_line = data[:8]
            print('\t'.join(vcf_line), file=fout)
            continue

        if data[0] == 'chrY' and ns < 1233 * ns_cut:
            # print(f'Skip low ns line {ns}: {line}')
            new_info = {'FILTERED': 'SITE_LOW_GT_RATE', 'NS_GT': ns}
            data[7] = json.dumps(new_info).replace(' ', '')
            vcf_line = data[:8]
            print('\t'.join(vcf_line), file=fout)
            continue

        # loop per alt allele
        for i, alt in enumerate(data[4].split(',')):
            if alt == '*' or alt == '<NON_REF>':
                new_info = {'FILTERED': 'ALT_NON_REF_OR_STAR'}
                data[4] = alt
                data[7] = json.dumps(new_info).replace(' ', '')
                vcf_line = data[:8]
                print('\t'.join(vcf_line), file=fout)
                continue

            if alt not in info['anno']:
                raise Exception(
                    f'*ERROR* alt {alt} not found in annotation. Line = {line}')

            ac = alt_count[i]

            if ac == 0:
                # no sample has this alt allele
                new_info = {'FILTERED': 'ALT_ZERO_AC', 'AC': ac}
                data[4] = alt
                data[7] = json.dumps(new_info).replace(' ', '')
                vcf_line = data[:8]
                print('\t'.join(vcf_line), file=fout)
                continue

            new_info = {}

            # sample with alt allele
            if ac == 1:
                freq = 'singleton'
            elif ac <= an*0.01:
                # singleton is also rare, but flagged separately
                freq = 'rare'
            elif ac <= an*0.05:
                freq = 'intmed'
            else:
                freq = 'common'
            af = ac/an

            new_info['AC'] = ac
            new_info['AN'] = an
            new_info['AF'] = round(af, 6)
            new_info['FreqBin'] = freq

            # normalize function annotation when multiple is available
            func = info['anno'][alt]['func']
            if func == {}:
                func = {'unknown': 1}
            new_info['FuncAnno'] = func

            # update key with normalize/left-aligned key
            known = False
            vkey = f'{data[0]}-{data[1]}-{data[3]}-{alt}'
            vid = info['anno'][alt]['vid']
            if vid != '.':
                # variant is known or needs normalization
                for v in vid.split(','):
                    if v.startswith('rs'):
                        known = True
                    else:
                        vkey = v

            c, p, r, a = vkey.split('-')
            vtype = get_vtype(r, a)
            known = 'known' if known == True else 'novel'

            new_info['Novelty'] = known
            new_info['VarType'] = vtype

            # bin by annotation (allow double counting per alt allele)
            tag = f'{vtype}:{known}:{freq}'
            new_info['Tag'] = tag

            # counter has no double counting, it is per alt
            if tag not in counter:
                counter[tag] = 0
            counter[tag] += 1

            # counter2 has double counting per alt due to multiple annotation per alt
            for k, v in func.items():
                tag2 = tag + f':{k}'
                if tag2 not in counter2:
                    counter2[tag2] = 0
                counter2[tag2] += 1

            data[4] = alt
            data[7] = json.dumps(new_info).replace(' ', '')
            vcf_line = data[:8]
            print('\t'.join(vcf_line), file=fout)

        # end loop per allele

    cmd = f'bgzip -f {out_vcf_fname} && tabix -f {out_vcf_fname}.gz'
    print(f'Running {cmd}')
    os.system(cmd)

    print(f'Writing {out_json_fname}')
    os.makedirs(dirname(realpath(out_json_fname)), exist_ok=True)
    with open(out_json_fname, 'w') as f:
        print(json.dumps({
            'counter': counter,
            'counter2': counter2
        }, indent=4, sort_keys=True), file=f)


in_fname, out_vcf_fname, out_json_fname, pass_only, ns_cut = sys.argv[1:6]
pass_only = True if pass_only == 'true' else False
ns_cut = float(ns_cut)
get_allele_freq_novel_func(in_fname, out_vcf_fname,
                           out_json_fname, pass_only, ns_cut)
