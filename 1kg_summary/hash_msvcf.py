import sys
import json
import subprocess
import os


def hash_msvcf(msvcf_url, out_vcf, region=None):

    if out_vcf.endswith('.gz'):
        out_vcf = out_vcf[:-3]

    fout = open(out_vcf, 'w')
    if region == None:
        cmd = f'zcat {msvcf_url}'
    else:
        cmd = f'bcftools view --with-header --regions {region} --regions-overlap 0 {msvcf_url}'
    proc = subprocess.Popen(cmd, shell=True, universal_newlines=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    print(f'Running {cmd}')
    print(f'Writing {out_vcf}')
    sample_list = []
    for line in proc.stdout:
        line = line.strip()
        if line.startswith('##'):
            print(line, file=fout)
            continue
        data = line.strip().split()
        if line.startswith('#'):
            sample_list = data[9:]
            # ##DRAGENCommandLine=<ID=
            print(
                f'##SampleInfo=<SampleCount={len(sample_list)},SampleIdList="{",".join(sample_list)}">', file=fout)
            print(f'##INFO=<ID=GTMAP,Number=.,Type=String,Description="GT sample ids map (ID is 0 based, excluded GT=0,0/0,0|0)">', file=fout)
            print('\t'.join(data[:8]), file=fout)
            continue

        # process GT
        fmt_dict = {x: i for i, x in enumerate(data[8].split(':'))}
        gt_map = {}
        for i, x in enumerate(data[9:]):
            gt = x.split(':')[0]
            if gt in {'0/0', '0|0'}:
                continue

            # ----- update GT when sample has no data -----
            # no data can not be distinguished from GT: . or ./. can be no data or hom-ref
            # use DP/GQ to indicate no data, when DP=0 or GQ="." no data
            # other fields are either not common among dragen and gatk versions, or allele dependent
            # dragen v4.0.3  GT:GQ:AD:FT:LPL:LAA:QL
            # dragen v3.7.8  GT:FT:GQ:DP:AD:SB:GP:PL:PRI
            # dragen v3.5.7b GT:FT:GQ:DP:AD:SB:GP:PL:PRI

            # gatk 1kg 3202  GT:AD:DP:GQ:PGT:PID:PL
            # gatk 1kg 2504  GT:AB:AD:DP:GQ:PL
            # note gatk release 2504 has variable number of fields per sample and not matching FORMAT field (see
            # tabix down/CCDG_13607_B01_GRM_WGS_2019-02-19_chr2.recalibrated_variants.vcf.gz chr2:308905-308905 | grep 308905
            # FORMAT = GT:AB:AD:DP:GQ:PGT:PID:PL, but some sample has ./.:.:0,0,0,0,0

            fields = x.split(':')

            # fix GQ=. dragen chrM
            if data[0] == 'chrM' and fields[fmt_dict['GQ']] == '.':
                fields[fmt_dict['GQ']] = '30'

            try:
                if 'DP' in fmt_dict:
                    if fields[fmt_dict['DP']] in {'0', '.'}:
                        gt = 'X'
                elif 'GQ' in fmt_dict:
                    if fields[fmt_dict['GQ']] == '.':
                        gt = 'X'
            except:
                # failed to get DP and GQ from the field due to GATK FORMAT field bug, treated as no data
                gt = 'X'

            if gt not in gt_map:
                gt_map[gt] = []
            gt_map[gt].append(i)

        gt_map_str = ",".join(
            [f'{gt}?{":".join([str(s) for s in slist])}' for gt, slist in gt_map.items()])
        out = data[:8]
        out[7] += f';GTMAP={gt_map_str}'
        print('\t'.join(out), file=fout)

    fout.close()

    stderr_list = []
    for x in proc.stderr:
        stderr_list.append(x)
        raise Exception(f'WARNING: non-empty stderr: ', ''.join(stderr_list))

    cmd = f'bgzip -f {out_vcf} && tabix -f {out_vcf}.gz'
    print(cmd)
    os.system(cmd)


msvcf_url, out_vcf = sys.argv[1:3]
hash_msvcf(msvcf_url, out_vcf)
