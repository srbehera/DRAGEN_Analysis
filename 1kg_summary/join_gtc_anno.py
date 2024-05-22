import os
import sys
import json
import gzip
import subprocess


def load_contig_size(shard_fname, chrom):

    length = -1
    for line in open(shard_fname):
        data = line.split('\t')
        if not data[1].startswith(chrom + ':') or ',' in data[1]:
            continue

        pend = int(data[1].split('-')[-1])
        if pend > length:
            length = pend

    if length < 0:
        raise Exception(f'{chrom} is not found in {shard_fname}')

    return length


def load_anno(fname):

    cmd = f'zcat {fname}'
    print(f'Running {cmd}')
    proc = subprocess.Popen(cmd, shell=True, universal_newlines=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    anno_dict = dict()
    for line in proc.stdout:
        line = line.strip()
        if line[0] == '#':
            continue
        data = line.split('\t')
        chrom = data[0]
        pos = int(data[1])
        key = f'{data[0]}:{data[1]}:{data[3]}:{data[4]}'
        anno = dict()
        info = data[7].split('=')[1]
        if info != '.':
            for x in info.split(','):
                k, v = x.split(':')
                anno[k] = int(v)
        anno_dict[key] = {'vid': data[2], 'func': anno}

    return anno_dict


def join_gtc_anno(gtc_fname, anno_fname, out_fname):

    anno_dict = load_anno(anno_fname)

    os.makedirs(os.path.dirname(os.path.realpath(out_fname)), exist_ok=True)
    if out_fname.endswith('.gz'):
        out_fname = out_fname[:-3]

    cmd = f'zcat {gtc_fname}'
    print(f'Running {cmd}')
    proc = subprocess.Popen(cmd, shell=True, universal_newlines=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    fout = open(out_fname, 'w')
    print(f'##fileformat=VCFv4.2', file=fout)
    print('\t'.join(f'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO'.split()),
          file=fout)
    for line in proc.stdout:
        line = line.strip()
        if line[0] == '#':
            continue
        data = line.split('\t')
        chrom = data[0]
        pos = int(data[1])
        anno = dict()
        for alt in data[4].split(','):
            if alt in {'*', '<NON_REF>'}:
                continue
            key = f'{data[0]}:{data[1]}:{data[3]}:{alt}'
            if key not in anno_dict:
                raise (f'key = {key} not found in annotation {anno_fname}')
            anno[alt] = anno_dict[key]

        info = json.loads(data[7])
        info = {'anno': anno, 'gtc': info}
        info = json.dumps(info).replace(' ', '')

        data[7] = info
        print('\t'.join(data), file=fout)

    fout.close()

    print(f'Compressing and indexing output VCF: {out_fname}')
    os.system(f'bgzip -f {out_fname}')
    os.system(f'tabix -f {out_fname}.gz')


gtc_fname, anno_fname, out_fname = sys.argv[1:4]
join_gtc_anno(gtc_fname, anno_fname, out_fname)
