import sys
import json


def print_freq_novel_func(in_json, out_csv):

    print(f'reading {in_json}')
    data = json.loads(open(in_json).read())

    # "del:known:common:3_prime_UTR_variant": 12764,

    # resummarize keys by
    # del + ins = indel
    # all, rare + singleton, novel

    print(f'writing {out_csv}')
    fout = open(out_csv, 'w')

    for k in sorted(data['counter2']):
        v = data['counter2'][k]
        print(k, v, sep='\t', file=fout)

    fout.close()

    # print freq novel table
    delimit = '\t'
    cols = ['singleton', 'rare', 'intmed', 'common']
    print(delimit.join(['DRAGEN v4.2.7'] + cols), file=fout)

    for row in ['snv:known', 'snv:novel', 'indel:known', 'indel:novel']:
        line = [row]
        for col in cols:
            if row.startswith('snv:'):
                k = f'{row}:{col}'
                line.append(str(data['counter'][k]))
            else:
                nvar = 0
                novel = row.split(':')[1]
                for vtype in ['ins', 'del']:
                    k = f'{vtype}:{novel}:{col}'
                    nvar += data['counter'][k]
                line.append(str(nvar))
        print(delimit.join(line), file=fout)

    fout.close()


in_json, out_csv = sys.argv[1:3]
print_freq_novel_func(in_json, out_csv)
