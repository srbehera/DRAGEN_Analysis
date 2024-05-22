import sys
import json


def summary_freq_novel_func(in_folder_path, nscut):

    shard_list = [f'shard-{i}' for i in range(1, 101)]
    data = {
        'counter': {},
        'counter2': {},
    }
    for shard in shard_list:
        json_fname = f'{in_folder_path}/{shard}/dragen.freq-novel-func.nscut-{nscut}.json'
        print(f'reading {json_fname}')
        _data = json.loads(open(json_fname).read())
        for kk in ['counter', 'counter2']:
            for k, v in _data[kk].items():
                if k not in data[kk]:
                    data[kk][k] = 0
                data[kk][k] += v

    out_fname = f'summary_freq_novel_func.nscut-{nscut}.json'
    print(f'writing {out_fname}')
    with open(out_fname, 'w') as f:
        print(json.dumps(data, indent=4), file=f)


in_folder_path = '../../data/var_stats'
nscut = sys.argv[1]
summary_freq_novel_func(in_folder_path, nscut)
