import os
import sys
import json


def summarize_var_count(in_folder_path, out_json_path):

    summary_dict = dict()

    pop_list = ['all', 'afr', 'eur', 'sas', 'eas', 'amr']
    shard_list = [f'shard-{i}' for i in range(1, 101)]
    for shard in shard_list:
        fname = f'{in_folder_path}/{shard}/dragen.count_var_cohort.json'
        print(f'Reading {fname}')
        data = json.loads(open(fname).read())
        for k, v in data.items():
            if k not in summary_dict:
                summary_dict[k] = {}
            for kk, vv in v.items():
                if kk not in summary_dict[k]:
                    summary_dict[k][kk] = 0
                summary_dict[k][kk] += vv

    print(f'Writing output to {out_json_path}')
    os.makedirs(os.path.dirname(
        os.path.realpath(out_json_path)), exist_ok=True)
    with open(out_json_path, 'w') as f:
        print(json.dumps(summary_dict, indent=4), file=f)


in_folder_path, out_json_path = sys.argv[1:3]
summarize_var_count(in_folder_path, out_json_path)
