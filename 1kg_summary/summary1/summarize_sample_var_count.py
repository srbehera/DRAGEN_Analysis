import os
import sys
import json


def load_sample_cohort(cohort_def_list):

    sample_cohort_map = dict()
    pop_sample_size = dict()
    for line in open(cohort_def_list):
        line = line.strip()
        pop, count, slist = line.split('\t')
        if not (pop.startswith('1kg3202.') and pop.endswith('.both')):
            continue
        pop = pop.split('.')[1]
        if pop not in pop_sample_size:
            pop_sample_size[pop] = 0
        for s in slist.split(','):
            if s not in sample_cohort_map:
                sample_cohort_map[s] = []
            sample_cohort_map[s].append(pop)
            pop_sample_size[pop] += 1

    for s in sample_cohort_map:
        if len(sample_cohort_map[s]) != 2:
            raise Exception(
                f'sample {s} has more than two pop labels: {sample_cohort_map[s]}')

    return sample_cohort_map, pop_sample_size


def summarize_var_count(in_folder_path, out_json_path):

    # load sample cohort tag
    sample_cohort_map, pop_sample_size = load_sample_cohort(
        '../metadata/cohort_def.list')
    summary_dict = dict()

    pop_list = ['all', 'afr', 'eur', 'sas', 'eas', 'amr']
    shard_list = [f'shard-{i}' for i in range(1, 101)]
    for shard in shard_list:
        fname = f'{in_folder_path}/{shard}/dragen.count_var_sample.json'
        print(f'Reading {fname}')
        data = json.loads(open(fname).read())
        for s, v in data.items():
            for p in sample_cohort_map[s]:
                if p not in summary_dict:
                    summary_dict[p] = {}
                for kk, vv in v.items():
                    if kk not in summary_dict[p]:
                        summary_dict[p][kk] = 0
                    summary_dict[p][kk] += vv

    # normalize by sample count
    for p in summary_dict:
        for kk in summary_dict[p]:
            summary_dict[p][kk] = int(
                summary_dict[p][kk]/pop_sample_size[p] + 0.5)

    print(f'Writing output to {out_json_path}')
    os.makedirs(os.path.dirname(
        os.path.realpath(out_json_path)), exist_ok=True)
    with open(out_json_path, 'w') as f:
        print(json.dumps(summary_dict, indent=4), file=f)


in_folder_path, out_json_path = sys.argv[1:3]
summarize_var_count(in_folder_path, out_json_path)
