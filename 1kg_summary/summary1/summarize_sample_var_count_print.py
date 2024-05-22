import json
import os


def print_tsv():

    data = dict()
    data['DRAGEN v4.2'] = json.loads(
        open('summarize_sample_var_count.json').read())

    call_list = ['DRAGEN v4.2']

    for pop in ['all', 'afr', 'eur', 'sas', 'eas', 'amr']:
        # header
        out = [pop.upper()] + call_list
        print('\t'.join(out))

        for vtype in ['loci', 'alleles', 'snv', 'ins', 'del', 'mnv']:
            # print data
            out = [vtype]
            for call in call_list:
                value = data[call][f'{pop}'][vtype]
                out.append(str(value))
            print('\t'.join(out))
        print()


print_tsv()
