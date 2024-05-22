import json
import os


def print_tsv():

    data = dict()

    data['DRAGEN v4.2'] = json.loads(
        open('summarize_cohort_var_count.json').read())

    call_list = ['DRAGEN v4.2']

    for pop in ['all', 'afr', 'eur', 'sas', 'eas', 'amr']:
        # header
        out = [pop.upper()] + call_list
        print('\t'.join(out))

        for vtype in ['loci', 'alleles', 'snv', 'ins', 'del', 'mnv']:
            # print data
            out = [vtype]
            for call in call_list:
                if vtype != 'alleles':
                    value = data[call][f'1kg3202.{pop}'][vtype]
                else:
                    value = sum([data[call][f'1kg3202.{pop}'][v] for v in [
                                'snv', 'ins', 'del', 'mnv']])
                out.append(str(value))
            print('\t'.join(out))
        print()


print_tsv()
