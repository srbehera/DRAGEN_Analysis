import os
import sys
import gzip
import json


def get_anno(input_json_path, output_vcf_path):

    if output_vcf_path.endswith('.gz'):
        output_vcf_path = output_vcf_path[:-3]

    print(f'Reading {input_json_path}')
    os.makedirs(os.path.dirname(
        os.path.realpath(output_vcf_path)), exist_ok=True)

    print(f'Writing {output_vcf_path}')

    fout = open(output_vcf_path, 'w')
    print(f'##fileformat=VCFv4.2', file=fout)
    print(f'##INFO=<ID=FA,Number=.,Type=String,Description="Nirvana functional annotation consequences">', file=fout)
    print('\t'.join(f'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO'.split()),
          file=fout)

    for line in gzip.open(input_json_path):
        line = line.decode().strip()
        if not line.startswith('{"chromosome"'):
            continue

        line = line.rstrip(',')
        # Nirvana bug at chr2:20447683: "phenotypes":[""See Cases""],
        line = line.replace('""', '"')
        data = json.loads(line)

        alt_alleles = [x for x in data["altAlleles"]
                       if x not in {'<NON_REF>', '*'}]
        nalt = len(alt_alleles)
        nvar = len(data["variants"])
        if nalt != nvar:
            out = {
                'issue': 'altAlleles and variants field length mismatch',
                'data': data
            }
            print(json.dumps(out))
            exit(0)

        # write vcf line per alt allele
        for alt, var_anno in zip(alt_alleles, data["variants"]):
            # _ref, _alt = normalize_var1(data["refAllele"], alt)
            qual = '.'
            filt = '.'

            conseq_dict = dict()
            if 'transcripts' in var_anno:
                for x in var_anno['transcripts']:
                    for y in x['consequence']:
                        if y not in conseq_dict:
                            conseq_dict[y] = 0
                        conseq_dict[y] += 1

            conseq_str = []
            for conseq, count in sorted(conseq_dict.items(), key=lambda x: x[1], reverse=True):
                conseq_str.append(f'{conseq}:{count}')
            conseq_str = ','.join(conseq_str) if conseq_str else '.'
            info = f'FA={conseq_str}'

            rsid_set = set()
            if 'dbsnp' in var_anno:
                for x in var_anno['dbsnp']:
                    rsid_set.add(x)

            if 'chr'+var_anno['vid'] != f'{data["chromosome"]}-{data["position"]}-{data["refAllele"]}-{alt}':
                rsid_set.add(var_anno['vid'])

            rsid = ','.join(sorted(rsid_set)) if rsid_set else '.'

            out = [
                data['chromosome'], str(data['position']), rsid,
                data["refAllele"], alt,
                qual, filt, info
            ]
            out = '\t'.join(out)

            # if conseq_dict:
            #     print(out)
            print(out, file=fout)

    fout.close()

    print(f'Compressing and indexing output VCF: {output_vcf_path}')
    os.system(f'bgzip -f {output_vcf_path}')
    os.system(f'tabix -f {output_vcf_path}.gz')


if __name__ == '__main__':

    input_json_path, output_vcf_path = sys.argv[1:3]
    get_anno(input_json_path, output_vcf_path)
