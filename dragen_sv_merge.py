#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import functools
import itertools

import pysam
import truvari
from truvari.bench import pick_single_matches

class RepHandler():
    """
    This is a light wrapper on a pysam.VariantFile that will edit entries as needed for merging
    """
    def __init__(self, vcf_fn, header):
        self.vcf_fn = vcf_fn
        self.vcf = pysam.VariantFile(self.vcf_fn)
        self.header = header

    def fetch(self, chrom, start=None, end=None):
        """
        Gets inbetween pysam.VariantFile.fetch
        """
        for entry in self.vcf.fetch(chrom, start, end):
            entry.translate(self.header)
            m_lens = []
            m_tys = []
            # ignore monozygotic reference for now
            if not entry.alts:
                continue
            alt = entry.alts[0]
            length = int(alt[4:-1])
            svlen = length * len(entry.info["RU"])
            if svlen > entry.info["RL"]:
                svtype = "INS"
            else:
                svtype = "DEL"
            entry.info["SVLEN"] = svlen
            entry.info["SVTYPE"] = svtype
            yield entry

class SvHandler():
    """
    This is a light wrapper on a pysam.VariantFile that will edit entries as needed for merging
    """
    def __init__(self, vcf_fn, header):
        self.vcf_fn = vcf_fn
        self.vcf = pysam.VariantFile(self.vcf_fn)
        self.header = header

    def fetch(self, chrom, start=None, end=None):
        """
        Gets inbetween pysam.VariantFile.fetch
        """
        for entry in self.vcf.fetch(chrom, start, end):
            entry.translate(self.header)
            yield entry

class CnvHandler():
    """
    This is a light wrapper on a pysam.VariantFile that will edit entries as needed for merging
    """
    def __init__(self, vcf_fn, header):
        self.vcf_fn = vcf_fn
        self.vcf = pysam.VariantFile(self.vcf_fn)
        self.header = header

    def fetch(self, chrom, start=None, end=None):
        """
        Gets inbetween pysam.VariantFile.fetch
        """
        for entry in self.vcf.fetch(chrom, start, end):
            # ignore monozygotic reference for now
            if not entry.alts:
                continue
            entry.translate(self.header)
            entry.info["SVTYPE"] = entry.alts[0][1:-1]
            yield entry

class MultiVCFIter():
    """
    In order to not make assumptions about sort order, we need a class that will
    allow us to call pysam.VariantFile.fetch on multiple vcfs and in the same contig order
    """
    def __init__(self, contigs):
        """
        Contigs is a list of lists of contigs. All the input contigs
        """
        self.contigs = set()
        for ctgs_a, ctgs_b in itertools.combinations(enumerate([set(_) for _ in contigs]), 2):
            d = ctgs_a[1] - ctgs_b[1]
            if d:
                logging.warning("VCF #%d and #%d have %d contig header differences", ctgs_a[0], ctgs_b[0], len(d))
            d = ctgs_b[1] - ctgs_a[1]
            if d:
                logging.warning("VCF #%d and #%d have %d contig header differences", ctgs_b[0], ctgs_a[0], len(d))

        self.contigs = list(set(itertools.chain(*contigs)))

    def iterate(self, vcf_file):
        """
        Given a VariantFile (or *Handler of one), fetch entries in order
        """
        for chrom in sorted(self.contigs):
            try:
                for entry in vcf_file.fetch(chrom):
                    yield entry
            except ValueError:
                # We've already warned if some contigs are vcf specific
                pass

def consolidate_genotypes(entry, other):
    """
    Consolidate the genotype field from `other` VariantRecords
    into the entry VariantRecord.
    Todos: 
    - Consolidate multiple format fields
    - Consolidate with multiple others (refactor truvari to keep it dry)
    - Add Hemi for partial missing..?
    """
    replace_gts = ["UNK", "REF", "NON"]
    for sample in entry.samples:
        m_gt = truvari.get_gt(entry.samples[sample]["GT"]).name
        o_gt = truvari.get_gt(other.samples[sample]["GT"]).name
        if m_gt in replace_gts and o_gt not in replace_gts:
            entry.samples[sample]["GT"] = other.samples[sample]["GT"]

def null_collapser(match):
    """
    The most basic collapser. Just sets that the base variant has been merged.
    Will return the VCF entry that should be written
    """
    match.base.info["TruMerged"] = True
    consolidate_genotypes(match.base, match.comp)
    return match.base

def collapse_calls(base, comp, chunk_id, bencher, collapser=null_collapser):
    """
    Use the bencher to find variants in comp that should be collapsed into base

    When we get to actually consolidating information, a collapser will be run to
    put comp information into the base. For now, the collapser will do nothing

    Returns the new lists of variants
    """
    match_matrix = bencher.build_matrix(base, comp, chunk_id)
    ret_b = []
    ret_c = []
    # intramerge only needs pick single, inter will need its own picking logic
    for match in pick_single_matches(match_matrix):
        if not match.state:
            if match.base:
                ret_b.append(match.base)
            if match.comp:
                ret_c.append(match.comp)
        if match.state:
            ret_b.append(collapser(match))
    return ret_b, ret_c

def compare_chunk(chunk, rep_sv_bench, sv_cnv_bench):
    """
    Compares all calls within a chunk and yields entries as they should be made
    """
    calls, chunk_id = chunk

    # Since SV are the middle, without them in a chunk we
    # just pump them all forward
    # For intermerge, we'll need to do a sv/sv comparison
    if not calls['sv']:
        return itertools.chain(calls['rep'], calls['cnv'])

    if calls['rep'] and calls['sv']:
        calls['rep'], calls['sv'] = collapse_calls(calls['rep'], calls['sv'], chunk_id, rep_sv_bench)

    # CNV will also need a self comparison when this becomes intermerge
    if calls['sv'] and calls['cnv']:
        calls['sv'], calls['cnv'] = collapse_calls(calls['sv'], calls['cnv'], chunk_id, sv_cnv_bench)

    return itertools.chain(*calls.values())

def consolidate_header(rep_fn, sv_fn, cnv_fn):
    """
    Consolidate header information. Return header
    """
    rep_header = pysam.VariantFile(rep_fn).header
    sv_header = pysam.VariantFile(sv_fn).header
    cnv_header = pysam.VariantFile(cnv_fn).header

    def consol(into, outof):
        for i in outof.records:
            if not ((i.type == "INFO" and i["ID"] in into.info.keys()) or
                    (i.type == "FORMAT" and i["ID"] in into.formats.keys())):
                into.add_line(str(i))

    new_header = rep_header.copy()
    consol(new_header, sv_header)
    consol(new_header, cnv_header)
    return new_header

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="dragen_truvari_merge", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("repeats", metavar="REP", type=str,
                        help="repeats VCF file")
    parser.add_argument("sv", metavar="SV", type=str,
                        help="sv VCF file")
    parser.add_argument("cnv", metavar="CNV", type=str,
                        help="sv VCF file")
    # Do we want to priortize piping or compression/indexing of output?
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output VCF file")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args

def check_params(args):
    """
    Make sure we can merge these files
    """
    check_fail = False
    seen_samples = None
    for fn in [args.repeats, args.sv, args.cnv]:
        if not os.path.exists(fn):
            logging.error("File %s does not exist", fn)
            check_fail = True
        if not fn.endswith(".gz"):
            logging.error("VCF %s does not end with .gz. Must be bgzip'd", fn)
            check_fail = True
        if not os.path.exists(fn + '.tbi'):
            logging.error("VCF index %s.tbi does not exist. Must be indexed", fn)
            check_fail = True
        v = pysam.VariantFile(fn)
        if seen_samples is None:
            seen_samples = list(v.header.samples)
        elif len(seen_samples) != len(v.header.samples):
            logging.error("VCFs have differing number of samples")
            check_fail = True
        elif [_ for _ in zip(seen_samples, v.header.samples) if _[0] != _[1]]:
            logging.error("VCFs sample columns do not line up")
            check_fail = True
            
    return check_fail

def dragen_merge_main(args):
    """
    Main entrypoint
    """
    args = parse_args(args)
    if check_params(args):
        logging.error("Unable to merge. Please fix parameters")
        sys.exit(1)

    header = consolidate_header(args.repeats, args.sv, args.cnv)
    header.add_line(('##INFO=<ID=TruMerged,Number=0,Type=Flag,'
                     'Description="Variant was merged by Truvari">'))

    rep = RepHandler(args.repeats, header)
    sv = SvHandler(args.sv, header)
    cnv = CnvHandler(args.cnv, header)

    m_iter = MultiVCFIter([rep.vcf.header.contigs, sv.vcf.header.contigs, cnv.vcf.header.contigs])

    matcher = truvari.Matcher()
    matcher.params.sizemin = 0
    matcher.params.sizefilt = 0
    matcher.params.sizemax = sys.maxsize
    matcher.params.check_monref = True # For now, we're going to be dropping these sites
    matcher.params.check_multi = True # For now, expecting split multi-allelics

    m_chunks = truvari.chunker(matcher,
        ('rep', m_iter.iterate(rep)),
        ('sv', m_iter.iterate(sv)),
        ('cnv', m_iter.iterate(cnv))
    )

    rep_sv_matcher = truvari.Matcher()
    rep_sv_matcher.params.pctseq = 0
    rep_sv_matcher.params.pctsize = 0
    rep_sv_matcher.params.refdist = 0
    rep_sv_matcher.params.sizemin = 0
    rep_sv_matcher.params.sizemax = 500
    rep_sv_bench = truvari.Bench(rep_sv_matcher)

    sv_cnv_matcher = truvari.Matcher()
    sv_cnv_matcher.params.pctseq = 0
    sv_cnv_matcher.params.pctsize = 0.70
    sv_cnv_matcher.params.sizemin = 800
    sv_cnv_matcher.params.sizefilt = 800
    sv_cnv_matcher.dup_to_ins = True
    sv_cnv_matcher.params.sizemax = sys.maxsize
    sv_cnv_bench = truvari.Bench(sv_cnv_matcher)

    m_compare = functools.partial(compare_chunk,
                                  rep_sv_bench=rep_sv_bench,
                                  sv_cnv_bench=sv_cnv_bench)

    out = pysam.VariantFile(args.output, 'w', header=header)

    num_merged = 0
    num_output = 0
    for variant in itertools.chain.from_iterable(map(m_compare, m_chunks)):
        # Eventually would like to make a counter object that will check
        # Multiple info fields that describe the changes made, collect them
        # Over all variants, and then can be output. For now, its just a single int
        if variant.info["TruMerged"]:
            num_merged += 1
        num_output += 1
        out.write(variant)
    logging.info("Merged %d calls", num_merged)
    logging.info("Output %d calls", num_output)
    logging.info("Finished")

if __name__ == '__main__':
    dragen_merge_main(sys.argv[1:])
