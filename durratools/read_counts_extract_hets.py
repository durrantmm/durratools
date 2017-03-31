import sys

def extract_hets(input_read_counts_stream, input_vcf_stream, out_stream, vcf_sample_name):
    hets = set()
    for rec in input_vcf_stream:
        gt = rec.genotype(vcf_sample_name)
        print(gt)