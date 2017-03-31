import sys

def extract_hets(input_read_counts_stream, input_vcf_stream, vcf_sample_name, out_stream):

    hets = set()
    for rec in input_vcf_stream:
        gt = rec.genotype(vcf_sample_name)['GT']
        if gt in ['0/1', '1/0', '0|1', '1/0']:
            hets.add(gt)

    for line in input_read_counts_stream:
        line = line.strip().split()
        print(line)