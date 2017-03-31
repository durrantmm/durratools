import sys

def extract_hets(input_read_counts_stream, input_vcf_stream, out_stream):
    for rec in input_vcf_stream:
        out_stream.write(str(rec))