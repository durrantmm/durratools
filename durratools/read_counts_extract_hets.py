import sys

def extract_hets(input_read_counts_stream, input_vcf_stream, vcf_sample_name, out_stream):

    hets = set()
    counter = 1
    for rec in input_vcf_stream:
        gt = rec.genotype(vcf_sample_name)['GT']
        if gt in ['0/1', '1/0', '0|1', '1/0']:
            hets.add((rec.CHROM, rec.POS))
            print(counter, len(hets))
        counter += 1

    header = input_read_counts_stream.readline().strip().split()
    for line in input_read_counts_stream:
        line = line.strip().split()
        line = {header[i]: line[i] for i in range(len(line))}
        print(line)