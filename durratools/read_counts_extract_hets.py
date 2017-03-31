import sys

def extract_hets(input_read_counts_stream, input_vcf_stream, vcf_sample_name, out_stream):

    hets = set()
    counter = 1
    for rec in input_vcf_stream:
        gt = rec.genotype(vcf_sample_name)['GT']
        if gt in ['0/1', '1/0', '0|1', '1/0']:
            hets.add((rec.CHROM, rec.POS))
        counter += 1

    header = input_read_counts_stream.readline()
    for line in input_read_counts_stream:
        line = line.strip().split()
        chrom = line[0]
        pos = line[1]
        if (chrom, pos) in hets:
            out_stream.write('\t'.join(line)+'\n')