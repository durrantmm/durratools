import sys
import vcf
from Bio import SeqIO

reference_format = "fasta"
vcf_header_start = "#CHROM"

def correct_vcf(vcf_path, reference_path, out_vcf_path=None):
    if not out_vcf_path:
        out_vcf_path = vcf_path.split('.')[0] + '.CORRECTED.vcf'

    reference = SeqIO.parse(reference_path, reference_format)
    current_chrom = next(reference)

    vcf_reader = vcf.Reader(filename=vcf_path)

    for var in vcf_reader:

        if len(var.ALT) > 1:
            raise TypeError("The vcf has multi-allelic sites. NOT COOL.")

        current_chrom = _get_chrom(current_chrom, reference, var.CHROM)
        query_info = [var.CHROM, var.POS, var.ID, var.REF, str(var.ALT[0])]
        true_ref = _get_true_ref(current_chrom, var.POS)

        _adjust_variant_record_simple(var, query_info, true_ref)


def _adjust_variant_record_simple(variant_rec, query_info, true_ref):
    query_chrom, query_pos, query_id, query_ref, query_alt = query_info

    GT = variant_rec.samples[0]['GT']
    if query_ref == true_ref:

        print(_output_line(query_chrom, query_pos, query_id, query_ref, query_alt, GT, GT))

    elif query_alt == true_ref:
        print(_output_line(query_chrom, query_pos, query_id, query_alt, query_ref, _flip_genotype(GT), GT))

    elif query_ref == _complement_base(true_ref):
        query_ref = _complement_base(query_ref)
        query_alt = _complement_base(query_alt)
        print(_output_line(query_chrom, query_pos, query_id, query_ref, query_alt, GT, GT))

    elif query_alt == _complement_base(true_ref):
        query_ref = _complement_base(query_ref)
        query_alt = _complement_base(query_alt)
        print(_output_line(query_chrom, query_pos, query_id, query_alt, query_ref, _flip_genotype(GT), GT))

    else:
        print("UNEXPECTED")
        print(query_info, true_ref)


def _comparison_variants_are_equal(comp_info_rsid, comp_info_loc):
    if comp_info_rsid == comp_info_loc:
        return True
    else:
        return False


def _alleles_match(ref1, alt1, alleles):
    if ref1 in alleles and alt1 in alleles:
        return True
    else:
        return False

def _alleles_match_complement(ref1, alt1, alleles):
    if _complement_base(ref1) in alleles and _complement_base(alt1) in alleles:
        return True
    else:
        return False

def _positions_match(chrom1, pos1, chrom2, pos2):
    if chrom1 == chrom2 and pos1 == pos2:
        return True
    else:
        return False

def _complement_base(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    else:
        return base

def _interpret_genotype(gt):
    if gt in ['0/1', '1/0']:
        return "HET"
    elif gt == '1/1':
        return "HOM_ALT"
    elif gt == '0/0':
        return "HOM_REF"

def _flip_genotype(gt):
    if gt == '1/1':
        return '0/0'
    elif gt == '0/0':
        return '1/1'
    else:
        return gt

def _output_line(query_chrom, query_pos, query_id, query_ref, query_alt, gt, gt_orig):
    gt, gt_orig = _interpret_genotype(gt), _interpret_genotype(gt_orig)
    return '\t'.join(map(str, [query_chrom, query_pos, query_ref, query_alt, query_id, gt]))


def _get_chrom(current_chrom, reference, selected_chrom):
    if str(current_chrom.id) == str(selected_chrom):
        return current_chrom
    else:
        while str(current_chrom.id) != str(selected_chrom):
            current_chrom = next(reference)
        return current_chrom


def _get_true_ref(current_chrom, pos):
    return current_chrom.seq[pos-1].upper()


if __name__ == '__main__':
    vcf_path = sys.argv[1]
    snp_reference_path = sys.argv[2]
    genome_reference_path = sys.argv[3]

    correct_vcf(vcf_path, snp_reference_path, genome_reference_path)
