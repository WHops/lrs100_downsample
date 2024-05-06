###
import cyvcf2
import sys
import csv
import pandas as pd
from pathlib import Path
import pdb
import glob
import warnings
import json
from collections import Counter

run_basedir = "/ifs/data/research/revio/work"
snv_pos_leniency = 1
ins_pos_leniency = 10
ins_size_leniency = 5
del_pos_leniency = 10
del_size_leniency = 5
bnd_pos_leniency = 10
cnv_overlap_minpct = 0.5
str_pos_leniency = 1000
str_size_leniency = 0.5
para_del_dup_inv_size_leniency = 5


def get_and_check_path(run_basedir, target_variant, source):
    if source == "snv":
        basestring = "{}/{}/GR*/SNV*/{}*.vcf.gz"
    elif source == "hificnv":
        basestring = "{}/{}/GR*/CNV*/{}*.vcf"
    elif source == "pbsv":
        basestring = "{}/{}/GR*/SV*/{}*.vcf"
    elif source == "para":
        basestring = "{}/{}/GR*/Par*/{}.general.variants.sorted.vcf"
    elif source == "str":
        basestring = "{}/{}/GR*/STR*/{}.sorted.vcf"
    elif source == "para_json":
        basestring = "{}/{}/GR*/Par*/{}.json"

    vcf_paths = glob.glob(
        basestring.format(
            run_basedir, target_variant["sample"], target_variant["sample"]
        )
    )

    if len(vcf_paths) == 0:
        print("VCF file not found for variant:", target_variant)
        return "File missing"
    if len(vcf_paths) > 1:
        print("Multiple VCF files found:", vcf_paths)
        return "Multiple files found"

    return vcf_paths


def choose_pos_leniency(target_vartype):
    if target_vartype == "SNV":
        return snv_pos_leniency
    elif target_vartype == "INS":
        return ins_pos_leniency
    elif target_vartype == "DEL":
        return del_pos_leniency
    elif target_vartype == "DUP":
        return del_pos_leniency

def get_all_str_motif_permutations(seq):
    results = []
    for i in range(len(seq)):
        shifted = seq[i:] + seq[:i]
        reverse_complement = shifted.translate(str.maketrans("ATCG", "TAGC"))[::-1]
        results.extend([shifted, reverse_complement])
    return results


def reciprocal_overlap(start1, end1, start2, end2):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start >= overlap_end:
        return 0  # No overlap
    overlap_length = overlap_end - overlap_start
    length1 = end1 - start1
    length2 = end2 - start2
    return max(overlap_length / length1, overlap_length / length2)


def condition_same_snp_substitution(target_ref, target_alt, variant_ref, variant_alt):
    if target_ref == variant_ref and target_alt == variant_alt:
        return True
    return False


def condition_start_within_leniency(target_start, variant_start, leniency):
    if abs(target_start - variant_start) <= leniency:
        return True
    return False


def condition_start_end_within_leniency(
    target_start, target_end, variant_start, variant_end, leniency
):
    if (
        abs(target_start - variant_start) <= leniency 
        and abs(target_end - variant_end)<= leniency
    ):
        return True
    return False


def condition_svlen_within_leniency_using_info_svlen(
    target_start, target_end, variant, leniency
):
    if (
        abs(abs(int(variant.INFO.get("SVLEN"))) - (target_end - target_start))
        <= leniency
    ):
        return True
    return False


def condition_inslen_within_leniency(target_inslen, variant, leniency):
    if abs(len(variant.ALT[0]) - target_inslen) <= leniency:
        return True
    return False


def condition_svlen_within_leniency_using_ref_alt(
    target_start, target_end, variant, leniency
):
    if (
        abs(abs(len(variant.REF) - len(variant.ALT[0])) - (target_end - target_start))
        <= leniency
    ):
        return True
    return False


# Variant processing functions (these need to be defined properly based on actual requirements)
def search_snv(
    target_variant,
    run_basedir,
    snv_pos_leniency,
    ins_pos_leniency,
    ins_size_leniency,
    del_pos_leniency,
    del_size_leniency,
):
    vcf_paths = get_and_check_path(run_basedir, target_variant, "snv")
    if (vcf_paths == "File missing") or (vcf_paths == "Multiple files found"):
        return vcf_paths

    vcf_reader = cyvcf2.VCF(vcf_paths[0])

    target_chrom = target_variant["region"].split(":")[0]
    target_start = int(target_variant["region"].split(":")[1].split("-")[0])
    target_vartype = target_variant["vartype"].upper()
    target_interval = "{}:{}-{}".format(
        target_chrom,
        target_start - choose_pos_leniency(target_vartype),
        target_start + choose_pos_leniency(target_vartype),
    )
    if target_vartype == "SNV":
        target_ref, target_alt = target_variant["specific_info"].split(">")[:2]

        for variant in vcf_reader(target_interval):
            if condition_same_snp_substitution(
                target_ref, target_alt, variant.REF, variant.ALT[0]
            ):
                print("Found: {}".format(variant))
                return True

        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype == "INS":
        target_inslen = int(target_variant["specific_info"])

        for variant in vcf_reader(target_interval):
            if variant.var_subtype == "ins" and condition_inslen_within_leniency(
                target_inslen, variant, ins_size_leniency
            ):
                print("Found: {}".format(variant))
                return True
        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype == "DEL|DUP":
        if len(target_variant["region"].split(":")[1].split("-")) == 1:
            target_end = target_start
        else:
            target_end = int(target_variant["region"].split(":")[1].split("-")[1])

        for variant in vcf_reader(target_interval):

            if (
                (variant.var_subtype == "del" or variant.var_subtype == "dup")
                and condition_svlen_within_leniency_using_ref_alt(
                    target_start, target_end, variant, del_size_leniency
                )
            ):
                print("Found: {}".format(variant))
                return True
        print("Missing: {}".format(target_variant))
        return False


def search_hificnv(target_variant, run_basedir, cnv_overlap_minpct):
    # cnv_overlap_minpct = 0.5

    vcf_paths = get_and_check_path(run_basedir, target_variant, "hificnv")
    if (vcf_paths == "File missing") or (vcf_paths == "Multiple files found"):
        return vcf_paths

    vcf_reader = cyvcf2.VCF(vcf_paths[0])

    # get the chromosome and start/end positions from the variant
    target_chrom = target_variant["region"].split(":")[0]
    target_start = int(target_variant["region"].split(":")[1].split("-")[0])
    target_end = int(target_variant["region"].split(":")[1].split("-")[1])
    target_vartype = target_variant["vartype"].upper()

    # Iterate over all variants in the VCF
    for variant in vcf_reader:
        if variant.CHROM == target_chrom and variant.var_subtype == target_vartype:
            rec_overlap = reciprocal_overlap(
                target_start, target_end, variant.start, variant.end
            )
            if rec_overlap >= cnv_overlap_minpct:
                print("Found: {}".format(variant))
                return True
    print("Missing: {}".format(target_variant))
    return False


def search_pbsv(
    target_variant,
    run_basedir,
    snv_pos_leniency,
    ins_pos_leniency,
    ins_size_leniency,
    del_pos_leniency,
    del_size_leniency,
):
    vcf_paths = get_and_check_path(run_basedir, target_variant, "pbsv")
    if (vcf_paths == "File missing") or (vcf_paths == "Multiple files found"):
        return vcf_paths

    vcf_reader = cyvcf2.VCF(vcf_paths[0])

    target_chrom = target_variant["region"].split(":")[0]
    target_vartype = target_variant["vartype"].upper()

    if target_vartype == "INS":
        target_start = int(target_variant["region"].split(":")[1].split("-")[0])
        target_inslen = int(target_variant["specific_info"])

        for variant in vcf_reader:
            if (
                variant.CHROM == target_chrom
                and variant.var_subtype == target_vartype
                and condition_start_within_leniency(
                    target_start, variant.start, ins_pos_leniency
                )
                and condition_inslen_within_leniency(
                    target_inslen, variant, ins_size_leniency
                )
            ):
                print("Found: {}".format(variant))
                return True
        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype in ["DEL", "DUP", "INV"]:
        target_start = int(target_variant["region"].split(":")[1].split("-")[0])
        target_end = int(target_variant["region"].split(":")[1].split("-")[1])
        for variant in vcf_reader:
            if (
                variant.CHROM == target_chrom
                and variant.var_subtype == target_vartype
                and condition_start_end_within_leniency(
                    target_start,
                    target_end,
                    variant.start,
                    variant.end,
                    del_pos_leniency,
                )
                and condition_svlen_within_leniency_using_info_svlen(
                    target_start, target_end, variant, del_size_leniency
                )
            ):
                print("Found: {}".format(variant))
                return True
        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype == "BND":
        target_chrom1 = target_variant["region"].split("-")[0].split(":")[0]
        target_pos1 = int(target_variant["region"].split("-")[0].split(":")[1])
        target_chrom2 = target_variant["region"].split("-")[1].split(":")[0]
        target_pos2 = int(target_variant["region"].split("-")[1].split(":")[1])

        for variant in vcf_reader:
            if (
                variant.var_subtype == "complex"
                and variant.CHROM == target_chrom1
                and condition_start_within_leniency(
                    target_pos1, variant.start, bnd_pos_leniency
                )
                and (
                    variant.ALT[0].replace("]", "[").split(":")[0].split("[")[1].lower()
                    == target_chrom2
                )
                # Hope this here is ok? We check if the 2nd bp is also in the leniency range
                and condition_start_within_leniency(
                    target_pos2,
                    int(variant.ALT[0].replace("]", "[").split(":")[1].split("[")[0]),
                    bnd_pos_leniency,
                )
            ):
                print("Found: {}".format(variant))
                return True

        print("Missing: {}".format(target_variant))
        return False


def search_para(
    target_variant,
    run_basedir,
    snv_pos_leniency,
    ins_pos_leniency,
    del_pos_leniency,
    para_del_dup_inv_size_leniency,
):
    vcf_paths = get_and_check_path(run_basedir, target_variant, "para")
    if (vcf_paths == "File missing") or (vcf_paths == "Multiple files found"):
        return vcf_paths

    vcf_reader = cyvcf2.VCF(vcf_paths[0])

    target_chrom = target_variant["region"].split(":")[0]
    target_start = int(target_variant["region"].split(":")[1].split("-")[0])
    target_vartype = target_variant["vartype"].upper()

    if target_vartype == "SNV":
        target_ref, target_alt = target_variant["specific_info"].split(">")[:2]
        for variant in vcf_reader:
            if (
                variant.CHROM == target_chrom
                and variant.var_type == "snp"
                and condition_start_within_leniency(
                    target_start, variant.start, snv_pos_leniency
                )
                and condition_same_snp_substitution(
                    target_ref, target_alt, variant.REF, variant.ALT[0]
                )
            ):
                print("Found: {}".format(variant))
                return True
        return False

    elif target_vartype == "INS":
        target_inslen = int(target_variant["specific_info"])

        for variant in vcf_reader():
            if (
                variant.CHROM == target_chrom
                and variant.var_subtype == "ins"
                and condition_start_within_leniency(
                    target_start, variant.start, ins_pos_leniency
                )
                and condition_inslen_within_leniency(
                    target_inslen, variant, ins_size_leniency
                )
            ):
                print("Found: {}".format(variant))
                return True
        return False

    elif target_vartype == "DEL":
        if len(target_variant["region"].split(":")[1].split("-")) == 1:
            target_end = target_start
        else:
            target_end = int(target_variant["region"].split(":")[1].split("-")[1])

        for variant in vcf_reader():
            if variant.ALT == []:
                variant.ALT = ["*"]

            if len(variant.ALT[0]) > len(variant.REF):
                var_type = "ins"
            elif len(variant.ALT[0]) < len(variant.REF):
                var_type = "del"
            else:
                var_type = "snp"

            if (
                variant.CHROM == target_chrom
                and var_type == "del"
                and condition_start_end_within_leniency(
                    target_start,
                    target_end,
                    variant.start,
                    variant.end,
                    del_pos_leniency,
                )
                and condition_svlen_within_leniency_using_ref_alt(
                    target_start, target_end, variant, para_del_dup_inv_size_leniency
                )
            ):
                print("Found: {}".format(variant))
                return True

        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype == "INV":
        target_end = int(target_variant["region"].split(":")[1].split("-")[1])
        for variant in vcf_reader():
            if (
                variant.CHROM == target_chrom
                and variant.var_subtype == "INV"
                and condition_start_end_within_leniency(
                    target_start,
                    target_end,
                    variant.start,
                    variant.end,
                    del_pos_leniency,
                )
                and condition_svlen_within_leniency_using_info_svlen(
                    target_start, target_end, variant, para_del_dup_inv_size_leniency
                )
            ):
                print("Found: {}".format(variant))
                return True
        print("Missing: {}".format(target_variant))
        return False


def search_str(target_variant, run_basedir, str_pos_leniency, str_size_leniency):
    vcf_paths = get_and_check_path(run_basedir, target_variant, "str")
    if (vcf_paths == "File missing") or (vcf_paths == "Multiple files found"):
        return vcf_paths

    vcf_reader = cyvcf2.VCF(vcf_paths[0])

    target_chrom = target_variant["region"].split(":")[0]
    target_start = int(target_variant["region"].split(":")[1].split("-")[0])
    target_vartype = target_variant["vartype"].upper()
    target_motif_options = get_all_str_motif_permutations(
        target_variant["specific_info"].split(":")[0]
    )

    for variant in vcf_reader:
        if variant.ALT == []:
            variant.ALT = ["*"]

        target_end = target_start + len(variant.ALT[0])

        effective_len = len(variant.ALT[0]) - len(variant.REF)
        target_sizes = target_variant["specific_info"].split(":")[1]
        target_motif = target_variant["specific_info"].split(":")[0]
        target_len = int(target_sizes.split(">")[1]) - int(target_sizes.split(">")[0])

        frac_ol = min(effective_len, target_len) / max(effective_len, target_len)

        if (
            variant.CHROM == target_chrom
            and variant.var_subtype.upper() == target_vartype
            and condition_start_end_within_leniency(
                target_start, target_end, variant.start, variant.end, str_pos_leniency
            )
            and frac_ol >= str_size_leniency
            and variant.INFO.get("MOTIFS") in target_motif_options
        ):
            print("Found: {}".format(variant))
            return True

    print("MISSING: {}".format(target_variant))
    return False


def search_para_json(target_variant, run_basedir):
    json_genes_underscores = ["smn1", "pms2"]
    json_genes_missing_names = ["rccx"]

    json_paths = get_and_check_path(run_basedir, target_variant, "para_json")
    if (json_paths == "File missing") or (json_paths == "Multiple files found"):
        return json_path

    with open(json_paths[0], "r") as file:
        data = json.load(file)

    target_main_gene = target_variant["region"]
    target_cns_list = target_variant["specific_info"].split(",")
    target_gene_dict = {
        item.split(":")[0]: int(item.split(":")[1]) for item in target_cns_list
    }

    haplotypes = list(data[target_main_gene]["final_haplotypes"].values())

    if target_main_gene == "rccx":
        haplotypes = [element.replace("hap", "rccx_hap") for element in haplotypes]

    if target_main_gene in json_genes_underscores:
        haplotypes = [element.replace("hap", "_hap") for element in haplotypes]

    gene_counts = Counter(haplotype.split("_")[0] for haplotype in haplotypes)

    # Check for zero-value entries in target_gene_dict and add them to gene_counts if not present
    for key, value in target_gene_dict.items():
        if value == 0 and key not in gene_counts:
            gene_counts[key] = 0

    is_subset = all(item in gene_counts.items() for item in target_gene_dict.items())

    if is_subset:
        print("Found: {}".format(target_variant))
        return True

    target_total_n = sum(map(int, target_gene_dict.values()))
    actual_total_n = sum(map(int, gene_counts.values()))

    if target_total_n == actual_total_n:
        print("Found: {}".format(target_variant))
        return "CN correct"

    return False


# Main function to execute
def main(input_variants):
    # Load the variants
    variants = pd.read_csv(
        input_variants,
        sep="\t",
        header=None,
        names=["sample", "source", "vartype", "region", "specific_info"],
    )
    # Results should have all columns from the input plus additional columns for the search results: 'found'
    all_res = variants
    all_res["found"] = "NA"

    for index, variant in variants.iterrows():
        source = variant["source"].lower()

        if source == "snv":
            result = search_snv(
                variant,
                run_basedir,
                snv_pos_leniency,
                ins_pos_leniency,
                ins_size_leniency,
                del_pos_leniency,
                del_size_leniency,
            )
        elif source == "pbsv":
            result = search_pbsv(
                variant,
                run_basedir,
                bnd_pos_leniency,
                ins_pos_leniency,
                ins_size_leniency,
                del_pos_leniency,
                del_size_leniency,
            )
        elif source == "hificnv":
            result = search_hificnv(variant, run_basedir, cnv_overlap_minpct)
        elif source == "para":
            result = search_para(
                variant,
                run_basedir,
                snv_pos_leniency,
                ins_pos_leniency,
                del_pos_leniency,
                para_del_dup_inv_size_leniency,
            )
        elif source == "para_json":
            result = search_para_json(variant, run_basedir)
        elif source == "str":
            result = search_str(
                variant, run_basedir, str_pos_leniency, str_size_leniency
            )

        else:
            print(
                f"Error: Source type '{source}' is not supported. Supported are: snv, pbsv, hificnv, para, para_json, str"
            )
            continue

        # Update the 'found' column
        all_res.at[index, "found"] = result

    print(all_res)

    # save it to a file
    output_file = 'out.txt'
    all_res.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_variants_file>")
        sys.exit(1)
    main(sys.argv[1])
