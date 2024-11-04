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
from typing import List, Tuple


snv_pos_leniency = 1

indel_pos_leniency = 10
indel_size_leniency = 5

sv_overlap_minpct = 0.5
cnv_overlap_minpct = 0.5

bnd_pos_leniency = 50

str_pos_leniency = 1000
str_cn_fraction_leniency = 20 
mt_percentpoints_leniency = 20


def get_and_check_path_downsample(
    run_basedir, permutation_index, target_variant, source
):
    if source == "snv":
        basestring = "{}/{}/i{}/SNV*/{}*.vcf.gz"
    elif source == "hificnv":
        basestring = "{}/{}/i{}/CNV*/{}*.vcf"
    elif source == "pbsv":
        basestring = "{}/{}/i{}/SV*/{}*.vcf"
    elif source == "para":
        basestring = "{}/{}/i{}/Par*/{}.general.variants.sorted.vcf"
    elif source == "str":
        basestring = "{}/{}/i{}/STR*/{}.sorted.vcf*"
    elif source == "para_json":
        basestring = "{}/{}/i{}/Par*/{}.json"
    elif source == "mt":
        basestring = "{}/{}/i{}/MT*/{}*hcdiffs.txt"

    vcf_paths = glob.glob(
        basestring.format(
            run_basedir,
            target_variant["sample"],
            permutation_index,
            target_variant["sample"],
        )
    )

    # Filter out any index files (such as .csi)
    vcf_paths = [
        path
        for path in vcf_paths
        if path.endswith((".vcf", ".vcf.gz", ".json", "hcdiffs.txt"))
    ]

    if len(vcf_paths) == 0:
        print("VCF file not found for variant:", target_variant)
        return "File-missing"
    if len(vcf_paths) > 1:
        print("Multiple VCF files found:", vcf_paths)
        return "Multiple-files-found"

    return vcf_paths


def get_and_check_path(run_basedir, target_variant, source):
    sourcedir = 'full' #perm
    
    if sourcedir == 'full':
        if source == "snv":
            basestring = "{}/{}/G*/SNV*/{}*.vcf.gz"
        elif source == "hificnv":
            basestring = "{}/{}/G*/CNV*/{}*.vcf"
        elif source == "pbsv":
            basestring = "{}/{}/G*/SV*/{}*.vcf"
        elif source == "para":
            basestring = "{}/{}/G*/Par*/{}.general.variants.sorted.vcf"
        elif source == "str":
            basestring = "{}/{}/G*/STR*/{}.sorted.vcf*"
        elif source == "para_json":
            basestring = "{}/{}/G*/Par*/{}.json"
        elif source == "mt":
            basestring = "{}/{}/N*/MT*/{}*hcdiffs.txt"
    elif sourcedir=='perm':
        if source == "snv":
            basestring = "{}/{}/SNV*/{}*.vcf.gz"
        elif source == "hificnv":
            basestring = "{}/{}/CNV*/{}*.vcf"
        elif source == "pbsv":
            basestring = "{}/{}/SV*/{}*.vcf"
        elif source == "para":
            basestring = "{}/{}/Par*/{}.general.variants.sorted.vcf"
        elif source == "str":
            basestring = "{}/{}/STR*/{}.sorted.vcf*"
        elif source == "para_json":
            basestring = "{}/{}/Par*/{}.json"
        elif source == "mt":
            basestring = "{}/{}/MT*/{}*hcdiffs.txt"

    vcf_paths = glob.glob(
        basestring.format(
            run_basedir, target_variant["sample"], target_variant["sample"]
        )
    )

    # Filter out any index files (such as .csi)
    vcf_paths = [
        path
        for path in vcf_paths
        if path.endswith((".vcf", ".vcf.gz", ".json", "hcdiffs.txt"))
    ]

    if len(vcf_paths) == 0:
        print("VCF file not found for variant:", target_variant)
        return "File-missing"
    if len(vcf_paths) > 1:
        print("Multiple VCF files found:", vcf_paths)
        return "Multiple-files-found"

    return vcf_paths


def choose_pos_leniency(target_vartype):
    if target_vartype in ["SNV", "Substitution"]:
        return snv_pos_leniency
    elif target_vartype in ["INS", "DEL", "DUP"]:
        return indel_pos_leniency


    # check if target_vartype is ins, del or dup

def get_all_str_motif_permutations(seq):
    results = []
    for i in range(len(seq)):
        shifted = seq[i:] + seq[:i]
        reverse_complement = shifted.translate(str.maketrans("ATCG", "TAGC"))[::-1]
        results.extend([shifted, reverse_complement])
    return results

# Function to extract motif counts for a predefined motif
def extract_specific_motif_count(variant, motif_of_interest: str) -> List[int]:
    
    motifs = variant.INFO.get("MOTIFS").split(',')
    motif_counts = variant.format('MC')[0]
    if (motif_counts == '.'):
        return []
    if motif_of_interest in motifs:
        motif_index = motifs.index(motif_of_interest)
        mc_fields = motif_counts.split(',')
        counts = [int(field.split('_')[motif_index]) for field in mc_fields]
        return counts
    return []

# Function to check if counts are within leniency range
def check_counts_within_leniency(counts: List[int], desired: Tuple[int, int], leniency: float) -> bool:
    

    condition_minimum_1 = desired[0] < 0
    condition_minimum_2 = desired[1] < 0
    
    counts = [abs(count) for count in counts]

    if len(counts) != 2:
        return False

    # Calculate the acceptable range for each desired count
    desired1_min = desired[0] * (1 - leniency / 100)
    desired1_max = desired[0] * (1 + leniency / 100)
    desired2_min = desired[1] * (1 - leniency / 100)
    desired2_max = desired[1] * (1 + leniency / 100)

    # Check if both counts are within the acceptable range
    if condition_minimum_1:
        count1_within_range = desired1_min <= counts[0]
    else:
        count1_within_range = desired1_min <= counts[0] <= desired1_max

    if condition_minimum_2:
        count2_within_range = desired2_min <= counts[1]
    else:
        count2_within_range = desired2_min <= counts[1] <= desired2_max
    
    return count1_within_range and count2_within_range


def reciprocal_overlap(start1, end1, start2, end2):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start >= overlap_end:
        return 0  # No overlap
    overlap_length = overlap_end - overlap_start
    length1 = end1 - start1
    length2 = end2 - start2
    return min(overlap_length / length1, overlap_length / length2)

def merge_cnv_intervals(intervals, distance):
    if not intervals:
        return []

    # Sort intervals by start position
    intervals.sort()

    merged_intervals = [intervals[0]]

    for current_start, current_end in intervals[1:]:
        last_start, last_end = merged_intervals[-1]

        # Check if the current interval overlaps or is within the distance threshold
        if current_start <= last_end + distance:
            merged_intervals[-1] = (last_start, max(last_end, current_end))
        else:
            merged_intervals.append((current_start, current_end))

    return merged_intervals

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
        and abs(target_end - variant_end) <= leniency
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


def condition_inslen_within_leniency_using_ref_alt(target_inslen, variant, leniency):
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


def search_mt(target_variant, vcf_paths, snv_pos_leniency, mt_percentpoints_leniency):

    tsv = pd.read_csv(vcf_paths[0], sep="\t", header=0)
    target_chrom = target_variant["region"].split(":")[0]
    target_start = int(target_variant["region"].split(":")[1].split("-")[0])
    target_vartype = target_variant["vartype"].upper()
    if target_vartype.upper() == "SNV":
        target_vartype == "Substitution"

    target_ref, target_alt = (
        target_variant["specific_info"].split(":")[0].split(">")[:2]
    )
    target_percent = float(target_variant["specific_info"].split(":")[1])
    for index, row in tsv.iterrows():
        if (
            row["Chromosome"] == target_chrom
            and condition_start_within_leniency(
                target_start, row["Start position"], snv_pos_leniency
            )
            and abs(row["% variation"] - target_percent) <= mt_percentpoints_leniency
            and condition_same_snp_substitution(
                target_ref, target_alt, row["Reference"], row["Variant"]
            )
        ):
            print("Found: {}".format(row))
            return True

    return False


def search_snv(
    target_variant,
    vcf_paths,
    snv_pos_leniency,
    indel_pos_leniency,
    indel_size_leniency,
    sv_overlap_minpct
):
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

            vartype_condition = variant.var_subtype.upper() == target_vartype

            if target_inslen <= 50:
                position_condition = condition_inslen_within_leniency_using_ref_alt(target_inslen, variant, indel_size_leniency)
            else:
                position_condition = reciprocal_overlap(target_start, target_start+target_inslen, variant.start, variant.start+len(variant.ALT[0])) >= sv_overlap_minpct
                
            if (vartype_condition and position_condition):
                print("Found: {}".format(variant))
                return True
            
        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype in ["DEL", "DUP"]:
        if len(target_variant["region"].split(":")[1].split("-")) == 1:
            target_end = target_start
        else:
            target_end = int(target_variant["region"].split(":")[1].split("-")[1])

        for variant in vcf_reader(target_interval):
            
            vartype_condition = variant.var_subtype.upper() == target_vartype
            
            if target_end - target_start <= 50:
                position_condition = condition_svlen_within_leniency_using_ref_alt(target_start, target_end, variant, indel_size_leniency)
            else:
                position_condition = reciprocal_overlap(target_start, target_end, variant.start, variant.end) >= sv_overlap_minpct

            if (vartype_condition and position_condition):
                print("Found: {}".format(variant))
                return True
                
        print("Missing: {}".format(target_variant))
        return False

def search_hificnv(target_variant, vcf_paths, cnv_overlap_minpct):
    target_chrom, target_coords = target_variant["region"].split(":")
    target_start, target_end = map(int, target_coords.split("-"))

    vcf_reader = cyvcf2.VCF(vcf_paths[0])
    intervals = [
        (variant.start, variant.end)
        for variant in vcf_reader
        if variant.CHROM == target_chrom
    ]
    merged_intervals = merge_cnv_intervals(intervals, 1000)

    for start, end in merged_intervals:
        if reciprocal_overlap(target_start, target_end, start, end) >= cnv_overlap_minpct:
            print(f"Found: {target_chrom}:{start}-{end}")
            return True

    print(f"Missing: {target_variant}")
    return False




def search_pbsv(
    target_variant,
    vcf_paths,
    snv_pos_leniency,
    indel_pos_leniency,
    indel_size_leniency,
    sv_overlap_minpct
):
    vcf_reader = cyvcf2.VCF(vcf_paths[0])

    target_chrom = target_variant["region"].split(":")[0]
    target_vartype = target_variant["vartype"].upper()

    if target_vartype == "INS":
        target_start = int(target_variant["region"].split(":")[1].split("-")[0])
        target_inslen = int(target_variant["specific_info"])

        for variant in vcf_reader:
            
            vartype_condition = variant.var_subtype.upper() == target_vartype
            
            if target_inslen <= 50:
                startpos_condition = condition_start_within_leniency(target_start, variant.start, indel_pos_leniency)
            else:
                startpos_condition = condition_start_within_leniency(target_start, variant.start, target_inslen)
                
            if not (variant.CHROM == target_chrom and vartype_condition and startpos_condition):
                continue
            
            if target_inslen <= 50:
                svlen_condition = condition_inslen_within_leniency_using_ref_alt(target_inslen, variant, indel_size_leniency)
            else:
                svlen_condition = reciprocal_overlap(target_start, target_start+target_inslen, variant.start, variant.start+len(variant.ALT[0])) >= sv_overlap_minpct
            
            if (svlen_condition):
                print("Found: {}".format(variant))
                return True
            
        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype == "DEL":
        target_start = int(target_variant["region"].split(":")[1].split("-")[0])
        target_end = int(target_variant["region"].split(":")[1].split("-")[1])
        
        for variant in vcf_reader:
            
            vartype_condition = variant.var_subtype.upper() == target_vartype

            if target_end-target_start <= 50:
                start_end_pos_condition = condition_start_end_within_leniency(target_start,target_end,variant.start,variant.end,indel_pos_leniency)
            else:
                start_end_pos_condition = condition_start_end_within_leniency(target_start,target_end,variant.start,variant.end,target_end-target_start)
                                  
            if not (variant.CHROM == target_chrom and vartype_condition and start_end_pos_condition):
                continue
            
            if target_end-target_start <= 50:
                svlen_condition = condition_svlen_within_leniency_using_ref_alt(target_start, target_end, variant, indel_size_leniency)
            else:
                svlen_condition = reciprocal_overlap(target_start, target_end, variant.start, variant.end) >= sv_overlap_minpct
                
            if (svlen_condition):
                print("Found: {}".format(variant))
                return True    
            
        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype in ["DUP", "INV"]:
        target_start = int(target_variant["region"].split(":")[1].split("-")[0])
        target_end = int(target_variant["region"].split(":")[1].split("-")[1])
        
        for variant in vcf_reader:
            
            vartype_condition = variant.var_subtype.upper() == target_vartype

            if target_end-target_start <= 50:
                start_end_pos_condition = condition_start_end_within_leniency(target_start,target_end,variant.start,variant.end,indel_pos_leniency)
            else:
                start_end_pos_condition = condition_start_end_within_leniency(target_start,target_end,variant.start,variant.end,target_end-target_start)
                    
            if not (variant.CHROM == target_chrom and vartype_condition and start_end_pos_condition):
                continue
            
            if target_end-target_start <= 50:
                svlen_condition = condition_svlen_within_leniency_using_info_svlen(target_start, target_end, variant, indel_size_leniency)
            else:
                svlen_condition = reciprocal_overlap(target_start, target_end, variant.start, variant.end) >= sv_overlap_minpct
                
            if (svlen_condition):
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
                    == target_chrom2.lower()
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
    vcf_paths,
    snv_pos_leniency,
    indel_pos_leniency,
    indel_size_leniency
):
    vcf_reader = cyvcf2.VCF(vcf_paths[0])

    target_chrom = target_variant["region"].split(":")[0]
    target_start = int(target_variant["region"].split(":")[1].split("-")[0])
    target_vartype = target_variant["vartype"].upper()
    if target_vartype == "SNV":
        target_ref, target_alt = target_variant["specific_info"].split(">")[:2]
        for variant in vcf_reader:
            if (
                variant.CHROM == target_chrom
                and variant.var_type.upper() in ["SNP", "SNV"]
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
                and variant.var_type.upper() == target_vartype
                and condition_start_within_leniency(
                    target_start, variant.start, indel_pos_leniency
                )
                and condition_inslen_within_leniency_using_ref_alt(
                    target_inslen, variant, indel_size_leniency
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

            if target_end - target_start <= 50:
                position_condition = condition_start_end_within_leniency(target_start,target_end,variant.start,variant.end,indel_pos_leniency)
            else:
                position_condition = condition_start_end_within_leniency(target_start,target_end,variant.start,variant.end,target_end - target_start)
            
            if not (variant.CHROM == target_chrom and var_type == "del" and position_condition):
                continue
            
            if target_end - target_start <= 50:
                svlen_condition = condition_svlen_within_leniency_using_ref_alt(target_start, target_end, variant, indel_size_leniency)
            else:
                svlen_condition = reciprocal_overlap(target_start, target_end, variant.start, variant.end) >= sv_overlap_minpct
            
            if (svlen_condition):
                print("Found: {}".format(variant))
                return True

        print("Missing: {}".format(target_variant))
        return False

    elif target_vartype == "INV":
        target_end = int(target_variant["region"].split(":")[1].split("-")[1])
        for variant in vcf_reader():
            
            if (variant.var_subtype.upper() != "INV"):
                continue

            if target_end - target_start <= 50:
                position_condition = condition_start_end_within_leniency(target_start,target_end,variant.start,variant.start + int(variant.INFO.get("SVLEN")),indel_pos_leniency)
            else:
                position_condition = condition_start_end_within_leniency(target_start,target_end,variant.start,variant.start + int(variant.INFO.get("SVLEN")),target_end - target_start)
                
        
            if not (variant.CHROM == target_chrom and variant.var_subtype.upper() == "INV" and position_condition):
                continue
                
            if target_end - target_start <= 50:
                svlen_condition = condition_svlen_within_leniency_using_info_svlen(target_start, target_end, variant, indel_size_leniency)
            else:
                svlen_condition = reciprocal_overlap(target_start, target_end, variant.start, variant.start + int(variant.INFO.get("SVLEN"))) >= sv_overlap_minpct
            
            if (svlen_condition):
                print("Found: {}".format(variant))
                return True
            
        print("Missing: {}".format(target_variant))
        return False


def search_str(target_variant, vcf_paths, str_pos_leniency, str_cn_fraction_leniency):
    vcf_reader = cyvcf2.VCF(vcf_paths[0])

    target_chrom = target_variant["region"].split(":")[0]
    target_start = int(target_variant["region"].split(":")[1].split("-")[0])
    
    target_sizes = target_variant["specific_info"].split(":")[1]
    
    # Little magic to handle the 'plus' ('p') in sizes
    target_sizes = ','.join(['-' + x[:-1] if 'p' in x else x for x in target_sizes.split(',')])
    
    target_motif = target_variant["specific_info"].split(":")[0]
    
    target_motifcounts = [int(target_sizes.split(",")[0]), int(target_sizes.split(",")[1])]
    for variant in vcf_reader:

        variant_motifcounts = extract_specific_motif_count(variant, target_motif)
        if len(variant_motifcounts) == 0:
            continue
        
        if (
            variant.CHROM == target_chrom
            and condition_start_within_leniency(
                target_start, variant.start, str_pos_leniency
                )
            and check_counts_within_leniency(variant_motifcounts, target_motifcounts, str_cn_fraction_leniency)
            ):
                print("Found: {}".format(variant))
                return True

    print("MISSING: {}".format(target_variant))
    return False

def search_para_json(target_variant, json_paths):
    with open(json_paths[0], "r") as file:
        data = json.load(file)

    json_genes_underscores = ["smn1", "pms2"]
    json_genes_missing_names = ["rccx"]

    target_main_gene = target_variant["region"]
   
    # Sometimes the json has no values at all. Not sure why.
    # See labbook entry [para_json_null] from 8th May 2024
    if all(value is None for value in data[target_main_gene].values()):
        return False 
    
    # If there is a '+' in target_variant["specific_info"], split by it
    if (target_variant["vartype"] == 'opn1lw_as'):
        target_gene,target_asconfig = target_variant["specific_info"].split(':')
        if (target_asconfig == 'LIVVA'):
            target_string = target_gene
        else:
            target_string = "{}_{}".format(target_gene, target_asconfig)
            
        if (target_string in data[target_main_gene]['annotated_haplotypes'].values()):
            return True
        else:
            return False
            
    if "+" in target_variant["specific_info"]:
        target_cns_list = target_variant["specific_info"].split("+")
        return_combined_cns = True
    else:
        target_cns_list = target_variant["specific_info"].split(",")
        return_combined_cns = False
        
    target_gene_dict = {
        item.split(":")[0]: int(item.split(":")[1]) for item in target_cns_list
    }


    haplotypes = list(data[target_main_gene]["final_haplotypes"].values())

    if target_main_gene in json_genes_missing_names:
        haplotypes = [
            element.replace("hap", "{}_hap".format(target_main_gene))
            for element in haplotypes
        ]

    if target_main_gene in json_genes_underscores:
        haplotypes = [element.replace("hap", "_hap") for element in haplotypes]

    gene_counts = Counter(haplotype.split("_")[0] for haplotype in haplotypes)

    # for smn1 we have separate rules because the json is inconsistent as hell.
    if target_main_gene == "smn1":
        gene_counts = {
            "smn1": data[target_main_gene]["smn1_cn"],
            "smn2": data[target_main_gene]["smn2_cn"],
        }

    # where the values of gene_counts are 'null', convert to 0
    gene_counts = {k: 0 if v is None else v for k, v in gene_counts.items()}
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
        if (return_combined_cns):
            return True
        else:
            return "CN-correct"

    return False


# Main function to execute
def main(input_variants, run_basedir, output_file, permutation_index=None):
    # Load the variants
    variants = pd.read_csv(
        input_variants,
        sep="\t",
        header=None,
        names=["sample", "variantID", "source", "vartype", "region", "specific_info"],
    )
    # Results should have all columns from the input plus additional columns for the search results: 'found'
    all_res = variants
    all_res["found"] = "NA"
    for index, variant in variants.iterrows():
        source = variant["source"].lower()
        
        if permutation_index:
            infile_paths = get_and_check_path_downsample(
                run_basedir, permutation_index, variant, source
            )
        else:
            infile_paths = get_and_check_path(run_basedir, variant, source)

        if (infile_paths == "File-missing") or (infile_paths == "Multiple-files-found"):
            # Update the 'found' column
            all_res.at[index, "found"] = infile_paths
            continue

        if source == "snv":
            result = search_snv(
                variant,
                infile_paths,
                snv_pos_leniency,
                indel_pos_leniency,
                indel_size_leniency,
                sv_overlap_minpct,
            )
        elif source == "pbsv":
            result = search_pbsv(
                variant,
                infile_paths,
                bnd_pos_leniency,
                indel_pos_leniency,
                indel_size_leniency,
                sv_overlap_minpct
            )
        elif source == "hificnv":
            result = search_hificnv(variant, infile_paths, cnv_overlap_minpct)
        elif source == "para":
            result = search_para(
                variant,
                infile_paths,
                snv_pos_leniency,
                indel_pos_leniency,
                indel_size_leniency
            )
        elif source == "para_json":
            result = search_para_json(variant, infile_paths)
        elif source == "str":
            result = search_str(
                variant, infile_paths, str_pos_leniency, str_cn_fraction_leniency
            )
        elif source == "mt":
            result = search_mt(
                variant, infile_paths, snv_pos_leniency, mt_percentpoints_leniency
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
    all_res.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    if not len(sys.argv) in {4, 5}:
        print(
            "Usage: python find_lrs100_variants.py <input_variants_file> <run_basedir> <output_file> <permutation_index (optional)>"
        )
        sys.exit(1)

    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        sys.argv[4] if len(sys.argv) == 5 else None,
    )


