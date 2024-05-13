
# Automated Re-Discovery of Disease-Associated Variants

## Introduction
This Python script checks for the presence or absence of disease-associated variants in various file formats (VCF, TXT, JSON) based on a given list of samples and their variants. Developed for a downsampling experiment, this tool estimates how many of the Revio100 'disease' variants can be re-discovered at various coverages (10X, 15X, 20X).

## Dependencies
This script requires the following dependencies:

- `Python >= 3.7.1`
- `cyvcf2`
- `pandas`


## Input Format

The script expects a TSV (Tab Separated Values) file with specific columns. Below are examples for each supported SV class

| Sample | Source | SVType | Coordinates | Specific Info |
|--------|--------|--------|-------------|---------------|
| P1-A1  | snv    | ins    | chrX:25013653 | 33 |
| P1-A1  | snv    | snv    | chr16:21757209 | C>T |
| P1-A1  | snv    | ins    | chr2:233760233 | 2 |
| P1-A1  | snv    | del    | chr16:21715122 | none |
| P1-A1  | pbsv   | bnd    | chr11:121205665-chr7:100147315 | none |
| P1-A1  | pbsv   | ins    | chrX:147912048 | 138 |
| P1-A1  | pbsv   | del    | chr15:67063774-67067017 | none |
| P1-A1  | pbsv   | inv    | chr6:118762532-118815256 | none |
| P1-A1  | hificnv| dup    | chr11:116268000-121206000 | none |
| P1-A1  | hificnv| del    | chr17:45618001-46134000 | none |
| P1-A1  | para   | snv    | chr6:32038297 | C>T |
| P1-A1  | para_json | cnv  | opn1lw | opn1lw:0 |
| P1-A1  | para_json | cnv  | opn1lw | opn1lw:1,opn1mw:2 |
| P1-A1  | str    | ins    | chrX:25013649:25013698 | GCC:49>82 |
| P1-A1  | str    | ins    | chr3:129172576-129172732 | AGGC:156>12726 |

**Supported Sources and SVTypes:**
- `snv` [snv, ins, del]
- `pbsv` [bnd, ins, del, inv]
- `hificnv` [dup, del]
- `para` [ins]
- `para_json` [cnv]
- `str` [ins]

**Contents of Specific.info Field (default: none):**
- Source 'snv/pbsv/para', SVType 'snv': Substitution type (e.g., C>T)
- Source 'snv/pbsv/para', SVType 'ins': Insertion length (bp)
- Source 'str', SVType 'ins': Motif: original length (bp) > new length (bp)
- Source 'para_json', SVtype 'cnv': expected copy-numbers of defined pseudogenes

**Special case para_json**

- Coordinates: name of the non-pseudo version of the gene whose cnv we want to check. 
- Specific.Info: comma-separated list of expected copy-numbers for genes and pseudogenes that are listed under the 'main' gene in the json file. 

## Usage
To run the tool, navigate to the script directory and execute the following command:

```bash
python find_lrs100_variants.py <input_variants_file> /path/to/your/samplefolders/ <output_file>
```

**Note:** If running on Turbo, ensure the `get_and_check_path` function is correctly pointing to the sample folders, as the folder structure may differ based on the downsampling runs.

## Output
The script produces an output file reporting the genotyping status. The additional column takes the form: 

- TRUE: Variant was re-discovered
- FALSE: Variant is missing
- FILE-MISSING: No variant file was found to search for the variant (e.g., false vcf paths)
- CN-correct: For para_json: the gene/pseudogene CN were wrong, but the sum of CN was correct. This means that a gene del/dup was spotted, but wrongly assignned to another copy of the gene/pseudogene. These cases require careful investigation. 

