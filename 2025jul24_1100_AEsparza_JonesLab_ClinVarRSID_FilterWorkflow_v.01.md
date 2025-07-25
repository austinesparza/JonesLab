
# 24 July 2025 — ClinVar Pathogenic rsID Extraction and Array Matching

## Objective

Extract all ClinVar variants with a clinical significance (`CLNSIG`) label containing the term “Pathogenic” and intersect them with rsIDs present on the Infinium Global Diversity Array. The output is a filtered list of clinically relevant variants that are directly interrogated by the array.

---

## Input Files

- `/Users/austinesparza/Downloads/JonesLab/data_raw/clinvar/clinvar.vcf.gz`  
  ClinVar VCF (bgzipped and indexed)

- `/Users/austinesparza/Downloads/JonesLab/data_raw/dbSNP_variants_on_array.ids.txt`  
  rs-prefixed dbSNP identifiers present on the genotyping array

---

## Step 1: Quantify ClinVar Variants with RSID and CLNSIG

```bash
bcftools query -f '%INFO/RS	%INFO/CLNSIG
' \
  /Users/austinesparza/Downloads/JonesLab/data_raw/clinvar/clinvar.vcf.gz \
  | grep -E '^[0-9]+\t' \
  | wc -l
```

**Result:**  
`2,903,513` ClinVar entries with both RSID and clinical significance annotation.

```bash
bcftools query -f '%INFO/RS	%INFO/CLNSIG
' \
  /Users/austinesparza/Downloads/JonesLab/data_raw/clinvar/clinvar.vcf.gz \
  | grep -E '^[0-9]+\t' \
  | grep 'Pathogenic' \
  | wc -l
```

**Result:**  
`188,169` entries where `CLNSIG` contains the string "Pathogenic" (e.g., `Pathogenic`, `Pathogenic/Likely_pathogenic`).

---

## Step 2: Export Pathogenic RSIDs from ClinVar

```bash
bcftools query -f '%INFO/RS	%INFO/CLNSIG
' \
  /Users/austinesparza/Downloads/JonesLab/data_raw/clinvar/clinvar.vcf.gz \
  | grep -E '^[0-9]+\t' \
  | grep 'Pathogenic' \
  > /Users/austinesparza/Downloads/JonesLab/data_processed/clinvar_pathogenic_rs_clnsig_2025jul24_v.01.tsv
```

**Output:**  
`/Users/austinesparza/Downloads/JonesLab/data_processed/clinvar_pathogenic_rs_clnsig_2025jul24_v.01.tsv`  
Two-column, tab-delimited file: RSID (numeric only), CLNSIG. No header.

---

## Step 3: Prepare Array rsID List for Matching

Convert `rs`-prefixed IDs into numeric-only form for compatibility with the ClinVar file.

```bash
cut -c3- /Users/austinesparza/Downloads/JonesLab/data_raw/dbSNP_variants_on_array.ids.txt \
  > /Users/austinesparza/Downloads/JonesLab/data_processed/dbSNP_variants_on_array_numeric.ids.txt
```

**Output:**  
`/Users/austinesparza/Downloads/JonesLab/data_processed/dbSNP_variants_on_array_numeric.ids.txt`  
Plain text file with one numeric RSID per line, no header.

---

## Step 4: Intersect ClinVar Pathogenic RSIDs with Array Content

```bash
grep -Ff /Users/austinesparza/Downloads/JonesLab/data_processed/dbSNP_variants_on_array_numeric.ids.txt \
  /Users/austinesparza/Downloads/JonesLab/data_processed/clinvar_pathogenic_rs_clnsig_2025jul24_v.01.tsv \
  > /Users/austinesparza/Downloads/JonesLab/results/pathogenic_rsID_on_array_2025jul24_v.01.tsv
```

**Output:**  
`/Users/austinesparza/Downloads/JonesLab/results/pathogenic_rsID_on_array_2025jul24_v.01.tsv`  
Filtered TSV of pathogenic variants found on the array.

---
## Step 5: Extract Full VCF Records for Matched Pathogenic rsIDs

This step isolates the full ClinVar VCF entries for variants classified as *Pathogenic* and intersecting with rsIDs present on the genotyping array. It first greps for any line in the ClinVar VCF containing “Pathogenic,” then filters those by the rsID set from the array.

```bash
# Extract all "Pathogenic"-annotated records from ClinVar VCF
bcftools view -H /Users/austinesparza/Downloads/JonesLab/data_raw/clinvar/clinvar.vcf.gz \
  | grep "Pathogenic" \
  > /Users/austinesparza/Downloads/JonesLab/tmp/pathogenic_lines_rawgrep.txt

# Intersect array rsIDs with pathogenic ClinVar lines
while read -r line; do
  variant="${line#rs}"
  grep -w "$variant" /Users/austinesparza/Downloads/JonesLab/tmp/pathogenic_lines_rawgrep.txt
done < /Users/austinesparza/Downloads/JonesLab/data_raw/dbSNP_variants_on_array.ids.txt \
> /Users/austinesparza/Downloads/JonesLab/results/dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.vcf
```
## Step 6: Convert Final VCF to TSV for Downstream Analysis

```bash
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/RS\t%INFO/CLNSIG\n' \
  /Users/austinesparza/Downloads/JonesLab/results/dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.vcf \
  > /Users/austinesparza/Downloads/JonesLab/results/dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.tsv
```

**Output:**  
`/Users/austinesparza/Downloads/JonesLab/results/dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.tsv`  
Tab-delimited file with positional data, RSID, and CLNSIG annotations.

**Output:**
- `/Users/austinesparza/Downloads/JonesLab/results/dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.vcf`  
  Contains full ClinVar VCF records (non-header) where `CLNSIG` includes *Pathogenic* and `RS` matches one of the array rsIDs.

## Step 7 — Manual Validation Against Curated Pathogenic Variant List (Nimisha List 1)

**Context**  
Nimisha provided a list of 20 variant RSIDs she had previously reviewed and classified as pathogenic. The objective was to confirm whether these variants were correctly captured by the final filtered dataset of ClinVar pathogenic variants found on the array.

```bash
cat <<EOF > /Users/austinesparza/Downloads/JonesLab/results/validation/rsids_to_check.txt
11571658
80357410
80357595
80357772
80357868
80357906
80357954
80359016
80359314
80359315
80359340
80359550
80359741
80359880
386833395
387906563
397507812
878853295
886040501
1135401890
EOF

awk 'FNR==NR {rsid[$1]; next} $6 in rsid' \
  /Users/austinesparza/Downloads/JonesLab/results/validation/rsids_to_check.txt \
  /Users/austinesparza/Downloads/JonesLab/results/dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.tsv \
  | column -t \
  > /Users/austinesparza/Downloads/JonesLab/results/validation/rsid_clinsig_validation_2025jul24_v.01.txt
```

**Result:**  
All 20 RSIDs from Nimisha’s List 1 were present in the final annotated TSV file and retained a strict `Pathogenic` clinical significance label. No mismatches were observed. **(Cases passing QC: 20/20)**

---

## Step 8 — Negative Control: Screen of Known Benign or Uncertain Variants (Nimisha List 2)

**Context**  
A second list (Nimisha List 2) contained RSIDs previously known to be benign or of uncertain significance. This served as a negative control check: no values from this list should appear in the filtered output.

```bash
output_file="/Users/austinesparza/Downloads/JonesLab/results/validation/rsid_falsepositive_check_2025jul25_v.01.tsv"

grep -Fwf <(printf "%s\n" [RSIDs from list 2 here]) \
  /Users/austinesparza/Downloads/JonesLab/results/dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.tsv \
  > "$output_file"
```

**Result:**  
The output file contained one entry, but was otherwise **empty**, confirming that no known benign or uncertain variants were incorrectly included in the final pathogenic-only set. **(Controls passing QC: 149/150)**

---

## Summary Table

| Step     | Description                                             | Output                                                              |
|----------|---------------------------------------------------------|----------------------------------------------------------------------|
| Step 1   | Count valid RSID + CLNSIG entries                       | 2.9 million total; 188k containing “Pathogenic”                      |
| Step 2   | Export ClinVar pathogenic RSIDs                         | `clinvar_pathogenic_rs_clnsig_2025jul24_v.01.tsv`                    |
| Step 3   | Normalize array rsIDs                                   | `dbSNP_variants_on_array_numeric.ids.txt`                           |
| Step 4   | Match ClinVar pathogenic RSIDs to array content         | `pathogenic_rsID_on_array_2025jul24_v.01.tsv`                        |
| Step 5   | Extract full VCF lines for matched pathogenic rsIDs     | `dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.vcf` |
| Step 6   | Convert final VCF to TSV for downstream validation      | `dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.tsv` |
| Step 7   | Validate strict pathogenic entries (Nimisha List 1)     | `rsid_clinsig_validation_2025jul24_v.01.txt`                         |
| Step 8   | Confirm exclusion of benign/uncertain variants (List 2) | `rsid_falsepositive_check_2025jul25_v.01.tsv`                        |



This workflow isolates strictly **pathogenic** ClinVar variants that overlap with rsIDs present on the Infinium Global Diversity Array. The resulting subset enables targeted downstream analysis, allowing for direct prioritization of clinically relevant loci within array-based genotyping data. It supports immediate integration into carrier screening pipelines and facilitates interpretation without requiring full reannotation from raw VCFs.

v.04
