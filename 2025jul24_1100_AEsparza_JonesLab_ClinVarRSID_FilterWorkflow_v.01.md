
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



## Summary

| Step     | Description                                             | Output                                                              |
|----------|---------------------------------------------------------|----------------------------------------------------------------------|
| Step 1   | Count valid RSID + CLNSIG entries                       | 2.9 million total; 188k containing “Pathogenic”                     |
| Step 2   | Export ClinVar pathogenic RSIDs                         | `clinvar_pathogenic_rs_clnsig_2025jul24_v.01.tsv`                   |
| Step 3   | Normalize array rsIDs                                   | `dbSNP_variants_on_array_numeric.ids.txt`                           |
| Step 4   | Match ClinVar pathogenic RSIDs to array content         | `pathogenic_rsID_on_array_2025jul24_v.01.tsv`                        |
| Step 5   | Extract full VCF entries for matched pathogenic rsIDs   | `dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.vcf` |
| Step 6   | Convert filtered VCF entries to TSV format              | `dbSNP_variants_on_array_pathogenic_clinvar_grepmatch_2025jul24_v.01.tsv` |


This workflow enables efficient identification of ClinVar-classified pathogenic variants that are directly assayed by the Infinium Global Diversity Array. By intersecting array rsIDs with ClinVar's annotated pathogenic entries, it provides a rapid mechanism to flag clinically significant markers already represented in the array design, supporting downstream variant prioritization, carrier identification, and assessment of clinical relevance from array-based genotyping data.

v.03
