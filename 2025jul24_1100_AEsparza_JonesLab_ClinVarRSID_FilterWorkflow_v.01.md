
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

## Summary

| Step     | Description                                     | Output                                                   |
|----------|-------------------------------------------------|----------------------------------------------------------|
| Step 1   | Count valid RSID + CLNSIG entries               | 2.9 million total; 188k containing “Pathogenic”          |
| Step 2   | Export ClinVar pathogenic RSIDs                 | `clinvar_pathogenic_rs_clnsig_2025jul24_v.01.tsv`        |
| Step 3   | Normalize array rsIDs                           | `dbSNP_variants_on_array_numeric.ids.txt`               |
| Step 4   | Match ClinVar pathogenic RSIDs to array content | `pathogenic_rsID_on_array_2025jul24_v.01.tsv`            |

This workflow establishes a filtered, array-aware list of pathogenic germline variants based on publicly available ClinVar data.
