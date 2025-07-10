# Jones Lab: ClinVar Annotation via rsID Matching

**Author**: Austin Esparza\
**Date**: 2025-07-09\
**Project**: Germline Variant Annotation for Infinium Global Diversity Array

---

## Objective

This workflow annotates germline SNPs from the Infinium Global Diversity Array using rsID-based matching to ClinVar. The goal is to capture clinical significance classifications (CLNSIG) for each array variant by aligning dbSNP rsIDs with entries in the ClinVar database. Unlike strict allele-level matching (CHROM, POS, REF, ALT), this approach uses rsIDs to achieve broader annotation coverage while accepting the trade-off of potential ambiguity due to rsID-level aggregation.

---

## Summary of Workflow

### Input Data

- `dbSNP_variants_on_array_cleaned.tsv`: Cleaned variant list from the array; contains CHROM, POS, REF, ALT, and rsID.
- `clinvar_20250706_chr13_14_16_17_with_CLNSIG.tsv`: Subset of ClinVar (GRCh37) records from chr13, 14, 16, 17 with fields: CHROM, POS, REF, ALT, ID, CLNSIG.

### Output

- `AEsparza_JonesLab_ArrayAnnotated_RSOnly_2025jul09_v.01.vcf`: Final annotated VCF with rsID-level CLNSIG.
- `AEsparza_JonesLab_ClinVar_RSOnly_ClassSummary_2025jul09_v.01.tsv`: Clinical classification frequency table.

---

## Supplement: File Reference Key

| File Path                                                                      | Description                                                                       |
| ------------------------------------------------------------------------------ | --------------------------------------------------------------------------------- |
| `~/Downloads/dbSNP_variants_on_array_cleaned.tsv`                              | Cleaned input variant table from the array with CHROM, rsID, POS, REF, ALT.       |
| `~/Downloads/dbsnp_rsID_map.tsv`                                               | Reformatted mapping of CHROM, POS, and rsID from array TSV for use with bcftools. |
| `~/Downloads/dbsnp_rsID_map.tsv.gz`                                            | Compressed version of the rsID map for tabix indexing.                            |
| `~/Downloads/AEsparza_JonesLab_ArrayConvertedVCF_2025jul09_v.01.vcf`           | Initial VCF generated from array TSV, placeholder ID field.                       |
| `~/Downloads/header_with_contigs.tmp`                                          | Custom VCF header including contig declarations for compatibility with bcftools.  |
| `~/Downloads/array_body.tmp`                                                   | Body of the array-converted VCF, used to rebuild a valid VCF file.                |
| `~/Downloads/AEsparza_JonesLab_ArrayConvertedVCF_2025jul09_v.04.vcf`           | Corrected VCF with header plus array body.                                        |
| `~/Downloads/AEsparza_JonesLab_ArrayWithRSIDs_2025jul09_v.01.vcf`              | VCF with rsIDs injected into the ID field using bcftools.                         |
| `~/Downloads/clinvar_20250706_chr13_14_16_17_with_CLNSIG.tsv`                  | ClinVar subset file with CHROM, POS, REF, ALT, ID, and CLNSIG fields.             |
| `~/Downloads/clinvar_rsID_annotation.tsv`                                      | Reformatted table for annotation with CHROM, POS, rsID, and CLNSIG.               |
| `~/Downloads/clinvar_rsID_annotation_noheader.tsv.gz`                          | Compressed, headerless version of ClinVar annotation table.                       |
| `~/Downloads/header_clnsig_rsID_only.tmp`                                      | bcftools-compatible header that defines the CLNSIG\_from\_rsID INFO tag.          |
| `~/Downloads/AEsparza_JonesLab_ArrayAnnotated_RSOnly_2025jul09_v.01.vcf`       | Final output VCF with CLNSIG\_from\_rsID annotations.                             |
| `~/Downloads/AEsparza_JonesLab_ClinVar_RSOnly_ClassSummary_2025jul09_v.01.tsv` | Summary table of CLNSIG values across all array variants.                         |

---

## Step-by-Step Workflow

### 1. Create Standardized rsID Annotation Map from dbSNP Array

```bash
awk 'NR>1 {gsub(/^chr/, "", $1); print $1"\t"$3"\t"$2}' \
  ~/Downloads/dbSNP_variants_on_array_cleaned.tsv \
  > ~/Downloads/dbsnp_rsID_map.tsv

bgzip -c ~/Downloads/dbsnp_rsID_map.tsv > ~/Downloads/dbsnp_rsID_map.tsv.gz
tabix -s1 -b2 -e2 ~/Downloads/dbsnp_rsID_map.tsv.gz
```

### 2. Generate Valid VCF from Array TSV

```bash
awk 'BEGIN {OFS="\t"; print "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}
     NR>1 {gsub(/^chr/, "", $1); print $1, $3, ".", $4, $5, ".", "PASS", "."}' \
  ~/Downloads/dbSNP_variants_on_array_cleaned.tsv \
  > ~/Downloads/AEsparza_JonesLab_ArrayConvertedVCF_2025jul09_v.01.vcf
```

### 3. Inject rsIDs into VCF ID Field

```bash
printf "##fileformat=VCFv4.2\n##contig=<ID=13>\n##contig=<ID=14>\n##contig=<ID=16>\n##contig=<ID=17>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" > ~/Downloads/header_with_contigs.tmp

grep -v "^#" ~/Downloads/AEsparza_JonesLab_ArrayConvertedVCF_2025jul09_v.01.vcf > ~/Downloads/array_body.tmp

cat ~/Downloads/header_with_contigs.tmp ~/Downloads/array_body.tmp > ~/Downloads/AEsparza_JonesLab_ArrayConvertedVCF_2025jul09_v.04.vcf

bcftools annotate \
  -a ~/Downloads/dbsnp_rsID_map.tsv.gz \
  -c CHROM,POS,ID \
  -o ~/Downloads/AEsparza_JonesLab_ArrayWithRSIDs_2025jul09_v.01.vcf \
  -O v \
  ~/Downloads/AEsparza_JonesLab_ArrayConvertedVCF_2025jul09_v.04.vcf
```

### 4. Create ClinVar Annotation Table Using rsIDs

```bash
cut -f1,2,5,6 ~/Downloads/clinvar_20250706_chr13_14_16_17_with_CLNSIG.tsv \
  | awk 'NR==1 {print "CHROM\tPOS\tID\tCLNSIG"; next} {print $1"\t"$2"\trs"$3"\t"$4}' \
  > ~/Downloads/clinvar_rsID_annotation.tsv

tail -n +2 ~/Downloads/clinvar_rsID_annotation.tsv \
  | bgzip -c > ~/Downloads/clinvar_rsID_annotation_noheader.tsv.gz
tabix -s1 -b2 -e2 ~/Downloads/clinvar_rsID_annotation_noheader.tsv.gz
```

### 5. Annotate VCF with CLNSIG from ClinVar

```bash
printf "##fileformat=VCFv4.2\n##contig=<ID=13>\n##contig=<ID=14>\n##contig=<ID=16>\n##contig=<ID=17>\n##INFO=<ID=CLNSIG_from_rsID,Number=.,Type=String,Description=\"ClinVar clinical significance from rsID match\">\n" > ~/Downloads/header_clnsig_rsID_only.tmp

bcftools annotate \
  -a ~/Downloads/clinvar_rsID_annotation_noheader.tsv.gz \
  -c CHROM,POS,ID,INFO/CLNSIG_from_rsID \
  -h ~/Downloads/header_clnsig_rsID_only.tmp \
  -o ~/Downloads/AEsparza_JonesLab_ArrayAnnotated_RSOnly_2025jul09_v.01.vcf \
  -O v \
  ~/Downloads/AEsparza_JonesLab_ArrayWithRSIDs_2025jul09_v.01.vcf
```

---

## Output Summary Table (CLNSIG Classification Counts)

Run:

```bash
bcftools query -f '%INFO/CLNSIG_from_rsID\n' \
  ~/Downloads/AEsparza_JonesLab_ArrayAnnotated_RSOnly_2025jul09_v.01.vcf \
  | grep -v '^$' | sort | uniq -c | sort -nr
```

### Results:

```
1091 Conflicting_classifications_of_pathogenicity
 925 Pathogenic
 612 Likely_benign
 601 Uncertain_significance
 240 .
 167 Benign
 131 not_provided
  63 Pathogenic/Likely_pathogenic
  59 Benign/Likely_benign
  42 Likely_pathogenic
```

These results total 3,931, matching the number of input array variants.

---

## Interpretation Notes

- This rsID-based annotation strategy provides broad coverage but should not be used for strict pathogenic variant calls.
- Multiple CLNSIG terms per rsID arise from ClinVar's aggregation of conflicting submissions.
- Compared to allele-level annotation (\~68 hits), this method prioritizes recall over precision and is useful for exploratory filtering, enrichment analysis, or manual review.

---

## Reproducibility Notes

- All commands are portable to Unix/Linux/macOS and assume `bcftools >=1.20` and `tabix` from HTSlib.
- Ensure correct formatting of tab-delimited files (no headers unless explicitly allowed).
- Avoid injecting non-standard lines (e.g. `#CHROM` headers) into `-h` files for `bcftools annotate`.

This annotation method complements strict allele-based workflows and allows for flexible downstream prioritization strategies.

