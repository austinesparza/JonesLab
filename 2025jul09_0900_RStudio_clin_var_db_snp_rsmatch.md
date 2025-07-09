# ClinVar-dbSNP Harmonization by rsID (CHROM, POS, ID)

**Author:** Austin Esparza\
**Project:** JonesLab – Germline Variant Annotation Pipeline\
**Date:** 2025-07-09\
**Genome Build:** GRCh37 / hg19

---

## Objective

Establish a harmonized join between ClinVar and dbSNP datasets based on exact matches at three genomic fields: `CHROM`, `POS`, and `rsID`. This provides a verified list of variant positions on the Infinium Global Diversity Array that are also represented in ClinVar with associated clinical significance annotations.

---

## Input Files

| Filename                                          | Description                                  |
| ------------------------------------------------- | -------------------------------------------- |
| `~/Downloads/clinvar_20250706.vcf.gz`             | Original ClinVar GRCh37 full VCF             |
| `~/Downloads/dbSNP_variants_on_array_cleaned.tsv` | Cleaned dbSNP export with array-mapped rsIDs |

---

## Output Files

| Filename                                                                 | Description                                            |
| ------------------------------------------------------------------------ | ------------------------------------------------------ |
| `~/Downloads/clinvar_20250709_chr13_14_16_17_with_RS_CLNSIG.tsv`         | ClinVar subset (chr13/14/16/17) with RSID and CLNSIG   |
| `~/Downloads/AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_2025jul09_v.01.tsv` | Final merged table of rsID-matched variants (N = 4017) |

---

## Workflow Summary

### 1. Index ClinVar VCF for random access

```bash
tabix -p vcf ~/Downloads/clinvar_20250706.vcf.gz
```

### 2. Subset to chromosomes 13, 14, 16, 17

```bash
bcftools view \
  -r 13,14,16,17 \
  -Oz \
  -o ~/Downloads/clinvar_20250706_chr13_14_16_17.vcf.gz \
  ~/Downloads/clinvar_20250706.vcf.gz
```

### 3. Extract CHROM, POS, REF, ALT, RSID, CLNSIG from ClinVar

```bash
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RS\t%CLNSIG\n' \
  ~/Downloads/clinvar_20250706_chr13_14_16_17.vcf.gz \
  > ~/Downloads/clinvar_20250709_chr13_14_16_17_with_RS_CLNSIG.tsv
```

---

## Data Loading and Cleaning (R)

```r
library(readr)
library(dplyr)
library(stringr)
```

### Load ClinVar subset with RSIDs

```r
clinvar_df <- read_tsv(
  "~/Downloads/clinvar_20250709_chr13_14_16_17_with_RS_CLNSIG.tsv",
  col_names = c("CHROM", "POS", "REF", "ALT", "RSID", "CLNSIG"),
  col_types = cols(
    CHROM = col_integer(),
    POS = col_integer(),
    REF = col_character(),
    ALT = col_character(),
    RSID = col_character(),
    CLNSIG = col_character()
  )
)

clinvar_slim <- clinvar_df %>%
  mutate(
    CHROM = as.character(CHROM),
    POS = as.integer(POS),
    ID = str_trim(RSID)
  ) %>%
  select(CHROM, POS, ID, REF, ALT, CLNSIG)
```

### Load and normalize dbSNP array-mapped variant list

```r
dbsnp_df <- read_tsv(
  "~/Downloads/dbSNP_variants_on_array_cleaned.tsv",
  col_types = cols(
    chr = col_integer(),
    dbSNP_id_on_array = col_character(),
    position = col_integer(),
    reference = col_character(),
    alternate = col_character()
  )
)

dbsnp_clean <- dbsnp_df %>%
  mutate(
    CHROM = as.character(chr),
    POS = as.integer(position),
    ID = str_remove(dbSNP_id_on_array, "^rs") %>% str_trim()
  ) %>%
  select(CHROM, POS, ID, REF = reference, ALT = alternate)
```

---

## Final Join by rsID Match

```r
matched_chr_pos_id <- inner_join(
  clinvar_slim,
  dbsnp_clean,
  by = c("CHROM", "POS", "ID")
)
```

### Export Final Match Table

```r
write_tsv(
  matched_chr_pos_id,
  "~/Downloads/AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_2025jul09_v.01.tsv"
)
```

---

## Clinical Significance Summary

### By Total Annotations

```r
clinsig_summary <- matched_chr_pos_id %>%
  count(CLNSIG, name = "Count") %>%
  arrange(desc(Count))
```

### By Unique rsID

```r
clinsig_by_rsid <- matched_chr_pos_id %>%
  distinct(ID, CLNSIG) %>%
  count(CLNSIG, name = "Unique_rsID_Count") %>%
  arrange(desc(Unique_rsID_Count))
```

### Result Table

```
| Clinical Significance                        | Unique rsID Count |
|:--------------------------------------------|------------------:|
| Pathogenic                                   |              1054 |
| Conflicting_classifications_of_pathogenicity |               885 |
| Uncertain_significance                       |               448 |
| Likely_benign                                |               359 |
| Benign                                       |               204 |
| not_provided                                 |               126 |
| Pathogenic/Likely_pathogenic                 |               109 |
| Benign/Likely_benign                         |                86 |
| Likely_pathogenic                            |                50 |
```

---

## Decision and Notes

- Matching was performed on `CHROM`, `POS`, and `ID` (rsID numeric only).
- REF and ALT fields were retained but not used for matching, due to inconsistent ALT encodings (`0`, missing minor alleles) in dbSNP.
- A total of 4,017 rsIDs were matched between ClinVar and the array-mapped dbSNP dataset.
- This file serves as the authoritative harmonized rsID set for downstream clinical significance filtering and reporting.
- No header or contig parsing errors remain.
- The current analysis matched 4,017 entries based on `CHROM`, `POS`, and `rsID`, without filtering on `REF` and `ALT` due to allele formatting inconsistencies in the variant data from Nimisha.

---

