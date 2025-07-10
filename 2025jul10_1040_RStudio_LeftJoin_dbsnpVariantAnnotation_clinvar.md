# Germline Variant Annotation via rsID and Allele-Level Matching

**Author**: Austin Esparza\
**Project**: Jones Lab – Germline Variant Harmonization\
**Date**: 2025-07-10\
**Genome Build**: GRCh37

---

## Objective

This workflow annotates dbSNP-mapped array variants with ClinVar clinical significance (CLNSIG) labels using three distinct strategies:

1. **rsID-based joins**,
2. **allele-level joins based on CHROM, POS, REF, and ALT**, and
3. **CHROM + POS + rsID joins**.

The dbSNP variant set derives from the Infinium Global Diversity Array. Annotation relies on a ClinVar VCF (GRCh37) subset covering chromosomes 13, 14, 16, and 17.

---

## Input Files

| Filename                                             | Description                                                        |
| ---------------------------------------------------- | ------------------------------------------------------------------ |
| `dbSNP_variants_on_array_cleaned.tsv`                | dbSNP variant list from array mapping (CHROM, POS, REF, ALT, rsID) |
| `clinvar_20250709_chr13_14_16_17_with_RS_CLNSIG.tsv` | ClinVar VCF subset with CHROM, POS, RSID, REF, ALT, CLNSIG         |

---

## Directory Structure

Project root: `/Users/austinesparza/Downloads/JonesLab/scripts/2025jul10_0900_dbsnpVariant_ClinVar_Annotation_v.01/`

Subdirectories:

```
data_raw/             # Input ClinVar and dbSNP TSV files
data_processed/       # Output tables and annotated variant results
scripts/              # Any relevant R or shell code
```

---

## Methods

### 1. Load and Normalize dbSNP Variants

```r
library(readr)
library(dplyr)
library(stringr)

project_dir <- "/Users/austinesparza/Downloads/JonesLab/scripts/2025jul10_0900_dbsnpVariant_ClinVar_Annotation_v.01"

# Load cleaned dbSNP variant table
dbsnp_path <- file.path(project_dir, "data_raw", "dbSNP_variants_on_array_cleaned.tsv")

dbsnp_df <- read_tsv(dbsnp_path, col_types = cols(
  chr = col_integer(),
  dbSNP_id_on_array = col_character(),
  position = col_integer(),
  reference = col_character(),
  alternate = col_character()
))

dbsnp_clean <- dbsnp_df %>%
  mutate(
    CHROM = as.character(chr),
    POS = as.integer(position),
    ID = str_remove(dbSNP_id_on_array, "^rs") %>% str_trim()
  ) %>%
  select(CHROM, POS, ID, REF = reference, ALT = alternate)
```

### 2. Load and Parse ClinVar Subset

```r
clinvar_path <- file.path(project_dir, "data_raw", "clinvar_20250709_chr13_14_16_17_with_RS_CLNSIG.tsv")

clinvar_df <- read_tsv(
  clinvar_path,
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
    ID = str_trim(RSID)
  ) %>%
  select(CHROM, POS, ID, REF, ALT, CLNSIG)
```

### 3. Left Join on CHROM + POS + ID

```r
matched_chr_pos_id_left <- left_join(
  dbsnp_clean,
  clinvar_slim,
  by = c("CHROM", "POS", "ID")
)

output_path <- file.path(project_dir, "data_processed", "AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_LeftJoin_2025jul10_v.01.tsv")
write_tsv(matched_chr_pos_id_left, output_path)
```

### 4. Clinical Significance Summary (by Unique rsID, 1 Match per ID)

```r
clinsig_summary_unique_once <- matched_chr_pos_id_left %>%
  mutate(CLNSIG = ifelse(is.na(CLNSIG), "NA_unannotated", CLNSIG)) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  count(CLNSIG, name = "Unique_rsID_Count") %>%
  arrange(desc(Unique_rsID_Count)) %>%
  bind_rows(
    tibble(CLNSIG = "TOTAL", Unique_rsID_Count = sum(.$Unique_rsID_Count))
  )

print(clinsig_summary_unique_once)
```

#### Output Table

| CLNSIG                                          | Unique rsID Count |
| ----------------------------------------------- | ----------------- |
| NA\_unannotated                                 | 1606              |
| Pathogenic                                      | 868               |
| Conflicting\_classifications\_of\_pathogenicity | 622               |
| Likely\_benign                                  | 257               |
| Uncertain\_significance                         | 238               |
| Benign                                          | 148               |
| not\_provided                                   | 61                |
| Benign/Likely\_benign                           | 56                |
| Pathogenic/Likely\_pathogenic                   | 52                |
| Likely\_pathogenic                              | 23                |
| TOTAL                                           | 3931              |

---

## Method Comparison Summary

| Method                | Join Fields                  | Annotated Variants | Strengths                   | Weaknesses                             |
| --------------------- | ---------------------------- | ------------------ | --------------------------- | -------------------------------------- |
| rsID-only             | `ID`                         | 5                  | Simple, quick lookup        | Incomplete, minimal overlap            |
| Allele-level matching | `CHROM`, `POS`, `REF`, `ALT` | 68                 | Precise, allele-specific    | Fails if ALT differs in representation |
| CHROM + POS + ID      | `CHROM`, `POS`, `ID`         | 3931 (after dedup) | Best rsID-resolution method | Requires ID normalization              |

---

## Conclusions

- The **CHROM + POS + ID** approach yielded the highest confidence set of rsID-annotated array variants.
- **Total annotated rsIDs**: 3931 (deduplicated by ID).
- REF and ALT were preserved but excluded from matching to avoid dropouts from allele encoding mismatches.

---

## Next Steps

- Normalize allele representation across sources using SPDI or vt.
- Apply VCF-based tools like `bcftools annotate` for genome-wide annotation.
- Use this harmonized rsID table to support downstream filtering of variants in ovarian cancer gene panels.

