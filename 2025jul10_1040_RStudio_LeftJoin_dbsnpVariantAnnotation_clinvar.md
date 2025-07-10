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


## Evaluation of Annotation Duplication in Joined Dataset

The table below summarizes key statistics used to evaluate annotation inflation in the joined dataset. Because individual rsIDs can map to multiple ClinVar entries, the total number of matched rows can significantly exceed the number of unique variant sites. We assessed this duplication to justify the use of a deduplicated summary for accurate interpretation of clinical significance distributions.

### Summary of rsID Duplication and Deduplication Metrics

| Metric                                                | Value | Description                                                                                                                                                                                          |
| ----------------------------------------------------- | ----- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Total joined rows** (`nrow`)                        | 5,623 | Total number of rows in the joined dataset. Each row represents a match between a dbSNP variant and a ClinVar entry. Multiple rows can exist per rsID due to overlapping or conflicting annotations. |
| **Unique rsID–CLNSIG pairs** (`distinct(ID, CLNSIG)`) | 4,927 | Number of distinct rsID and clinical significance label combinations. This removes exact duplicates but still includes multiple rows per rsID if the clinical labels differ.                         |
| **Unique rsIDs with one annotation** (`slice(1)`)     | 3,931 | Count of unique rsIDs after limiting to one representative annotation per ID. Used for summary-level reporting to avoid inflated totals.                                                             |
| **rsIDs with multiple matched rows**                  | 1,229 | Number of rsIDs that appear more than once in the joined dataset. These entries contribute directly to the inflation of the total row count.                                                         |
| **rsIDs with multiple distinct CLNSIG values**        | 866   | Subset of the above group. These rsIDs are linked to more than one distinct clinical significance label, indicating true interpretive variation across submissions.                                  |

### Interpretation

These results confirm that a substantial proportion of rsIDs (1,229) appear more than once in the joined data, with 866 of those linked to multiple clinical interpretations. Without deduplication, summary statistics would reflect repeated rows rather than unique loci, overstating annotation coverage and misrepresenting category frequencies. For this reason, we report one clinical significance label per rsID in summary tables while retaining the full joined output for downstream analyses that require allele-level detail.


# Clinical Significance Consistency per Variant Site

## Join Logic and Purpose

This section evaluates annotation consistency at the variant-site level, using a composite join key (`CHROM + POS + ID`) to ensure each genomic position is interpreted in the context of its full site. The goal is to determine whether all matched annotations at a given site agree on a "Pathogenic" interpretation.

We define **strictly pathogenic sites** as those where every ClinVar annotation entry is labeled `Pathogenic`, with no conflicting interpretations such as *Benign*, *Likely\_benign*, *Uncertain\_significance*, or `NA`.

---

## Grouping and Classification Logic

- **Join logic:** `CHROM + POS + ID` (applied during the dbSNP–ClinVar merge)
- **Grouping logic:** `group_by(CHROM, POS, ID)`
- **Classification criteria:**
  - **Strictly Pathogenic (singleton):** one matching row labeled `Pathogenic`
  - **Strictly Pathogenic (multi-row):** multiple rows, all consistently labeled `Pathogenic`
  - **Conflicted / Mixed:** one or more rows with a non-Pathogenic label

---

## R Code Implementation

```r
rsid_summary_pos <- matched_chr_pos_id_left %>%
  group_by(CHROM, POS, ID) %>%
  summarise(
    row_count = n(),
    clnsig_count = n_distinct(CLNSIG),
    all_pathogenic = all(CLNSIG == "Pathogenic", na.rm = TRUE)
  ) %>%
  ungroup()

exclusive_pathogenic_singleton <- rsid_summary_pos %>%
  filter(row_count == 1 & all_pathogenic == TRUE)

exclusive_pathogenic_multirow <- rsid_summary_pos %>%
  filter(row_count > 1 & all_pathogenic == TRUE)

conflicted_pathogenic <- rsid_summary_pos %>%
  filter(all_pathogenic == FALSE)

n_exclusive_single <- nrow(exclusive_pathogenic_singleton)  # 2105
n_exclusive_multi  <- nrow(exclusive_pathogenic_multirow)   # 139
n_conflicted       <- nrow(conflicted_pathogenic)           # 1687
n_total            <- nrow(rsid_summary_pos)                # 3931
```

---

## Classification Summary Table

| Category                        | Count | Description                                                                                |
| ------------------------------- | ----- | ------------------------------------------------------------------------------------------ |
| Strictly Pathogenic (singleton) | 2,105 | Variant sites with a single ClinVar match (one row) labeled `Pathogenic`.                  |
| Strictly Pathogenic (multi-row) | 139   | Variant sites with multiple ClinVar entries, all consistently labeled `Pathogenic`.        |
| Conflicted or mixed annotations | 1,687 | Variant sites with multiple distinct CLNSIG labels (e.g., `Pathogenic + Benign`, or `NA`). |
| **Total unique variant sites**  | 3,931 | All variant sites uniquely defined by `CHROM + POS + ID`.                                  |

---

## Interpretation

This analysis characterizes clinical significance consistency at the variant-site level. It shows that over half (57%) of evaluated sites are exclusively annotated as Pathogenic, while the remaining 43% display variation in interpretation. These results highlight the relevance of using site-level resolution for annotation review and suggest that stricter filtering may be warranted in workflows where clinical confidence is a priority.
---


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

# Supplemental Appendix: File Paths and Directory Structure

## Project Root
```r
setwd("/Users/austinesparza/Downloads/JonesLab/scripts/2025jul10_0900_dbsnpVariant_ClinVar_Annotation_v.01")
project_dir <- "/Users/austinesparza/Downloads/JonesLab/scripts/2025jul10_0900_dbsnpVariant_ClinVar_Annotation_v.01"
```
### File Structures and Critical File Pathways

## Directory Layout
```bash
2025jul10_0900_dbsnpVariant_ClinVar_Annotation_v.01/
├── data_raw/           # Input files (ClinVar + dbSNP TSVs)
│   ├── clinvar_20250706_chr13_14_16_17_with_RS_CLNSIG.tsv
│   └── dbSNP_variants_on_array_cleaned.tsv
├── data_processed/     # Output tables and annotated results
│   ├── AEsparza_JonesLab_ClinVar_Annotation_RSOnly_2025jul10_v.01.tsv
│   ├── AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_LeftJoin_2025jul10_v.01.tsv
│   └── AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_2025jul10_v.02.tsv
├── scripts/            # R scripts or shell utilities (if added)
├── AEsparza_JonesLab_RSessionLog_2025jul10_1113_v.01.txt
├── AEsparza_JonesLab_RCommandHistory_2025jul10_1113.Rhistory
```

## Input Files
| Filename                                                                 | Description                                                                 |
|--------------------------------------------------------------------------|-----------------------------------------------------------------------------|
| `dbSNP_variants_on_array_cleaned.tsv`                                    | Variant list from Global Diversity Array with `CHROM`, `POS`, `REF`, `ALT`, and `rsID` |
| `clinvar_20250706_chr13_14_16_17_with_RS_CLNSIG.tsv`                     | Extracted ClinVar VCF (GRCh37) subset from chromosomes 13, 14, 16, 17        |

## Output Files
| Filename                                                                  | Description                                                            |
|---------------------------------------------------------------------------|------------------------------------------------------------------------|
| `AEsparza_JonesLab_ClinVar_Annotation_RSOnly_2025jul10_v.01.tsv`         | rsID-only left join annotation table                                  |
| `AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_LeftJoin_2025jul10_v.01.tsv`    | CHROM + POS + ID left join annotation table                           |
| `AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_2025jul10_v.02.tsv`             | Inner join version of the above (diagnostic comparison)               |
| `AEsparza_JonesLab_RSessionLog_2025jul10_1113_v.01.txt`                  | Raw RStudio console log (after `sink()` command)                      |
| `AEsparza_JonesLab_RCommandHistory_2025jul10_1113.Rhistory`              | Full R command history saved via `savehistory()`                      |

## Path Initialization Example in R
```r
clinvar_path <- file.path(project_dir, "data_raw", "clinvar_20250706_chr13_14_16_17_with_RS_CLNSIG.tsv")
dbsnp_path <- file.path(project_dir, "data_raw", "dbSNP_variants_on_array_cleaned.tsv")
output_path <- file.path(project_dir, "data_processed", "AEsparza_JonesLab_ClinVar_Annotation_RSOnly_2025jul10_v.01.tsv")
```

