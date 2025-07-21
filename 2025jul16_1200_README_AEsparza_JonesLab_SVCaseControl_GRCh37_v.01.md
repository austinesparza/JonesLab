# Structural Variant Exploration – Germline Case-Control Comparison (GRCh37)

## Project Overview

This project aims to characterize germline structural variant (SV) patterns using DRAGEN-called data from a case-control cohort, provided by Nimisha. The data are aligned to the GRCh37 reference genome and focus on deletions, duplications, inversions, and unresolved breakends (BNDs). This analysis is conducted in support of downstream work associated with the Infinium Global Diversity Array and requires consistent preprocessing and standardized outputs.



## Assignment Scope

Nimisha provided two DRAGEN-called SV files representing case and control cohorts. The following objectives were specified:

1. Identify chromosomes with the highest and lowest total SV counts. Determine whether these observations are consistent with previously reported germline variant patterns.
2. Calculate the full length distribution of SVs, including deletions, duplications, and inversions. Compare length trends between case and control groups.
3. Include unresolved variants (breakends, or BNDs) in burden counts. Consider potential reasons for lack of resolution and evaluate chromosome-level clustering.
4. Recommend appropriate visualizations for summarizing the above trends.

All work must remain GRCh37-specific to ensure compatibility with the Global Diversity Array platform.

# Jones Lab SV Standardization Summary  
**Project:** *Discovery of Rare Germline Structural Variants in Epithelial Ovarian Cancer*  
**Date:** July 17, 2025  
**Analyst:** Austin Esparza  

# Requirements
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
seaborn>=0.11.0
scipy>=1.7.0
platformdirs>=2.0.0
python-dateutil>=2.8.0

## Environment Notes
All analysis was performed on a local workstation (macOS 10.15.7) using Python 3.9. Packages were managed via `pip`. Scripts are compatible with Unix-style paths and assume availability of standard command-line utilities (e.g., `cat`, `wc`, `sort`) for optional audit steps.

## Directory Structure

```
SV_Exploration_2025jul15/
├── data_raw/
│ ├── all_cases.dragen.sv.txt
│ └── all_controls.dragen.sv.txt
├── data_processed/
│ ├── AEsparza_JonesLab_CleanedCasesSV_.tsv
│ └── AEsparza_JonesLab_CleanedControlsSV_.tsv
├── results/
│ ├── AEsparza_JonesLab_ChromSVCounts_.tsv
│ ├── AEsparza_JonesLab_SVTypeSummary_.tsv
│ ├── plots/
│ │ └── AEsparza_JonesLab_SVLEN_ViolinPlot_.png
│ ├── tables/
│ │ ├── AEsparza_JonesLab_SVLengthRange_.tsv
│ │ ├── AEsparza_JonesLab_SVLengthStats_.tsv
│ │ ├── AEsparza_JonesLab_TopBottomChromSVs_.tsv
│ │ └── AEsparza_JonesLab_ChromSVTypeSummary_.tsv
│ └── dataframes/
│ ├── AEsparza_JonesLab_AllSVs_Cleaned_.tsv
│ ├── AllSVs_Case/
│ │ └── AEsparza_JonesLab_SVCatalog_case_chr*.tsv
│ └── AllSVs_Control/
│ └── AEsparza_JonesLab_SVCatalog_control_chr*.tsv
├── scripts/
│ ├── AEsparza_JonesLab_ParseCasesSV_.py
│ ├── AEsparza_JonesLab_ParseControlsSV_.py
│ ├── AEsparza_JonesLab_SVSummary_.py
│ ├── AEsparza_JonesLab_SVStats_Tables_.py
│ ├── AEsparza_JonesLab_SVLengthRange_.py
│ └── AEsparza_JonesLab_SplitSVs_ByGroupAndChrom_.py
└── logs/

11 directories, 85 files

```
---

## Parsing & Standardization Workflow

### Input Files
- `data_raw/all_cases.dragen.sv.txt`
- `data_raw/all_controls.dragen.sv.txt`

### Scripts Executed
- `scripts/AEsparza_JonesLab_ParseCasesSV_2025jul17_v.03.py`
- `scripts/AEsparza_JonesLab_ParseControlsSV_2025jul17_v.02.py`

### Transformations Applied
- Parsed columns: `sample_id`, `chrom`, `start_pos`, `end_pos`, `sv_type`, `sv_len`, `format_field`
- Computed fields:
  - `sv_len_numeric`: Numeric conversion of SV length
  - `sv_len_abs`: Absolute value of SV length
- Metadata fields added:
  - `caller = DRAGEN`
  - `build = GRCh37`
  - `source_group = case` or `control`
- Output format: Cleaned `.tsv` files with consistent column names and structure

### Post-Cleaning Transformations
- Merged cleaned case/control TSVs into a unified master file
- Calculated:
  - `log10_sv_len_abs` (base-10 log of absolute SV length)
  - Added `source_group` labels (case or control)
- Split master table by group → chromosome using:
  - `scripts/AEsparza_JonesLab_SplitSVs_ByGroupAndChrom_2025jul22_v.01.py`
- Result: One `.tsv` per chromosome × group in `results/dataframes/AllSVs_[Group]/`

## Output Naming Convention

All files follow this convention:

`AEsparza_JonesLab_[ContentDescription]_YYYYmonDD_v.##.[tsv|png|py]`

- Example: `AEsparza_JonesLab_SVCatalog_case_chr17_2025jul22_v.01.tsv`
- This ensures reproducibility, version tracking, and auditability.

---

## Sample Representation Summary

| Dataset                             | Output Filename                                              | Unique `sample_id` Count |
|-------------------------------------|--------------------------------------------------------------|---------------------------|
| **Cases**                           | `AEsparza_JonesLab_CleanedCasesSV_2025jul17_v.03.tsv`        | 1,016                     |
| **Controls**                        | `AEsparza_JonesLab_CleanedControlsSV_2025jul17_v.02.tsv`     | 2,945                     |

## Summary Statistics Tables

| Filename                                                                 | Description                                 |
|--------------------------------------------------------------------------|---------------------------------------------|
| `AEsparza_JonesLab_SVTypeSummary_2025jul15_v.01.tsv`                     | Raw type distribution summary               |
| `AEsparza_JonesLab_ChromSVTypeTotals_WithAndWithoutBND_2025jul17_v.01.tsv` | SV type totals per chromosome (with/without BNDs) |
| `AEsparza_JonesLab_ChromSVTypeSummary_2025jul17_v.01.tsv`               | Burden stratified by chromosome and type    |
| `AEsparza_JonesLab_SVLengthStats_AllTypes_2025jul17_v.01.tsv`           | Length stats (min, max, median, etc.)       |
| `AEsparza_JonesLab_SVLengthRange_DEL_DUP_INV_2025jul17_v.01.tsv`        | Raw and log-scaled length ranges            |
| `AEsparza_JonesLab_TopBottomChromSVs_2025jul17_v.01.tsv`                | Top/bottom 3 chromosomes by burden          |

---

## Analytical Objectives

The following tasks are planned based on the cleaned and harmonized datasets:

1. **Burden Analysis by Chromosome**
   - Total SV counts by chromosome, disaggregated by SV type.
   - Focus on top three and bottom three chromosomes by count.

2. **SV Length Characterization**
   - Length range, median, and distribution plots for each SV type.
   - Case versus control comparisons.

3. **Breakend (BND) Analysis**
   - Chromosome-specific burden of unresolved SVs.
   - Exploration of patterns that may explain lack of resolution.

4. **Visualization and Reporting**
   - Violin and histogram plots of SV lengths.
   - Bar plots of SV count per chromosome.
   - All visualizations stored in `results/plots/`.

## Next Steps (As of July 22, 2025)

- [x] Clean and standardize DRAGEN-called SV files
- [x] Merge case/control data with length and log transforms
- [x] Split full catalog by chromosome × group
- [x] Compute summary tables of SV burden and length ranges
- [x] Visualize SV length distributions using violin plots
- [ ] Visualize per-chromosome SV burden using bar plots
- [ ] Evaluate BND resolution failures (preliminary clustering analysis)
- [ ] Apply statistical comparisons between case/control distributions (e.g., KS test, Wilcoxon)


## Notes

- All coordinate-based work must remain tied to GRCh37.
- Visualization choices should prioritize interpretability and publication-readiness.
- Breakends should be handled as a distinct category in all summaries.
- Versioned outputs are required at every stage for reproducibility.

## Supplemental Context and Analytical Rationale

This analysis was initiated to investigate structural variant (SV) patterns in a germline case-control dataset generated using the DRAGEN variant caller and aligned to the GRCh37 reference genome. The provided SV callsets, one for cases and one for controls, contain chromosomal coordinates, SV type, and length data for events including deletions, duplications, inversions, insertions, and unresolved breakends (BNDs). The immediate goal is to produce a foundational characterization of SV burden and length distribution in support of downstream variant prioritization and visualization.

The initial analytical scope is defined by the following tasks:

- Identify chromosomes with the highest and lowest counts of SVs, and determine whether these distributions reflect previously reported germline variant patterns.
- Analyze the full length range of deletions and duplications, and assess whether their distributions differ between cases and controls.
- Extend this evaluation to include inversions, and determine whether similar trends are observed across groups.
- Examine the distribution of unresolved breakends (BNDs), treating these as a separate category. This includes quantifying their chromosomal burden and considering possible reasons for partial resolution during variant calling.

The aim of this work is to inform the design of scalable, reproducible SV analysis workflows that are compatible with other GRCh37-based datasets, including genotyping arrays such as the Infinium Global Diversity Array. All coordinate systems and annotations remain GRCh37-specific. Visualization strategies are selected to prioritize clarity, interpretability, and publication-quality output. This document will be expanded as additional input files are incorporated or as analytical priorities shift.


## Document History

| Version | Date       | Updates Summary                                                              |
|---------|------------|-------------------------------------------------------------------------------|
| v.01    | 2025-07-17 | Initial documentation of parsing workflow and objectives                     |
| v.02    | 2025-07-18 | Added merged master table and violin plot generation                         |
| v.03    | 2025-07-22 | Split by group/chromosome, added SV type summaries, log-scaled lengths, etc. |

**Current Version:** v.03  
**Maintainer:** Austin Esparza  
