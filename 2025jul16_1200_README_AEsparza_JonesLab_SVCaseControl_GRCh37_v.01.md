# Structural Variant Exploration – Germline Case-Control Comparison (GRCh37)

## Project Overview

This project aims to characterize germline structural variant (SV) patterns using DRAGEN-called data from a case-control cohort, provided by Nimisha. The data are aligned to the GRCh37 reference genome and focus on deletions, duplications, inversions, and unresolved breakends (BNDs). This analysis is conducted in support of downstream work associated with the Infinium Global Diversity Array and requires consistent preprocessing and standardized outputs.

## Directory Structure

```
SV_Exploration_2025jul15/
├── data_raw/
│   ├── all_cases.dragen.sv.txt
│   └── all_controls.dragen.sv.txt
├── data_processed/
│   ├── AEsparza_JonesLab_CleanedCasesSV_2025jul16_v.01.tsv
│   └── AEsparza_JonesLab_CleanedControlsSV_2025jul16_v.01.tsv
├── results/
│   └── plots/
├── scripts/
│   ├── AEsparza_JonesLab_ParseCasesSV_2025jul16_v.02.py
│   └── AEsparza_JonesLab_ParseControlsSV_2025jul16_v.01.py
└── logs/
```

## Assignment Scope

Nimisha provided two DRAGEN-called SV files representing case and control cohorts. The following objectives were specified:

1. Identify chromosomes with the highest and lowest total SV counts. Determine whether these observations are consistent with previously reported germline variant patterns.
2. Calculate the full length distribution of SVs, including deletions, duplications, and inversions. Compare length trends between case and control groups.
3. Include unresolved variants (breakends, or BNDs) in burden counts. Consider potential reasons for lack of resolution and evaluate chromosome-level clustering.
4. Recommend appropriate visualizations for summarizing the above trends.

All work must remain GRCh37-specific to ensure compatibility with the Global Diversity Array platform.

## Preprocessing Workflow

The case and control files were independently parsed and standardized using dedicated Python scripts. Both raw files were converted into clean tab-delimited tables with the following structure:

| Column Name     | Description                                            |
|------------------|--------------------------------------------------------|
| `sample_id`      | Name of the originating sample file                    |
| `chrom`          | Chromosome (e.g., chr1, chr2, etc.)                    |
| `start_pos`      | SV start coordinate in GRCh37                          |
| `end_pos`        | SV end coordinate in GRCh37                            |
| `sv_type`        | SV class (DEL, DUP, INV, INS, BND, etc.)              |
| `sv_len`         | Reported SV length (signed)                            |
| `sv_len_abs`     | Absolute SV length (for length-based analysis)         |
| `source_group`   | Case or control indicator                              |
| `caller`         | Variant caller (DRAGEN)                                |
| `build`          | Genome reference build (GRCh37)                        |

### Case File Processing

- Input: `data_raw/all_cases.dragen.sv.txt`
- Output: `data_processed/AEsparza_JonesLab_CleanedCasesSV_2025jul16_v.01.tsv`
- Script: `scripts/AEsparza_JonesLab_ParseCasesSV_2025jul16_v.02.py`

### Control File Processing

- Input: `data_raw/all_controls.dragen.sv.txt`
- Output: `data_processed/AEsparza_JonesLab_CleanedControlsSV_2025jul16_v.01.tsv`
- Script: `scripts/AEsparza_JonesLab_ParseControlsSV_2025jul16_v.01.py`

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

## Next Steps

- Implement core summary statistics for DEL, DUP, and INV SVs.
- Perform per-chromosome burden analysis and generate bar plots.
- Create violin plots for SV length distributions by group and SV type.
- Evaluate BND distribution separately and identify possible artifacts or resolution failures.

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

- Version: `v.01`
- Author: Austin Esparza
- Last Updated: 2025-07-16
- Mentor: Nimisha Mazumdar
