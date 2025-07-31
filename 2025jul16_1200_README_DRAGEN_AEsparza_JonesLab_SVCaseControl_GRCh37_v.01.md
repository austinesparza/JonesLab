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
   - Bar plots (stacked?) of SV count per chromosome.
   - All visualizations stored in `results/plots/`.


##  AEsparza_JonesLab_ChromSVBarplot_case_2025jul28_v.01.png
![SV Burden by Chrom and Type CASE](https://github.com/austinesparza/JonesLab/blob/AEsparza_WorkingFiles_Tables_Images/AEsparza_JonesLab_ChromSVBarplot_case_2025jul28_v.01.png)

##  AEsparza_JonesLab_ChromSVBarplot_control_2025jul28_v.01.png
![SV Burden by Chrom and Type CONTROL](https://github.com/austinesparza/JonesLab/blob/445b3cbc2925f3c8357d502348291b33566253a5/AEsparza_JonesLab_ChromSVBarplot_control_2025jul28_v.01.png)

   
## Next Steps (As of July 29, 2025)

- [x] Clean and standardize DRAGEN-called SV files
- [x] Merge case/control data with length and log transforms
- [x] Split full catalog by chromosome × group
- [x] Compute summary tables of SV burden and length ranges
- [x] Visualize SV length distributions using violin plots
- [x] Visualize per-chromosome SV burden using bar plots
- [ ] Evaluate BND resolution failures (preliminary clustering analysis)
- [x] Apply statistical comparisons between case/control distributions (e.g., KS test, Wilcoxon)


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

## KS Test Summary: Structural Variant Length Distributions (Cases vs Controls)
 
**Script:** `AEsparza_JonesLab_KSTest_SVLen_CaseControl_2025jul28_v.01.py`  
**Directory:** `scripts/`  
**Output:** `results/stats/AEsparza_JonesLab_KSTest_SVLen_CaseControl_2025jul28_v.01.tsv`  
**Purpose:** Assess whether the distribution of SV lengths differs between cases and controls using the two-sample Kolmogorov–Smirnov (KS) test.

---

### Objective

To statistically compare the empirical distributions of structural variant lengths (`sv_len_abs`) between case and control cohorts across three SV types: deletions (DEL), duplications (DUP), and inversions (INV). Breakends (BND) were excluded from this analysis due to lack of resolution and ambiguous lengths.

---

### Methods

1. **Data Input**
   - **Cases File:** `data_processed/AEsparza_JonesLab_CleanedCasesSV_2025jul17_v.03.tsv`
   - **Controls File:** `data_processed/AEsparza_JonesLab_CleanedControlsSV_2025jul17_v.02.tsv`
   - Input files must include:
     - `sv_type`: structural variant type (e.g., DEL, DUP, INV)
     - `sv_len_abs`: absolute length of the variant

2. **Filtering Logic**
   - For each SV type (`DEL`, `DUP`, `INV`):
     - Subset lengths using `df[df["sv_type"] == sv]["sv_len_abs"].dropna()`
     - Skip SV type if either group is empty

3. **KS Test Execution**
   - Performed with `scipy.stats.ks_2samp` to test the null hypothesis that the distributions of variant lengths are the same in cases and controls
   - Reports:
     - KS statistic (magnitude of distributional difference)
     - P-value (statistical significance)
     - Sample sizes per group

4. **Output**
   - Results saved in TSV format:
     - Columns: `sv_type`, `ks_statistic`, `p_value`, `n_case`, `n_control`
     - File: `results/stats/AEsparza_JonesLab_KSTest_SVLen_CaseControl_2025jul28_v.01.tsv`

---

### KS Output

| sv_type | ks_statistic | p_value     | n_case  | n_control |
|---------|--------------|-------------|---------|-----------|
| DEL     | 0.0193       | 0E+00       | 4840312 | 14760659  |
| DUP     | 0.0263       | 8.3785E-23  | 46903   | 181494    |


---

### Interpretation

- The extremely low p-values indicate a statistically significant difference in SV length distributions between cases and controls for both DEL and DUP.
- The KS statistics are small, suggesting modest effect sizes—appropriate given large sample sizes.
- This supports further characterization or stratification by genomic context or clinical phenotype, and suggests SV size distributions may offer subtle but real biological signal between groups.

---

## Mann–Whitney U Test Summary: SV Lengths by Group  
**Script**: `AEsparza_JonesLab_WilcoxonTest_SVLen_CaseControl_2025jul28_v.01.py`  
**Objective**:  
To compare the distribution of absolute structural variant (SV) lengths between cases and controls across three SV types (`DEL`, `DUP`, `INV`) using the non-parametric Mann–Whitney U test (Wilcoxon rank-sum).  

---

### Inputs  

| Filepath | Description |
|----------|-------------|
| `data_processed/AEsparza_JonesLab_CleanedCasesSV_2025jul17_v.03.tsv` | Cleaned SV call table for case samples |
| `data_processed/AEsparza_JonesLab_CleanedControlsSV_2025jul17_v.02.tsv` | Cleaned SV call table for control samples |

Both input tables must include the fields:  
- `sv_type`: categorical (e.g. DEL, DUP, INV)  
- `sv_len_abs`: numeric absolute SV length (precomputed, required for test)

---

### Workflow  

1. **Load data**:  
   Case and control files are read using `pandas.read_csv()`, with tab-delimited format.

2. **Loop through SV types**:  
   For each SV type of interest:
   - Subset lengths for that type in cases and controls
   - Drop missing (`NaN`) values
   - Skip test if either group is empty

3. **Perform test**:  
   Run `scipy.stats.mannwhitneyu()` with `alternative='two-sided'`.

4. **Record results**:  
   Save `u_statistic`, `p_value`, and sample sizes per group to a structured results dictionary.

5. **Export results**:  
   Output is written to:

   ```
   results/stats/AEsparza_JonesLab_WilcoxonTest_SVLen_CaseControl_2025jul28_v.01.tsv
   ```

---

### Output Columns  

| Column | Description |
|--------|-------------|
| `sv_type` | Structural variant type (DEL, DUP, INV) |
| `u_statistic` | Mann–Whitney U test statistic |
| `p_value` | P-value for difference in SV length distributions |
| `n_case` | Number of case SVs tested |
| `n_control` | Number of control SVs tested |

---

### Wilcoxon Output

| sv_type | u_statistic     | p_value    | n_case  | n_control |
|---------|------------------|------------|---------|-----------|
| DEL     | 35078845316775   | 0E+00      | 4840312 | 14760659  |
| DUP     | 4354776315       | 1.025E-14  | 46903   | 181494    |

---

### Notes  

- The Mann–Whitney U test does not assume normality and is suitable for highly skewed distributions such as SV lengths.
- SV types with no data in either group are skipped with an explicit log message.
- This test complements the Kolmogorov–Smirnov (KS) test by comparing rank distributions rather than cumulative distribution shapes.

---

## Note on Inversions (INV)

Although inversions (`INV`) were included as a structural variant type in both the Kolmogorov–Smirnov and Mann–Whitney U test workflows, no qualifying inversion events were present in the cleaned datasets for either cases or controls. This was verified using the chromosome-level SV summary table (`AEsparza_JonesLab_ChromSVTypeSummary_2025jul17_v.01.tsv`), where all values in the `num_INV` column were `0`. As such, statistical comparisons involving inversions were not performed, and this SV class is excluded from the final results.


**Script:** `AEsparza_JonesLab_SVComparison_NormalizedPlot_2025jul30_v.01.py`  
**Genome Build:** GRCh37  

---

## Objective

Quantify and compare the burden of structural variants (SVs) per chromosome between:

- **case_DRAGEN**: Case samples (n = 1016) processed via DRAGEN-based SV calling
- **gnomAD**: Publicly available structural variant dataset from population controls (n = 14,891)

SV counts are normalized by both sample size and chromosome length to yield:  
**SVs per sample per megabase (Mb)** — a metric allowing intra-cohort chromosomal burden comparison.

**Note:** Due to differences in SV detection algorithms, definitions, and quality filters between pipelines, comparisons across datasets are approximate and **should not be interpreted as absolute biological differences.**

---

## Input Files

**Primary TSV:**  
`AEsparza_JonesLab_SVComparison_NoChrM_CaseGnomAD_2025jul30_v.01.tsv`  
- Fields: `chrom`, `group`, `total_SVs`

**Chromosome Sizes:**  
Reference values for GRCh37 autosomes and sex chromosomes (Mb). Stored internally in script.

---

## Normalization Formula

```python
SVs_per_sample_per_Mb = total_SVs / (sample_count × chrom_size_mb)
```

- `total_SVs`: SV count per chromosome  
- `sample_count`: 14891 for gnomAD; 1016 for DRAGEN  
- `chrom_size_mb`: Chromosome length in megabases (from GRCh37)  
- Rounded to six decimal places for clarity

---

## Output Files

**Normalized SV Table:**  
`AEsparza_JonesLab_SVComparison_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.tsv`

**Plot Outputs:**  
- PNG: `AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.png`  
- PDF: `AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.pdf`

---

## Visualization Description

Two vertically stacked barplots were generated:

1. **Top Panel:**  
   - Y-axis: `SVs_per_sample_per_Mb`  
   - X-axis: Chromosomes 1–22, X, Y  
   - Grouped by dataset (gnomAD, case_DRAGEN)  
   - Color-coded (`Set2`)  

2. **Bottom Panel:**  
   - Chromosome size (Mb) from GRCh37  
   - Contextual reference for normalization

**Standard deviation (SD) error bars were not plotted**, as the input data was aggregated per chromosome and did not include replicate variance. SD bars can be incorporated if raw replicate-level data is introduced in future analyses.

---
##  AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.png
![AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.png](https://github.com/austinesparza/JonesLab/blob/AEsparza_WorkingFiles_Tables_Images/AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.png)


## Observations

- **DRAGEN cohort** shows a consistent ~4–6 SVs/sample/Mb across most chromosomes, with chromosome 19 peaking above 6.5.  
- **gnomAD burden** appears markedly lower, typically around ~0.04–0.05 SVs/sample/Mb.

This is expected due to:
- Conservative SV definitions in gnomAD  
- Filtering of low-confidence and rare variants  
- Cohort health status (controls)

**Interpretation should focus on relative chromosomal patterns rather than absolute SV burden across datasets.**

---

## File Structure Summary

```
/results/
├── plots/
│   ├── AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.png
│   └── AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.pdf
└── tables/
    ├── AEsparza_JonesLab_SVComparison_NoChrM_CaseGnomAD_2025jul30_v.01.tsv
    └── AEsparza_JonesLab_SVComparison_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.tsv
```

---

## Script Location

```bash
/scripts/AEsparza_JonesLab_SVComparison_NormalizedPlot_2025jul30_v.01.py
```

---

## Future Directions

- Recompute normalized rates with callable base-pair coverage instead of chromosome length  
- Add replicate-level data to compute error bars  
- Stratify by SV type (DEL, DUP, INV, INS) and size classes  
- Extend to additional cohorts for multi-arm comparison

## Objective

Assess the unresolved structural variant burden (BNDs) across chromosomes, comparing:

- **case_DRAGEN**: Case cohort (n = 1016) with SVs called via DRAGEN  
- **control_DRAGEN**: Control cohort (n = 2,945 ) with SVs called via DRAGEN

BNDs (breakends) represent imprecisely resolved structural variants lacking definitive breakpoint mapping. While often excluded from clinical interpretation, they may reflect technical artifacts or alignment instability across genomic regions. This figure evaluates chromosomal BND distribution in both groups to flag potential bias or quality control concerns.

---

## Input Files

**BND Summary Table:**  
`AEsparza_JonesLab_ChromSVTypeSummary_2025jul17_v.01.tsv`  
- Fields used: `chrom`, `group`, `SVTYPE`, `count`

> All chromosomes were included (chr1–22, chrX, chrY, chrM). Only rows where `SVTYPE == BND` were used in this analysis.

---

## Plot Output

**BND Burden Barplot:**  
`AEsparza_JonesLab_BNDChromBurdenPlot_2025jul30_v.01.png`  
- Y-axis: BND count  
- X-axis: Chromosome (ordered 1–22, X, Y, M)  
- Bars grouped by cohort (case, control)  
- Color scheme: seaborn default, with legend for group

##  AEsparza_JonesLab_BND_ChromBurden_CaseControl_2025jul30_v.01.png
![AEsparza_JonesLab_BND_ChromBurden_CaseControl_2025jul30_v.01.png](https://github.com/austinesparza/JonesLab/blob/AEsparza_WorkingFiles_Tables_Images/AEsparza_JonesLab_BND_ChromBurden_CaseControl_2025jul30_v.01.png)
---

## Visualization Description

Grouped barplot showing raw BND counts per chromosome:

- Prominent BND inflation in chr1, chr2, chr5, chr13, and chr20.
- Mitochondrial (chrM) and sex chromosomes (chrY) show negligible counts.
- No normalization was applied, as this was a raw burden check for artifact detection.

---

## Observations

| Finding | Interpretation |
|--------|----------------|
| **Controls exhibit 2–4× higher BND counts** 
| **Specific spike regions (chr1, chr2, chr5, chr13)** | May correspond to segmental duplications, centromeric gaps, or other regions prone to mapping ambiguity |
| **Cases show consistently lower BND rates** | DRAGEN pipeline may be more conservative in reporting unresolved SVs |

BNDs are typically **uninformative for biological association testing** due to lack of breakpoint resolution. This analysis supports **excluding BNDs from downstream burden analyses** or handling them separately as potential pipeline artifacts.

---

## File Structure Summary

```
/results/
├── plots/
│   └── AEsparza_JonesLab_BNDChromBurdenPlot_2025jul30_v.01.png
└── tables/
    └── AEsparza_JonesLab_ChromSVTypeSummary_2025jul17_v.01.tsv
```

---

## Script Location

```bash
/scripts/AEsparza_JonesLab_BNDChromBurdenPlot_2025jul30_v.01.py
```

---

## Next Steps

- Calculate BND **proportion per chromosome**: `BND count / total SVs`  
- Repeat analysis using **DEL, DUP, INS, INV** for comparative context  
- Consider masking BNDs in chromosomal SV burden models  
- Reassess BND distribution after filtering for high-confidence regions or mappability tracks

---

| **Purpose**                                           | **Expected Filename**                                                    | **Intended Location**               |
| ----------------------------------------------------- | ------------------------------------------------------------------------ | ----------------------------------- |
| Chromosome-level SV burden bar plot (all types)       | `AEsparza_JonesLab_ChromSVBarplot_case_2025jul28_v.01.png`               | `results/plots/`                    |
| Chromosome-level SV burden bar plot (all types)       | `AEsparza_JonesLab_ChromSVBarplot_control_2025jul28_v.01.png`            | `results/plots/`                    |
| BND-only SV bar plot (case)                           | `AEsparza_JonesLab_BND_Barplot_case_2025jul28_v.01.png`                  | `results/plots/`                    |
| BND-only SV bar plot (control)                        | `AEsparza_JonesLab_BND_Barplot_control_2025jul28_v.01.png`               | `results/plots/`                    |
| SV type-specific burden bar plots (BND, DEL, DUP, INV)| `AEsparza_JonesLab_<SVTYPE>_BurdenByChrom_2025jul28_v.01.png`            | `results/plots/svtype_comparison/`  |
| Cleaned SV data file (cases)                          | `AEsparza_JonesLab_CleanedCasesSV_2025jul17_v.03.tsv`                    | `data_processed/`                   |
| Cleaned SV data file (controls)                       | `AEsparza_JonesLab_CleanedControlsSV_2025jul17_v.02.tsv`                 | `data_processed/`                   |
| Kolmogorov-Smirnov test results (DEL, DUP)            | `AEsparza_JonesLab_KSTest_SVLen_CaseControl_2025jul28_v.01.tsv`          | `results/stats/`                    |
| Wilcoxon rank-sum test results (DEL, DUP)             | `AEsparza_JonesLab_WilcoxonTest_SVLen_CaseControl_2025jul28_v.01.tsv`    | `results/stats/`                    |
| Combined SV master table (case + control, log-scaled) | `AEsparza_JonesLab_AllSVs_Cleaned_2025jul18_v.01.tsv`                    | `results/dataframes/`               |
| Statistical summary of SV length comparisons          | `AEsparza_JonesLab_SVLengthStats_StatTests_2025jul28_v.01.tsv`           | `results/tables/`                   |
| Updated README with new figures and analysis summary  | `2025jul28_README_DRAGEN_AEsparza_JonesLab_SVCaseControl_GRCh37_v.04.md` | `./` (project root)                 |


## Document History

| Version | Date       | Updates Summary                                                                                         |
|---------|------------|---------------------------------------------------------------------------------------------------------|
| v.01    | 2025-07-17 | Initial documentation of parsing workflow and objectives                                                |
| v.02    | 2025-07-18 | Added merged master table                                                                               |
| v.03    | 2025-07-22 | Split by group/chromosome, added SV type summaries, log-scaled lengths, etc.                            |
| v.04    | 2025-07-28 | Generated SV burden bar plots (case/control, by SV type); performed KS and Wilcoxon tests on SV lengths |
| v.05    | 2025-07-29 | SV burden bar plots (case/control, by SV type) uploaded to JonesLab repo, pathways in report            |



**Current Version:** v.05  
**Maintainer:** Austin Esparza  
