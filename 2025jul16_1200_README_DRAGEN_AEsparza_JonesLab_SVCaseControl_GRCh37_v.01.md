# Structural Variant Exploration: Germline Case-Control Comparison (GRCh38)

## Project Overview

This project aims to characterize germline structural variant (SV) patterns using DRAGEN-called data from a case-control cohort, provided by Nimisha. The data are aligned to the GRCh38 reference genome and focus on deletions, duplications, inversions, and unresolved breakends (BNDs).



## Assignment Scope

Nimisha provided two DRAGEN-called SV files representing case and control cohorts. The following objectives were specified:

1. Identify chromosomes with the highest and lowest total SV counts. Determine whether these observations are consistent with previously reported germline variant patterns.
2. Calculate the full length distribution of SVs, including deletions, duplications, and inversions. Compare length trends between case and control groups.
3. Include unresolved variants (breakends, or BNDs) in burden counts. Consider potential reasons for lack of resolution and evaluate chromosome-level clustering.
4. Recommend appropriate visualizations for summarizing the above trends.


# Jones Lab SV Standardization Summary  
**Project:** *Discovery of Rare Germline Structural Variants in Epithelial Ovarian Cancer*  
**Date:** July 17, 2025  
**Author:** Austin Esparza  

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
SV_Exploration_2025jul15
├── 1kgp
│   ├── 1kgp_raw_extract.tsv
│   ├── ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz
│   ├── ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz.tbi
│   ├── ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
│   ├── ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi
│   └── README_phase3_sv_callset_20150224
├── Gnomad_SV_VCF
│   ├── BND_ChrBurden_CaseControl_DRAGEN_2025aug04.png
│   ├── gnomad.v4.1.sv.sites.vcf.gz
│   ├── gnomad.v4.1.sv.sites.working.vcf.gz
│   ├── gnomad_raw_extract.tsv
│   ├── gnomad_raw_extract.with_ids_2025jul30_v.01.tsv
│   └── gnomad_raw_extract.working.tsv
├── clinvar_pathogenic_variants_per_gene.png
├── data_processed
├── data_raw
├── logs
├── results
│   ├── dataframes
│   │   ├── AEsparza_JonesLab_AllSVs_Cleaned_2025jul18_v.01.tsv
│   │   ├── AllSVs_1kgp
│   │   │   ├── SVs_1kgp_chr10_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr11_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr12_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr13_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr14_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr15_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr16_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr17_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr18_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr19_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr1_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr20_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr21_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr22_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr2_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr3_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr4_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr5_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr6_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr7_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr8_2025jul29_v.01.tsv
│   │   │   ├── SVs_1kgp_chr9_2025jul29_v.01.tsv
│   │   │   └── SVs_1kgp_chrX_2025jul29_v.01.tsv
│   │   ├── AllSVs_Case
│   │   ├── AllSVs_Control
│   │   └── AllSVs_gnomad
│   │       ├── gnomad_chr10_sv.tsv
│   │       ├── gnomad_chr10_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr11_sv.tsv
│   │       ├── gnomad_chr11_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr12_sv.tsv
│   │       ├── gnomad_chr12_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr13_sv.tsv
│   │       ├── gnomad_chr13_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr14_sv.tsv
│   │       ├── gnomad_chr14_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr15_sv.tsv
│   │       ├── gnomad_chr15_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr16_sv.tsv
│   │       ├── gnomad_chr16_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr17_sv.tsv
│   │       ├── gnomad_chr17_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr18_sv.tsv
│   │       ├── gnomad_chr18_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr19_sv.tsv
│   │       ├── gnomad_chr19_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr1_sv.tsv
│   │       ├── gnomad_chr1_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr20_sv.tsv
│   │       ├── gnomad_chr20_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr21_sv.tsv
│   │       ├── gnomad_chr21_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr22_sv.tsv
│   │       ├── gnomad_chr22_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr2_sv.tsv
│   │       ├── gnomad_chr2_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr3_sv.tsv
│   │       ├── gnomad_chr3_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr4_sv.tsv
│   │       ├── gnomad_chr4_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr5_sv.tsv
│   │       ├── gnomad_chr5_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr6_sv.tsv
│   │       ├── gnomad_chr6_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr7_sv.tsv
│   │       ├── gnomad_chr7_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr8_sv.tsv
│   │       ├── gnomad_chr8_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chr9_sv.tsv
│   │       ├── gnomad_chr9_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chrX_sv.tsv
│   │       ├── gnomad_chrX_sv_2025jul29_v.02.tsv
│   │       ├── gnomad_chrY_sv.tsv
│   │       └── gnomad_chrY_sv_2025jul29_v.02.tsv

```
**Total:** 17 directories, 229 files


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

## AEsparza_JonesLab_SVLengthStats_AllTypes_2025jul17_v.01
![AEsparza_JonesLab_SVLengthStats_AllTypes_2025jul17_v.01] <img width="812" height="273" alt="image" src="https://github.com/user-attachments/assets/f6127cb1-00f5-4f9f-bd36-62e59dd1a9ca" />

## AEsparza_JonesLab_TopBottomChromSVs_2025jul17_v.01
![AEsparza_JonesLab_TopBottomChromSVs_2025jul17_v.01] 
<img width="220" height="373" alt="image" src="https://github.com/user-attachments/assets/034173ba-565e-45cd-a392-62e21278453b" />


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
- [x] Apply statistical comparisons between case/control distributions (e.g., KS test, Wilcoxon)


## Notes

- All coordinate-based work must remain tied to GRCh38.
- Visualization choices should prioritize interpretability and publication-readiness.
- Breakends should be handled as a distinct category in all summaries.
- Versioned outputs are required at every stage for reproducibility.

## Supplemental Context and Analytical Rationale

This analysis was initiated to investigate structural variant (SV) patterns in a germline case-control dataset generated using the DRAGEN variant caller and aligned to the GRCh38 reference genome. The provided SV callsets, one for cases and one for controls, contain chromosomal coordinates, SV type, and length data for events including deletions, duplications, inversions, insertions, and unresolved breakends (BNDs). The immediate goal is to produce a foundational characterization of SV burden and length distribution in support of downstream variant prioritization and visualization.

The initial analytical scope is defined by the following tasks:

- Identify chromosomes with the highest and lowest counts of SVs, and determine whether these distributions reflect previously reported germline variant patterns.
- Analyze the full length range of deletions and duplications, and assess whether their distributions differ between cases and controls.
- Extend this evaluation to include inversions, and determine whether similar trends are observed across groups.
- Examine the distribution of unresolved breakends (BNDs), treating these as a separate category. This includes quantifying their chromosomal burden and considering possible reasons for partial resolution during variant calling.

The aim of this work is to inform the design of scalable, reproducible SV analysis workflows that are compatible with other GRCh38-based datasets, including genotyping arrays such as the Infinium Global Diversity Array. All coordinate systems and annotations remain GRCh38-specific. Visualization strategies are selected to prioritize clarity, interpretability, and publication-quality output. This document will be expanded as additional input files are incorporated or as analytical priorities shift.

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
- This supports further characterization or stratification by genomic context or clinical phenotype.

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
**Genome Build:** GRCh38  

---

## Objective

Quantify and compare the burden of structural variants (SVs) per chromosome between:

- **case_DRAGEN**: Case samples (n = 1016) processed via DRAGEN-based SV calling
- **control_DRAGEN**: control samples (n = 2,945) processed via DRAGEN-based SV calling
SV counts are normalized by both sample size and chromosome length to yield:  
**SVs per sample per megabase (Mb)** — a metric allowing intra-cohort chromosomal burden comparison.

**Note:** Due to differences in SV detection algorithms, definitions, and quality filters between pipelines, comparisons across datasets are approximate and **should not be interpreted as absolute biological differences.**

---

**Chromosome Sizes:**  
Reference values for GRCh38 autosomes and sex chromosomes (Mb). Stored internally in script.

---

## Normalization Formula

```python
SVs_per_sample_per_Mb = total_SVs / (sample_count × chrom_size_mb)
```

- `total_SVs`: SV count per chromosome  
- `sample_count`: 1016 for DRAGEN cases; 2945 for controls
- `chrom_size_mb`: Chromosome length in megabases (from GRCh38)  
- Rounded to six decimal places for clarity

---

## Output Files

**Normalized SV Table:**  
`AEsparza_JonesLab_SVComparison_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.tsv`

**Plot Outputs:**  
- PNG: `AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.png`  
- PDF: `AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.pdf`

## Source and Generation Details

- **Source file:** `results/tables/AEsparza_JonesLab_SVLengthStats_AllTypes_2025jul17_v.01.tsv`  
- **Generation script:** `scripts/AEsparza_JonesLab_SVLengthStats_Tables_2025jul17_v.01.py`  
- **Description:** This file contains per-cohort structural variant (SV) length statistics for DRAGEN-called case and control datasets, stratified by SV type (BND, DEL, DUP, INS, and INV where present).  
  Statistics include count, minimum, median, mean, maximum, standard deviation, and quartile-based spread (Q1, Q3, IQR).

---

## Structural Variant Length Statistics — DRAGEN Case vs Control

The table below summarizes count, size distribution, and spread for each structural variant (SV) type
in the DRAGEN-called case and control cohorts.

```
| sv_type   | group   |    count |   min |   median |           mean |            max |              std |   Q1 |   Q3 |   IQR |
|:----------|:--------|---------:|------:|---------:|---------------:|---------------:|-----------------:|-----:|-----:|------:|
| BND       | case    |  2076463 |     1 |        1 |    1           |    1           |      0           |    1 |    1 |     0 |
| BND       | control |  6634637 |     1 |        1 |    1           |    1           |      0           |    1 |    1 |     0 |
| DEL       | case    |  4840312 |    50 |      190 | 4553.24        |    9.1229e+07  | 136421           |   76 |  340 |   264 |
| DEL       | control | 14760659 |    50 |      205 | 6389.88        |    9.1229e+07  | 167910           |   77 |  354 |   277 |
| DUP       | case    |    46903 |   998 |     4052 |    1.01797e+06 |    1.2595e+08  |      7.25616e+06 | 1698 | 7880 |  6182 |
| DUP       | control |   181494 |  1000 |     3936 |    1.46502e+06 |    9.95544e+07 |      8.72255e+06 | 1646 | 7381 |  5735 |
| INS       | case    |  4179534 |    50 |      101 |  156.508       | 1836           |    131.843       |  nan |  nan |   nan |
| INS       | control | 12710188 |    50 |      102 |  158.724       | 2079           |    135.921       |  nan |  nan |   nan |
```
## Source and Generation Details — Top and Bottom Chromosomes by SV Count

- **Source file:** `results/tables/AEsparza_JonesLab_TopBottomChromSVs_2025jul17_v.01.tsv`  
- **Generation script:** `scripts/AEsparza_JonesLab_SVSummary_2025jul17_v.01.py`  
- **Description:** This file lists the chromosomes with the highest and lowest total SV counts for each cohort
  (case and control) in the DRAGEN-called dataset.

---

## Chromosomes with Highest and Lowest SV Burden — DRAGEN Case vs Control

```
| group   | chrom   |   total_sv_count |
|:--------|:--------|-----------------:|
| case    | chr1    |          1026503 |
| case    | chr2    |          1014179 |
| case    | chr6    |           840116 |
| case    | chrM    |               28 |
| case    | chrY    |             2718 |
| case    | chr21   |           235695 |
| control | chr1    |          3176535 |
| control | chr2    |          3123131 |
| control | chr6    |          2579239 |
| control | chrM    |               53 |
| control | chrY    |             7935 |
| control | chr21   |           728390 |
```


---

## Visualization Description

Two vertically stacked barplots were generated:

1. **Top Panel:**  
   - Y-axis: `SVs_per_sample_per_Mb`  
   - X-axis: Chromosomes 1–22, X, Y  
   - Grouped by dataset (gnomAD, case_DRAGEN)  
   - Color-coded (`Set2`)  

2. **Bottom Panel:**  
   - Chromosome size (Mb) from GRCh38  
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


##  AEsparza_WorkingFiles_Tables_Images/AEsparza_JonesLab_BND_ChromBurden_Normalized_2025jul30_v.02.png
![AEsparza_WorkingFiles_Tables_Images/AEsparza_JonesLab_BND_ChromBurden_Normalized_2025jul30_v.02.png](https://github.com/austinesparza/JonesLab/blob/AEsparza_WorkingFiles_Tables_Images/AEsparza_JonesLab_BND_ChromBurden_Normalized_2025jul30_v.02.png)

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

- Consider masking BNDs in chromosomal SV burden models  
- Reassess BND distribution after filtering for high-confidence regions

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
^^Updated readme filed with incorrect build name. GRCh38 confirmed used

## Reproducibility Checklist — DRAGEN Case vs Control SV Analysis

This section documents the exact sequence of steps required to reproduce the analysis from raw DRAGEN SV calls to final tables and plots.  
All scripts referenced are located in `scripts/` unless otherwise specified.  
All file paths are relative to the project root: `SV_Exploration_2025jul15/`.

---

### 1. **Prepare Input Data**
**Goal:** Ensure DRAGEN-called case and control SV files are available in `data_raw/`.

**Files:**
- `data_raw/all_cases.dragen.sv.txt`
- `data_raw/all_controls.dragen.sv.txt`

**Notes:**
- These are unmodified outputs from the DRAGEN pipeline.
- Confirm integrity using stored MD5 checksums (`logs/md5_sums.txt`).

---

### 2. **Clean and Standardize SV Data**
**Goal:** Filter, format, and standardize the DRAGEN outputs for downstream parsing.

**Scripts:**
- `AEsparza_JonesLab_ParseCasesSV_2025jul17_v.01.py`
- `AEsparza_JonesLab_ParseControlsSV_2025jul17_v.01.py`

**Outputs:**
- `data_processed/AEsparza_JonesLab_CleanedCasesSV_2025jul17_v.01.tsv`
- `data_processed/AEsparza_JonesLab_CleanedControlsSV_2025jul17_v.01.tsv`

---

### 3. **Generate Summary Counts by Chromosome and SV Type**
**Goal:** Quantify SVs by chromosome and type (DEL, DUP, INV, BND, INS).

**Script:**
- `AEsparza_JonesLab_SVSummary_2025jul17_v.01.py`

**Outputs:**
- `results/tables/AEsparza_JonesLab_ChromSVTypeSummary_2025jul17_v.01.tsv`
- `results/tables/AEsparza_JonesLab_TopBottomChromSVs_2025jul17_v.01.tsv` ← *(Top/Bottom table)*

---

### 4. **Calculate Length Statistics for Each SV Type**
**Goal:** Summarize SV length distributions across case/control cohorts.

**Script:**
- `AEsparza_JonesLab_SVLengthStats_Tables_2025jul17_v.01.py`

**Output:**
- `results/tables/AEsparza_JonesLab_SVLengthStats_AllTypes_2025jul17_v.01.tsv` ← *(Stats table)*

---

### 5. **Produce Length Range Tables**
**Goal:** Capture min/max ranges per SV type, optionally per chromosome.

**Scripts:**
- `AEsparza_JonesLab_SVLengthRange_2025jul17_v.01.py`

**Outputs:**
- `results/tables/AEsparza_JonesLab_SVLengthRange_DEL_DUP_INV_2025jul17_v.01.tsv`
- `results/tables/AEsparza_JonesLab_SVLengthRange_ByChrom_DEL_DUP_INS_2025jul17_v.01.tsv`

---

### 6. **Normalize SV Counts by Sample Size and Chromosome Length**
**Goal:** Enable fair burden comparisons between cohorts.

**Scripts:**
- `AEsparza_JonesLab_SVComparison_NormalizedPlot_2025jul30_v.01.py` *(for gnomAD comparison)*
- `AEsparza_JonesLab_SVComparison_PerSamplePerMb_2025jul30_v.03.py` *(for DRAGEN case vs control)*

**Inputs:**
- `results/tables/AEsparza_JonesLab_SVQuantification_ByChromGroup_2025jul22_v.01.tsv`
- Chromosome size reference: GRCh37 (sourced from UCSC `chromInfo.txt`).

**Outputs:**
- Normalized TSV files in `results/tables/`
- Corresponding burden plots in `results/plots/`

---

### 7. **Generate Figures**
**Goal:** Visualize SV distribution and burden patterns.

**Scripts & Outputs:**
- Violin plots: `results/plots/AEsparza_JonesLab_SVLEN_ViolinPlot_2025jul17_v.01.png`
- Bar plots for normalized burden:  
  - `AEsparza_JonesLab_SVRate_Deletions_PerMb_CaseControl_2025aug12_v.01.png`  
  - `AEsparza_JonesLab_SVRate_Duplications_PerMb_CaseControl_2025aug12_v.01.png`

---

### 8. **Integrate Summary Tables into README**
**Goal:** Document final results with direct links to the generated tables and plots.

**Included Tables:**
- **SV Length Stats Table** → `AEsparza_JonesLab_SVLengthStats_AllTypes_2025jul17_v.01.tsv`
- **Top/Bottom Chromosomes Table** → `AEsparza_JonesLab_TopBottomChromSVs_2025jul17_v.01.tsv`

**Included Figures:**
- Normalized burden plots for deletions and duplications
- SV length violin plots

---

**Final Verification:**
1. All output files present in `results/tables/` and `results/plots/`
2. README cross-references match actual filenames and versions
3. MD5 checksums confirm no corruption since last run
4. All plots regenerate without error when scripts are re-run from raw files

---


## Document History

| Version | Date       | Updates Summary                                                                                         |
|---------|------------|---------------------------------------------------------------------------------------------------------|
| v.01    | 2025-07-17 | Initial documentation of parsing workflow and objectives                                                |
| v.02    | 2025-07-18 | Added merged master table                                                                               |
| v.03    | 2025-07-22 | Split by group/chromosome, added SV type summaries, log-scaled lengths, etc.                            |
| v.04    | 2025-07-28 | Generated SV burden bar plots (case/control, by SV type); performed KS and Wilcoxon tests on SV lengths |
| v.05    | 2025-07-29 | SV burden bar plots (case/control, by SV type) uploaded to JonesLab repo, pathways in report            |
| v.06    | 2025-08-14 | Added comprehensive reproducibility checklist and file provenance details to README                     |


**Current Version:** v.06  
**Maintainer:** Austin Esparza  
