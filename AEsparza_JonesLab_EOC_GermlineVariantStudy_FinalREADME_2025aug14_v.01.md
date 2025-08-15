# Discovery of Rare Germline Variants in Ovarian Cancer  
**Jones Lab | CS Internship Codebase and Project Repository**  
**Research Project in Germline Variant Annotation and Structural Variant Interpretation in EOC Risk Genes**  
**Author:** Austin Esparza | Summer 2025 UGROW Scholar  

---

## Project Summary

This repository supports a research training project focused on identifying rare germline variants in epithelial ovarian cancer (EOC) susceptibility genes. The project emphasizes reproducible workflows for filtering, annotating, and harmonizing pathogenic variants, particularly SNPs, indels, and structural variants (SVs), using GRCh37-based reference data.


**Current Phase Status (Summer 2025):**  
- **Phase I – SNP/Indel Harmonization and Array-Based Pathogenicity Profiling:**  
  Completed final harmonization and array intersection across a defined panel of EOC risk genes using GRCh37 coordinates.  
  Implemented a dbSNP–ClinVar join strategy to annotate array sites and classify pathogenicity.  
  In the European ancestry subset of the case cohort (_n_ = 155), **41 individuals (~26%)** carried at least one pathogenic variant per ClinVar.  

- **Phase II – Structural Variant (SV) Burden and Distribution Analysis:**  
  Completed per-chromosome SV burden quantification and normalization (SVs per sample per Mb) across cases and controls;  
  ongoing regional enrichment analysis and SV length distribution comparison.  


The project has two integrated goals:

1. **Harmonize pathogenic SNPs and indels** from public databases (ClinVar, dbSNP) within defined ovarian cancer risk gene regions, and assess which of these are represented on the **Infinium Global Diversity Array (GDA)**.
2. **Quantify and annotate structural variants (SVs)** from case and control cohorts (both DRAGEN-called) to characterize per-chromosome burden and identify overlaps with coding and regulatory regions of known EOC risk genes.

---

## Learning Objectives

- Understand the discovery and interpretation of genetic risk variants, from familial studies to GWAS.
- Gain exposure to sequencing platform modalities (array, short-read WGS, long-read WGS).
- Develop bioinformatics proficiency in `bcftools`, `bedtools`, `R`, and command-line tools.
- Learn how pathogenicity is defined across variant types and platforms.
- Investigate the clinical challenges of EOC, including the problem of "missing heritability."
- Navigate public variant databases (e.g., UCSC Genome Browser, gnomAD SVs) to identify disease-relevant features.

---

## Project Tasks

### Part 1: SNP/Indel Harmonization and Genotyping Array Intersection
- Extract gene regions of interest (GRCh37 coordinates) using GENCODE annotations.
- Subset ClinVar and dbSNP files to these regions.
- Normalize and match variants by chromosome, position, reference, and alternate allele.
- Retain entries labeled "Pathogenic" in ClinVar.
- Intersect with external genotyping array site lists to identify pathogenic variants captured by arrays.


**Phase I Updates (2025-07-03):**  
- **Join strategy:** Final strategy resolves sensitivity/specificity by joining on **CHROM + POS + rsID**, following evaluation of location-only and strict allele matches. Reference naming inconsistencies (UCSC `chr17` vs RefSeq `NC_000017.10`) and REF/ALT mismatches were handled using `bcftools norm`. Multiple ClinVar records per rsID were collapsed to a single, consensus label.  
- **Confirmed counts:** From dbSNP-mapped array positions, **3,931 total variants** evaluated; **638** classified as **strictly Pathogenic** by ClinVar.  
- **Array-detectable pathogenicity distribution:** **BRCA1/2 = 98.1%** of calls (**2,986/3,045**); **PALB2 = 36** calls (top non‑BRCA gene).  
- **Carrier prevalence in genotyped cases:** In two EOC case sets, **21–26%** of cases carried a known pathogenic variant.  
- **Outputs:** Curated array-detectable pathogenic set; per-gene pathogenic counts.

### Part 2: Structural Variant (SV) Annotation
- Annotate deletions, insertions, breakends, and duplications overlapping target gene regions.
- Identify which exons, UTRs, or intronic regions are affected.
- Create tables summarizing gene-region impacts for downstream variant interpretation efforts.


**Phase II Updates (2025-07-31):**  
- **Cohorts and callers:** Both case and control cohorts are from DRAGEN-called SV data (GRCh37). SV types include **DEL, DUP, INV, BND**.  
- **Normalization:**  
  **SVs per sample per Mb** = Total SVs ÷ (Sample count × Chromosome size in Mb)  
- **Observed patterns:** Genome-wide **SV burden per Mb is similar** between cases and controls; **deletions cover more genome** than duplications in both cohorts. 
- **Outputs:** Normalized SV burden tables and figures per chromosome and SV type; candidate regions list for follow-up.


---

## Results (Summer 2025)

**Phase I – Array Pathogenicity Profiling:**  
- Produced a **curated set of 638 strictly Pathogenic** sites from dbSNP-mapped array positions within EOC risk genes (GRCh37).  
- **BRCA1/2 dominate array-detectable pathogenic calls (98.1%, 2,986/3,045)**; **PALB2 (n = 36)** leads non‑BRCA genes.  
- In two genotyped EOC case sets, **21–26%** of cases carried a known pathogenic variant.  
- Delivered **carrier classification tables** and **per-gene pathogenic counts**

**Phase II – SV Burden and Distribution:**  
- Computed **SVs per sample per Mb** across chromosomes for cases (DRAGEN) and controls (gnomAD).  
- **Deletion coverage exceeds duplication coverage** across cohorts; **overall burden per Mb is similar** between cases and controls.  
- Identified **chr1, chr2, chr4, chr6, chr7** as top-burden chromosomes in cases; regional enrichment testing is ongoing.  
- Delivered **normalized burden tables and figures** plus a preliminary **candidate regions** list for follow-up.


---

## Genes of Interest (GRCh37 Coordinates)

These genes are involved in homologous recombination, DNA damage repair, or related cancer risk pathways:

| Gene   | Pathway                    | GRCh37 Coordinates                 |
|--------|-----------------------------|------------------------------------|
| BRCA1  | Homologous Recombination   | chr17:41196311-41322262            |
| BRCA2  | Homologous Recombination   | chr13:32890597-32972907            |
| BRIP1  | Fanconi Anemia / HRR       | chr17:59760656-59938900            |
| PALB2  | Homologous Recombination   | chr16:23614779-23652478            |
| RAD51C | HRR / Fanconi Anemia       | chr17:56770004-56811583            |
| RAD51D | Homologous Recombination   | chr17:33338986-33448312            |
| ATM    | DNA Damage Response        | chr11:108179697-108332286          |
| CHEK2  | DNA Damage Checkpoint      | chr22:29117587-29130709            |

---

## Working Directory Tree

```
JonesLab/
├── data_raw/                    # Untouched input files (VCFs, BEDs, etc.)
│   ├── clinvar/
│   ├── dbsnp/
│   └── BED/
├── data_processed/             # Filtered, cleaned, or annotated data
│   ├── clinvar/
│   └── dbsnp/
├── results/                    # Final outputs and matched variants
│   ├── coordmatch/
│   └── isec_coordfilter/
├── scripts/                    # Executable scripts (audit, filters, SV tools)
├── logs/                       # Terminal logs, audit reports, quarantine logs
├── docs/                       # Literature notes, SOPs, schematics
├── isec_output/                # bcftools isec intersection files
├── notebooks/                  # (Optional) Jupyter notebooks or RMarkdowns
├── quarantine/                 # Auto-collected zero-byte files for review
└── FileStructuresAndWorkflowGovernance.md  # Governance standards

```



---

## Tools and Resources

**Reference Genome:** GRCh37 (hg19)  
**Annotation Source:** GENCODE v19 (longest transcript per gene)  
**Key Tools:** `bcftools`, `bedtools`, `R`, UCSC Table Browser  

**Primary Databases**:  
- [dbSNP FTP VCFs](https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/)  
- [ClinVar Documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)  
- [gnomAD Structural Variants](https://gnomad.broadinstitute.org/downloads#svs)

---

## Educational References

- Mandiracioglu B, Ozden F, Kaynar G, Yilmaz MA, Alkan C, Cicek AE. (2024): *ECOLE: Learning to call copy number variants on whole exome sequencing data*. Nature Communications, 15(1). [https://doi.org/10.1038/s41467-023-44116-y](https://doi.org/10.1038/s41467-023-44116-y)

- Collins RL, Brand H, Karczewski KJ, et al. (2020): *A structural variation reference for medical and population genetics*. Nature, 581(7809), 444–451. [https://doi.org/10.1038/s41586-020-2287-8](https://doi.org/10.1038/s41586-020-2287-8)

- Dareng EO, Tyrer JP, Barnes DR, et al. (2022): *Polygenic risk modeling for prediction of epithelial ovarian cancer risk*. European Journal of Human Genetics, 30(3), 349–362. [https://doi.org/10.1038/s41431-021-00987-7](https://doi.org/10.1038/s41431-021-00987-7)

- De La Vega FM, Irvine SA, Anur P, et al. (2024): *Benchmarking of germline copy number variant callers from whole genome sequencing data for clinical applications*. Bioinformatics Advances, 5(1), vbaf071. [https://doi.org/10.1093/bioadv/vbaf071](https://doi.org/10.1093/bioadv/vbaf071)

- DeVries AA, Dennis J, Tyrer JP, et al. (2022): *Copy number variants are ovarian cancer risk alleles at known and novel risk loci*. Journal of the National Cancer Institute, 114(11), 1533–1544. [https://doi.org/10.1093/jnci/djac160](https://doi.org/10.1093/jnci/djac160)

- Jones MR, Kamara D, Karlan BY, Pharoah PDP, Gayther SA. (2017): *Genetic epidemiology of ovarian cancer and prospects for polygenic risk prediction*. Gynecologic Oncology, 147(3), 705–713. [https://doi.org/10.1016/j.ygyno.2017.10.001](https://doi.org/10.1016/j.ygyno.2017.10.001)

- Jones MR, Peng P-C, Coetzee SG, et al. (2020): *Ovarian cancer risk variants are enriched in histotype-specific enhancers and disrupt transcription factor binding sites*. American Journal of Human Genetics, 107(4), 622–635. [https://doi.org/10.1016/j.ajhg.2020.08.021](https://doi.org/10.1016/j.ajhg.2020.08.021)

- Koboldt DC. (2020): *Best practices for variant calling in clinical sequencing*. Genome Medicine, 12(1). [https://doi.org/10.1186/s13073-020-00791-w](https://doi.org/10.1186/s13073-020-00791-w)

- Rosenthal SH, Sun W, Zhang K, et al. (2020): *Development and validation of a 34-gene inherited cancer predisposition panel using next-generation sequencing*. BioMed Research International, 2020, Article ID 3289023. [https://doi.org/10.1155/2020/3289023](https://doi.org/10.1155/2020/3289023)

- Quinlan Lab: *Unix and bioinformatics teaching resources*. University of Utah. [http://quinlanlab.org/teaching.html](http://quinlanlab.org/teaching.html)

---

## Technical Notes and Guides

- Illumina: *Infinium HD Assay Ultra Protocol Guide*. Overview of genotyping workflow for Illumina SNP arrays.  
- Illumina: *iSelect Custom Genotyping Array Design Guide*. Technical specifications for SNP content and probe orientation.  
- ClinGen: *SOP v3.2 for germline variant classification using ACMG/AMP criteria*.  
---

## Final Outputs (Completed to Date and Planned)

- **Curated pathogenic variant table by gene region (Phase I)** — *completed*  
- **Array-detectable pathogenic variant match report** — *completed*  
- **BRCA(+)/BRCA(–) carrier classification tables** — *completed*  
- **Normalized SV burden plots and tables (Phase II)** — *completed*  
- **Annotated SV tables for the gene panel** — *completed*  
- **Versioned and documented bash/R/Python scripts for all steps** — *completed*  

---

## Archive Contents & Provenance

This section documents the curated set of files and directories archived with this repository snapshot (`2025aug14_v.01`), including their exact role in the project and phase association.  
All listed items are versioned and traceable to workflow steps described in **Project Tasks** and **Results**.

---

### 1. Documentation & Governance

| Path | Description & Purpose |
|------|-----------------------|
| `README.md` / `AEsparza_JonesLab_README_2025aug14_v.01.md` | Core project documentation defining scope, methodology, and outputs for both SNP/indel and SV phases, ensuring reproducibility and context for all downstream users. |
| `FileStructuresAndWorkflowGovernance.md` | Governance record describing directory standards, naming conventions, and audit policies to keep the project reproducible and auditable. |
| `docs/JonesLab_InfrastructureGovernance_2025jun27.md` | Extended governance notes for the Jones Lab environment, outlining system-wide reproducibility and QA standards. |

---

### 2. Reference Data

| Path | Description & Purpose |
|------|-----------------------|
| `data_raw/BED/AEsparza_JonesLab_EOCGenes_LongestIsoforms_GRCh37_2025jun30_v.01.bed` | Definitive BED file for target EOC gene coordinates (longest isoform, GRCh37) used in all variant extractions. |
| `BED/AEsparza_JonesLab_EOCGenes_GRCh37_Coordinates_2025jun27_v.01.bed` | Early BED version capturing initial coordinate set for extraction benchmarking. |
| `data_raw/dbsnp/GCF_000001405.25.gz` | Raw dbSNP GRCh37 reference VCF with RefSeq contigs, forming the base for all array site intersections. |
| `ClinVar_GRCh37_20250615/clinvar_GRCh37.vcf.gz` & `.tbi` | ClinVar GRCh37 variant set and index used for pathogenicity annotation of dbSNP and array sites. |

---

### 3. Pathogenicity Join Results (Phase I)

| Path | Description & Purpose |
|------|-----------------------|
| `scripts/2025jul10_0900_dbsnpVariant_ClinVar_Annotation_v.01/AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_LeftJoin_ColRenamed_2025jul10_v.01.tsv` | First integrated join table linking dbSNP rsIDs to ClinVar pathogenicity, critical for Phase I counts. |
| `scripts/2025jul10_0900_dbsnpVariant_ClinVar_Annotation_v.01/data_processed/2025jul11_1200_AEsparza_JonesLab_ClinVar_dbSNP_RSMatch_LeftJoin_Collapsed_v.03.tsv` | Collapsed version resolving multiple ClinVar records per rsID to a single consensus call. |

---

### 4. Phase I – Figures

| Path | Description & Purpose |
|------|-----------------------|
| `results/plots/AEsparza_JonesLab_ClinVarPathogenicCounts_2Panel_2025aug12_v.04.png` | Final barplot of pathogenic counts by gene, split BRCA vs non-BRCA, used in poster and presentations. |
| `results/plots/AEsparza_JonesLab_ClinVarPathogenic_BRCA_vs_Others_2025aug12_v.05.[pdf|png]` | Comparative figure highlighting BRCA dominance in array-detectable pathogenic sites. |
| `results/plots/AEsparza_JonesLab_GDA_FinalAnnotated_Pathogenic_BRCA_vs_Others_2025aug12_v.07.[pdf|png]` | Poster-ready visualization of final annotated pathogenic set for the Infinium GDA. |

---

### 5. Phase I – Tables & Validation

| Path | Description & Purpose |
|------|-----------------------|
| `results/dbSNP_Pathogenic_ByCanonicalGene_2025jul25_v.01.tsv` | Per-gene pathogenic site counts from dbSNP/ClinVar harmonization. |
| `results/AEsparza_JonesLab_CLNSIG_Stats_2025jul08_v.01.tsv` | ClinVar CLNSIG classification breakdown for extracted sites. |
| `results/validation/rsid_clinsig_validation_Cases_NM_2025jul24_v.01.txt` | Validation log confirming CLNSIG matches in case cohort. |
| `results/validation/rsid_falsepositive_Controls_NM_2025jul25_v.01.tsv` | Table of apparent false positives in controls for QC tracking. |

---

### 6. Phase II – SV Burden Tables

| Path | Description & Purpose |
|------|-----------------------|
| `SV_Exploration_2025jul15/results/tables/AEsparza_JonesLab_SVComparison_PerSamplePerMb_2025jul29_v.03.tsv` | First complete normalized SV burden table by chromosome |
| `SV_Exploration_2025jul15/results/tables/AEsparza_JonesLab_SVComparison_PerSamplePerMb_2025jul30_v.03.tsv` | Updated table with DRAGEN case cohort added for comparison. |
| `SV_Exploration_2025jul15/results/tables/AEsparza_JonesLab_SVComparison_PerSamplePerMb_Combined_CaseGnomAD_2025jul30_v.01.tsv` | Combined case/control table for side-by-side plotting. |
| `SV_Exploration_2025jul15/results/tables/AEsparza_JonesLab_SVNormalized_gnomAD_2025jul30_v.02.tsv` | gnomAD-only normalized SV rates for baseline comparison. |
| `SV_Exploration_2025jul15/results/tables/AEsparza_JonesLab_TopBottomChromSVs_2025jul17_v.01.tsv` | Ranked list of highest/lowest SV burden chromosomes in early runs. |

---

### 7. Phase II – Figures

| Path | Description & Purpose |
|------|-----------------------|
| `results/plots/AEsparza_JonesLab_SVRate_PerSamplePerMb_2025jul29_v.03.png` | SV burden plot per Mb for gnomAD baseline. |
| `results/plots/AEsparza_JonesLab_SVRate_PerSamplePerMb_CaseGnomAD_2025jul30_v.01.png` | Direct burden comparison of cases vs controls. |
| `results/plots/AEsparza_JonesLab_SVRate_Deletions_PerMb_CaseControl_2025aug12_v.01.png` | Poster figure showing normalized deletion burden by chromosome. |
| `results/plots/AEsparza_JonesLab_SVRate_Duplications_PerMb_CaseControl_2025aug12_v.01.png` | Poster figure showing normalized duplication burden by chromosome. |
| `results/plots/AEsparza_JonesLab_Top5Chrom_SVType_CaseControl_Normalized_2025aug05_v.01.png` | Top-burden chromosome summary by SV type (normalized). |
| `results/plots/AEsparza_JonesLab_Top5Chrom_SVType_Case_Control_ChrSizeNorm_2025aug05_v.02.png` | Alternate visualization normalizing by chromosome size. |
| `results/plots/svtype_comparison/AEsparza_JonesLab_DEL_BurdenByChrom_2025jul28_v.01.png` | DEL-specific burden profile across chromosomes. |
| `results/plots/svtype_comparison/AEsparza_JonesLab_DUP_BurdenByChrom_2025jul28_v.01.png` | DUP-specific burden profile. |
| `results/plots/svtype_comparison/AEsparza_JonesLab_INV_BurdenByChrom_2025jul28_v.01.png` | INV-specific burden profile. |

---

### 8. Scripts

| Path | Description & Purpose |
|------|-----------------------|
| `scripts/AEsparza_JonesLab_ClinVarPathogenicCounts_PerGene_2025aug12_v.04.py` | Generates Phase I per-gene pathogenicity count plots. |
| `SV_Exploration_2025jul15/scripts/AEsparza_JonesLab_SVComparison_NormalizedPlot_2025jul30_v.05.py` | Plots normalized SV burden for Phase II. |
| `SV_Exploration_2025jul15/scripts/AEsparza_JonesLab_SVRate_ByChrom_PerSamplePerMb_2025jul29_v.01.py` | Chromosome-level SV burden plotting. |
| `SV_Exploration_2025jul15/scripts/AEsparza_JonesLab_SVComparisonNormalized_2025jul29_v.01.py` | Combined case/control normalized burden plotter. |
| `SV_Exploration_2025jul15/scripts/compute_gnomad_sv_burden_2025jul30_v.01.py` | Calculates per-sample per-Mb burden from gnomAD VCF. |

---

### 9. Presentation & Poster

| Path | Description & Purpose |
|------|-----------------------|
| `AEsparza_UGROW_Poster_2025aug10.pdf` | Final summer poster including Phase I and Phase II figures. |
| `AEsparza_LabPresentation_2_2025Aug05_.pdf` | End-of-summer lab presentation summarizing interim results. |
| `INPROGRESS_AEsparza_JonesLab_LabMeetingDeck_2025jul21_v.01.pptx` | First lab meeting deck outlining early analyses and plans. |

---

### 10. Logs

| Path | Description & Purpose |
|------|-----------------------|
| `logs/ClinVar_SubsetExtract_20250701_203203.log` | Audit trail for ClinVar subset extraction. |
| `logs/dbSNP_SubsetExtract_20250701_203610.log` | Audit trail for dbSNP subset extraction (part 1). |
| `logs/dbSNP_SubsetExtract_20250701_204949.log` | Audit trail for dbSNP subset extraction (part 2). |
| `logs/validation_2025jun30.txt` | Validation notes from early coordinate matching tests. |

---
# Supplemental Project README – Phase Context, Documented Methods, and File Mapping  
**Jones Lab | U-GROW Summer 2025**  
**Author:** Austin Esparza  

---

## 1. Project Overview & Phase Goals (Documented)  

**Phase I – SNP/Indel Harmonization & Pathogenicity Profiling**  
- Objective: Filter and annotate array-based SNP/indel positions intersecting with EOC risk genes to identify known pathogenic variants.  
- Rationale (as stated in poster): This approach allows **carrier classification** into BRCA(+) and BRCA(–) groups from array data and defines the subset of variants most clinically interpretable.  

**Phase II – Structural Variant Burden and Distribution**  
- Objective: Quantify per-chromosome burden of deletions, duplications, inversions, and breakends (BND) in EOC cases and population controls.  
- Rationale (as stated in presentations): Understanding **SV distribution and chromosomal burden** may identify regions contributing to ovarian cancer risk beyond BRCA1/2-associated pathways.  

---

## 2. Documented Methods & Analysis Steps  

**From Poster & Presentations**  

**Phase I:**  
1. Input: Infinium Global Diversity Array site list, ClinVar VCF, dbSNP VCF.  
2. Subset ClinVar and dbSNP to longest transcript coordinates for 8 risk genes.  
3. Normalize coordinates and alleles to resolve contig naming differences.  
4. Join on **CHROM + POS + rsID** to map pathogenic ClinVar entries to array positions.  
5. Output: Per-gene pathogenic variant counts; carrier classification tables.  

**Phase II:**  
1. Input: DRAGEN-called SVs for ~1,016 EOC cases; ~2,945 controls
2. SV types: DEL, DUP, INV, BND.  
3. Chromosome size normalization applied: SVs per sample per megabase.  
4. Compare per-chromosome burden between cases and controls.  
5. Output: Burden plots (normalized and raw); top-burden chromosome identification; early candidate list for follow-up.  

---

## 3. Documented Results & Figures (Source Mapping)  

| Figure/File | Documented Source | Summary from Poster/Presentation |
|-------------|------------------|-----------------------------------|
| `AEsparza_JonesLab_ClinVarPathogenicCounts_2Panel_2025aug12_v.04.png` | Poster Fig. 4 | Per-gene pathogenic counts from Phase I; BRCA vs non-BRCA split. |
| `AEsparza_JonesLab_GDA_FinalAnnotated_Pathogenic_BRCA_vs_Others_2025aug12_v.07.png` | Poster Fig. 6 | Final annotated pathogenic set on GDA; visualization not normalized for gene size. |
| `AEsparza_JonesLab_SVRate_Deletions_PerMb_CaseControl_2025aug12_v.01.png` | Poster Fig. 7 | Normalized deletion burden by chromosome for Phase II. |
| `AEsparza_JonesLab_SVRate_Duplications_PerMb_CaseControl_2025aug12_v.01.png` | Poster Fig. 8 | Normalized duplication burden by chromosome for Phase II. |

---

## 4. File-to-Phase Mapping  

| Path | Phase | Role in Analysis |
|------|-------|------------------|
| `scripts/2025jul10_0900_dbsnpVariant_ClinVar_Annotation_v.01/...` | Phase I | Joins dbSNP and ClinVar pathogenic entries for array site annotation. |
| `results/dbSNP_Pathogenic_ByCanonicalGene_2025jul25_v.01.tsv` | Phase I | Per-gene pathogenic count table from final join. |
| `SV_Exploration_2025jul15/results/tables/AEsparza_JonesLab_SVComparison_PerSamplePerMb_2025jul30_v.03.tsv` | Phase II | Normalized SV burden table (cases + controls). |
| `results/plots/AEsparza_JonesLab_Top5Chrom_SVType_CaseControl_Normalized_2025aug05_v.01.png` | Phase II | Visualization of top-burden chromosomes per SV type. |

---

## 5. Data Sources & Versions (As Stated in Materials)  

- **Reference Genome:** GRCh37 (hg19).  
- **GENCODE v19:** Used for defining longest transcript coordinates of risk genes.  
- **ClinVar GRCh37 release:** 2025-06-15.  
- **dbSNP GRCh37 build:** GCF_000001405.25.  
- **gnomAD SVs:** Version not stated in poster, but confirmed as GRCh37-based in presentations.  
- **Infinium Global Diversity Array (GDA):** Array platform for SNP/indel coverage mapping.  

---


## Author


**Austin Esparza**  
Research Training Scholar  
[GitHub Profile](https://github.com/austinesparza)  
austin.esparza@cshs.org  
austin.m.esparza@gmail.com  

## Data Privacy and Compliance Statement

This repository includes only publicly available genomic data and reference annotations. No internal cohort data, personally identifiable information (PII), or protected health information (PHI) are included or processed. All analyses are conducted on public datasets such as ClinVar, dbSNP, and GENCODE, using the GRCh37/38 reference genome. This repository adheres to institutional data-sharing policies and maintains full compliance with research ethics and transparency standards.
