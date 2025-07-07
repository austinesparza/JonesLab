# Jones Lab – ClinVar Allele-Level Annotation Workflow

**Author:** Austin Esparza  
**Last updated:** 2025-07-07  

---

## Objective
Annotate Infinium Global Diversity Array variants with allele-specific ClinVar data (June 2025 release) and produce a compressed, indexed VCF suitable for downstream analyses. Exact matches on `CHROM | POS | REF | ALT` eliminate earlier mis-annotations caused by position-only matching.

---
## Software Versions

| Tool        | Version Tested     |
|-------------|--------------------|
| bcftools    | 1.20               |
| tabix       | 1.20 (via htslib)  |
| awk         | BSD/macOS default  |
| bash        | 5.x+               |


## Input Files

| File                                                         | Purpose                                                                   |
|--------------------------------------------------------------|---------------------------------------------------------------------------|
| `AEsparza_JonesLab_ArrayConvertedVCF_2025jul02_v.01.vcf`     | Array-derived variants (3 931 records)                                    |
| `clinvar_20250630_allele_level_annotations_clean.tsv.gz`     | Allele-level ClinVar table (`CHROM POS REF ALT RS CLNID CLNSIG`)          |
| `header_clean.tmp`                                           | Corrected VCF header (adds INFO tags and contig declarations)             |

---

## Key Problems Addressed
- **Position-only matching** collapsed distinct alleles (e.g., benign `rs799917` vs. pathogenic `rs80357962`).  
- Header contained an invalid placeholder `##INFO=<ID=.,…>` and lacked contig declarations, causing `Invalid tag name "."` and `Contig … not defined` errors.

---

## Step-by-Step Solution

# Extract and sanitize header from array-derived VCF to correct invalid tags and define missing fields.
> **Note:** bcftools annotate requires exact `CHROM, POS, REF, ALT` matches to transfer INFO fields. This eliminates spurious assignments caused by rsID or position-only matches.


### 1. Extract Original Header
```bash
bcftools view -h AEsparza_JonesLab_ArrayConvertedVCF_2025jul02_v.01.vcf \
  > header_clean.tmp       # Extract header

### 2. Remove Invalid INFO Tag
sed -i '' '/ID=\./d' header_clean.tmp    # macOS-compatible in-place edit

### 3. Append INFO Definitions
cat <<EOF >> header_clean.tmp
##INFO=<ID=RS_clinvar_v.20250630,Number=1,Type=String,Description="dbSNP RS ID">
##INFO=<ID=CLNID_clinvar_v.20250630,Number=1,Type=String,Description="ClinVar internal ID">
##INFO=<ID=CLNSIG_clinvar_v.20250630,Number=.,Type=String,Description="Clinical significance from ClinVar">
EOF

### 4. Generate Contig Declarations
bcftools view -H AEsparza_JonesLab_ArrayConvertedVCF_2025jul02_v.01.vcf \
  | cut -f1 | sort -u > contigs.tmp
awk '{print "##contig=<ID="$1">"}' contigs.tmp >> header_clean.tmp

### 5. Prepare ClinVar Allele-Level TSV
# Filter to PASS SNVs and index
bcftools view -Oz -f PASS -V indels \
  -o clinvar_20250630_corefields.vcf.gz \
  clinvar_20250630.vcf.gz
tabix -p vcf clinvar_20250630_corefields.vcf.gz

# Extract allele-specific fields
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%RS\t%ID\t%CLNSIG\t%CLNSIGCONF\n' \
  clinvar_20250630_corefields.vcf.gz \
  > clinvar_20250630_allele_level_annotations.tsv

sort -k1,1 -k2,2n clinvar_20250630_allele_level_annotations.tsv \
  | bgzip -c > clinvar_20250630_allele_level_annotations_clean.tsv.gz
tabix -s1 -b2 -e2 clinvar_20250630_allele_level_annotations_clean.tsv.gz


# Compress and index TSV
bgzip -f clinvar_20250630_allele_level_annotations.tsv
tabix -s1 -b2 -e2 clinvar_20250630_allele_level_annotations.tsv.gz

# Sort and reformat TSV for annotation compatibility
bcftools sort -t $'\t' -k1,1 -k2,2n clinvar_20250630_allele_level_annotations.tsv \
  | bgzip -c > clinvar_20250630_allele_level_annotations_clean.tsv.gz
tabix -s1 -b2 -e2 clinvar_20250630_allele_level_annotations_clean.tsv.gz


### 6. Annotate Array VCF
bcftools annotate \
  -a clinvar_20250630_allele_level_annotations_clean.tsv.gz \
  -c CHROM,POS,REF,ALT,\
INFO/RS_clinvar_v.20250630,\
INFO/CLNID_clinvar_v.20250630,\
INFO/CLNSIG_clinvar_v.20250630 \
  -h header_clean.tmp \
  -o AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf \
  -O v \
  AEsparza_JonesLab_ArrayConvertedVCF_2025jul02_v.01.vcf
# No fatal errors; only harmless legacy “.” warnings.

### 7. Compress and Index Annotated VCF
bgzip -f AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf
tabix -p vcf AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf.gz

### 8. Verify Annotation Count
bcftools view -H AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf.gz \
  | cut -f8 | grep -c "CLNSIG_clinvar_v.20250630"
# Expected output: 68

### Output Files
| File                                                         | Description                              |
| ------------------------------------------------------------ | ---------------------------------------- |
| `AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf.gz`     | Final allele-level ClinVar-annotated VCF |
| `AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf.gz.tbi` | Tabix index                              |

### QC Results
Total variants in array VCF: 3 931

Variants with ClinVar allele-level annotation: 68 (≈ 1.7 %)
– reflects strict allele matching on a population genotyping array.

No header or contig parsing errors remain.

### Proposed Next Actions
Optional filtering (pathogenic-only):
bcftools view -i 'INFO/CLNSIG_clinvar_v.20250630 ~ "Pathogenic"' \
  AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf.gz \
  > pathogenic.vcf

Add dbSNP fields with an additional bcftools annotate call.

Package the workflow into a reproducible script
(AEsparza_JonesLab_ClinVarAnnotate_ArrayVCF_2025jul07_v.01.sh).

### Supplemental Workflow Info

# Initial Array Audit Script

bash scripts/AEsparza_JonesLab_ArrayAudit_2025jul02_v.01.sh
(See audit summary in Initial Array Audit above.)

Pathogenic-Only TSV Extraction

# Export all annotated variants to TSV
bcftools query -H \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%RS_clinvar_v.20250630\t%CLNID_clinvar_v.20250630\t%CLNSIG_clinvar_v.20250630\n' \
  AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf.gz \
  > ArrayAnnotated_AllCLNSIG.tsv

# Filter strict “Pathogenic”
awk -F'\t' 'NR==1 || $8=="Pathogenic"' \
  ArrayAnnotated_AllCLNSIG.tsv \
  > ArrayAnnotated_PathogenicStrict.tsv

Reproducible Pipeline Script (Skeleton)
# AEsparza_JonesLab_ClinVarAnnotate_ArrayVCF_2025jul07_v.01.sh
# 1. Array audit bash scripts/AEsparza_JonesLab_ArrayAudit_2025jul02_v.01.sh
# 2. ClinVar TSV preparation (filter, query, compress, index)
# 3. Header cleanup and contig injection
# 4. bcftools annotate with allele-level ClinVar
# 5. bgzip + tabix
# 6. Export summary TSV + optional pathogenic filter

md5sum AEsparza_JonesLab_ArrayAnnotated_2025jul07_v.13.vcf.gz > final.md5
