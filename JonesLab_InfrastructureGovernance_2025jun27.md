# Jones Lab Infrastructure Standards and Workflow Governance

This document outlines the standards, directory structure, file naming conventions, data hygiene protocols, and audit mechanisms established for the Jones Lab project. It is intended to ensure reproducibility, clarity, and professional-grade project management across all workflows.

---

## File Structure

```
JonesLab/
├── scripts/                      # Executable scripts for audit, filtering, sanitization
├── data_raw/                    # Untouched input files (VCFs, BEDs, etc.)
│   ├── clinvar/
│   ├── BED/
│   └── dbsnp/
├── data_processed/             # Intermediate filtered, cleaned, or annotated data
│   ├── clinvar/
│   └── dbsnp/
├── results/                    # Final output files and filtered results
│   ├── coordmatch/
│   └── isec_coordfilter/
├── logs/                       # Logs of script runs, audits, and terminal sessions
├── docs/                       # Reference literature, SOPs, schematics, and planning
├── BED/                        # Supplemental or backup BED files (deprecated)
└── isec_output/                # bcftools isec outputs (used for variant intersection)
```

---

## File Naming Conventions

All files must include:

* Project identifier: `AEsparza_JonesLab`
* Clear descriptor of workflow/module
* ISO-style timestamp: `YYYYmonDD`
* Version: `v.##`

**Example:**

```
AEsparza_JonesLab_FileAudit_2025jun27_v.01.tsv
AEsparza_JonesLab_FileQuarantine_2025jun27_v.01.py
```

---

## Data Hygiene Protocols

### File Audits

* Run daily or post-task using `AEsparza_JonesLab_FileAudit_YYYYmonDD_v.##.py`
* Generates: `FileAudit_*.tsv` containing path, size, checksum, modification date, and status tag (`raw`, `processed`, `results`, `script`, `documentation`, `other`)

### Zero-byte Quarantine

* Executed via: `AEsparza_JonesLab_FileQuarantine_YYYYmonDD_v.##.py`
* Moves all `0 byte` files to `JonesLab/quarantine/`
* Log file saved to `logs/zero_byte_quarantine_log_YYYYmonDD.tsv`

### Duplicate File Detection

* Based on MD5 checksum match across file system
* Output files: `duplicate_md5_files.tsv` and `true_duplicate_files.tsv`
* Used to avoid redundancy and ensure intermediate file uniqueness

---

## Workflow Commands

### Aliases

To streamline frequent tasks:

```bash
alias JonesLab_Sanitize='python3 ~/Downloads/JonesLab/scripts/AEsparza_JonesLab_FileQuarantine_2025jun27_v.01.py'
```

Add to `~/.zshrc` or `~/.bashrc` and reload shell with `source ~/.zshrc`

### Audit Run

```bash
cd ~/Downloads/JonesLab/scripts
python3 AEsparza_JonesLab_FileAudit_2025jun27_v.01.py
```

### Post-Audit Filtering

```bash
python3 AEsparza_JonesLab_FileAudit_2025jun27_v.01b.py
```

Outputs:

* `zero_byte_files.tsv`
* `duplicate_md5_files.tsv`

### Sanitization Run

```bash
JonesLab_Sanitize
```

---

## Planning and Logging Standards

* Daily logs written to: `logs/terminal_log_YYYYmonDD.txt`
* All terminal sessions must be captured to preserve reproducibility
* Write-protect finalized outputs unless reprocessing is intentional

---

## Notes

* All code and data changes must be versioned.
* No manual file deletions without a logged rationale.
* Every result must trace back to raw data and a script of origin.

---

Maintained by: **Austin Esparza**
Updated: **2025-06-27**
