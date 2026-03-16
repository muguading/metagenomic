# Codebase Structure

**Analysis Date:** 2026-03-16

## Directory Layout

```
/Users/wuhhh/Desktop/徐老师/代码/metagenomic/
├── Bac_assemble_260112_newformat.py    # Main pipeline script (4,249 lines)
├── .planning/                          # Planning documentation
│   └── codebase/                       # Codebase analysis documents
│       ├── ARCHITECTURE.md             # Architecture overview
│       └── STRUCTURE.md                # This file
└── .git/                               # Git repository data
```

## File Structure

**Main Script:** `Bac_assemble_260112_newformat.py`
- **Size:** ~277 KB, 4,249 lines
- **Language:** Python 3
- **Encoding:** UTF-8 with Chinese comments
- **Shebang:** `#!/home/dell/miniconda3/envs/TB_ONT/bin/python`

## Script Organization

### Section Breakdown

| Lines | Section | Description |
|-------|---------|-------------|
| 1-140 | Header & Configuration | Imports, argparse, global variables |
| 140-200 | Utility Functions | File validation, logging helpers |
| 200-540 | QC Module | Quality control functions and dataclasses |
| 540-750 | QC Function | Main QC processing (`QC_func()`) |
| 750-785 | Polish Function | Assembly polishing (`polish_func()`) |
| 785-966 | Binning Functions | MAG binning for metagenomics |
| 966-1104 | Assembly Functions | De novo assembly logic |
| 1104-1384 | Reference Assembly | Reference-guided assembly |
| 1384-1738 | Main Assembly | Primary assembly orchestration |
| 1738-1957 | Annotation | Gene prediction and annotation |
| 1957-2040 | Drug/VF Functions | AMR and virulence detection |
| 2040-2231 | Species ID | Kraken2/Bracken classification |
| 2231-2415 | RGI & Helpers | Resistance gene identification |
| 2415-3100 | Serotyping | Species-specific typing functions |
| 3100-3237 | MLST Functions | Multi-locus sequence typing |
| 3237-3314 | Element Prediction | Mobile genetic elements |
| 3314-3478 | Result Combination | Output organization |
| 3478-3876 | Pipeline Functions | Main execution pipelines |
| 3876-3920 | Input Checking | File type detection |
| 3920-3976 | Basecalling | FAST5/POD5 conversion |
| 3976-4249 | Main Execution | Entry point and sample processing |

## Key File Locations

**Entry Point:**
- Line 3978: Main execution block begins
- Line 3845: `main_process()` function definition

**Configuration:**
- Lines 34-57: Argument parser setup
- Lines 72-140: Global variable initialization

**Core Logic:**
- Line 541: `QC_func()` - Quality control
- Line 967: `denovo_asb()` - De novo assembly
- Line 1133: `reassm_fun()` - Reference assembly
- Line 1385: `asb_func()` - Main assembly dispatcher
- Line 1743: `Annotate_func()` - Gene annotation
- Line 2041: `kk2()` - Species identification

**Testing:**
- Line 3583: `run_fake_pipeline()` - Test/demo mode

## Naming Conventions

**Functions:**
- snake_case for all functions
- Chinese comments describe purpose
- Verb-noun pattern: `QC_func()`, `Annotate_func()`, `combine_func()`

**Variables:**
- Mixed English and Chinese variable names
- Global constants: UPPERCASE (e.g., `Krdb`)
- Local variables: lowercase with underscores
- Common abbreviations: `Pre` (prefix), `nt` (num_threads), `ofn` (output filename)

**Files:**
- Script: `Bac_assemble_260112_newformat.py` (date-stamped)
- Outputs: `{Pre}.{extension}` pattern (e.g., `Sample1.final.fasta`)
- Log files: `{step}.log` pattern (e.g., `QC.log`, `asb.log`)

**Classes:**
- PascalCase for dataclasses
- Example: `QCTarget` (line 324)

## Where to Add New Code

**New Analysis Module:**
1. Define function in appropriate section (follow line organization)
2. Add to flow dispatch in `run_fastq_flow()` (line 3705) or `run_fasta_flow()` (line 3794)
3. Add command-line argument in argparse section (line 34)
4. Update documentation in help text

**New Serotyping Function:**
- Location: Lines 2415-3100 (serotyping section)
- Pattern: Follow existing serotype functions (e.g., `serotype_A()`, `serotype_B()`)
- Integration: Add to `mlst_serotype()` dispatcher (lines 3189-3235)

**New QC Step:**
- Location: Within `QC_func()` (lines 541-749)
- Pattern: Use subprocess.run() with shell=True
- Logging: Open dedicated log file

**New External Tool:**
- Requirements: Add to conda environment
- Wrapper: Create function following pattern at lines 162-163
- Call: Integrate into appropriate pipeline stage

## Special Directories

**Runtime Directories (Created During Execution):**

```
{output_dir}/
├── fastq_analysis/              # Main working directory
│   ├── Samplelist.txt           # Sample tracking
│   ├── {Sample}/                # Per-sample directory
│   │   ├── {Sample}.raw.fastq   # Input reads
│   │   ├── {Sample}.final.fastq # QC'd reads
│   │   ├── {Sample}.final.fasta # Assembly
│   │   ├── {Sample}_prokka/     # Gene annotation
│   │   ├── flye_output/         # Assembly info
│   │   ├── {Sample}_genome_complete_result/  # Organized results
│   │   └── *.log                # Step-specific logs
│   └── basecaller_outputs/      # Basecalling results (if applicable)
└── sample_result.txt            # Aggregate results
```

**Result Directory Structure:**

```
{Sample}_genome_complete_result/
├── 1.DataSum/                   # QC summaries
├── 2.Spereads/                  # Species identification
├── 3.Assemble/                  # Assembly files
├── 4.Repeat/                    # Repeat analysis
├── 5.Fun_Element/               # Functional elements
└── 7.Mlst/                      # MLST results
```

## Input/Output Patterns

**Supported Input Types:**
- `fqdir`: Directory with FASTQ files
- `fqfile`: Single FASTQ file
- `f5dir`: Directory with FAST5 files
- `pod5`: Directory with POD5 files
- `bardir`: Barcode-demultiplexed directory
- `fadir`: Directory with FASTA files
- `fafile`: Single FASTA file

**Output File Patterns:**
- Assembly: `{Pre}.final.fasta`
- QC: `{Pre}.QC.summary.tsv`
- Annotation: `{Pre}.prokka.tsv`
- AMR: `{Pre}.card.tsv`, `{Pre}.rgi.tsv`
- Virulence: `{Pre}.vfdb.tsv`
- MLST: `{Pre}.mlst_Stat.txt`
- Serotype: `{Pre}_serotype_result.tsv`

## Dependencies

**External Tools (40+):**
- Assembly: spades.py, flye, unicycler, canu, megahit
- QC: fastp, seqkit, fastqc, porechop, rasusa
- Classification: kraken2, bracken
- Annotation: prokka, bakta, abricate, rgi
- Alignment: minimap2, bwa, samtools
- Variant: freebayes, clair3, bcftools
- Typing: mlst, sistr, ectyper, kleborate
- Visualization: krona, cgview

**Python Libraries:**
- pandas: Data manipulation
- numpy: Numerical operations
- Bio (Biopython): Sequence parsing
- pytaxonkit: Taxonomic operations

**Conda Environments:**
- TB_ONT: Main environment
- medaka: Assembly polishing
- clair3: Variant calling
- gtdbtk: Taxonomic classification
- BASALT: Metagenomic binning
- RGI_new: Resistance gene identification
- And 10+ more specialized environments

## Configuration Files

**Hardcoded Paths:**
- `/data1/shanghai_pip/meta_genome/`: Main database directory
- `/home/dell/miniconda3/`: Conda installation
- `/data/deploy/`: Secondary tools and databases
- `/home/dell/biosoft/`: Custom software installations

**Database References:**
- Kraken2: `/home/dell/kraken2_custom_202101_24G`
- CheckM2: `/data1/shanghai_pip/meta_genome/uniref100.KO.1.dmnd`
- CARD, VFDB: Various locations under `/data/deploy/`

---

*Structure analysis: 2026-03-16*
