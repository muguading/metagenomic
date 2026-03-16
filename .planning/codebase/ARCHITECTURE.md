# Architecture

**Analysis Date:** 2026-03-16

## Pattern Overview

**Overall:** Monolithic Pipeline Script with Modular Functional Architecture

**Key Characteristics:**
- Single-file Python script (~4,200 lines) containing complete bacterial genome assembly pipeline
- Functional programming approach with 80+ specialized functions
- Sequential processing pipeline with branching logic for different data types
- Heavy reliance on external bioinformatics tools via subprocess calls
- Mixed Chinese/English code comments and documentation
- Dual-mode support: command-line tool and importable module

## Layers

**Configuration & CLI Layer:**
- Purpose: Argument parsing, global configuration, and environment setup
- Location: Lines 1-140
- Contains: argparse configuration, global variables, database paths
- Key Components: 25+ command-line arguments for pipeline customization

**Data Input/Validation Layer:**
- Purpose: File type detection, input validation, and format conversion
- Location: Lines 140-200, 3877-3920
- Contains: `is_fasta()`, `is_fastq()`, `check_input()` functions
- Supports: FASTQ, FASTA, FAST5, POD5, barcode directories, cfg files

**Basecalling Layer (Optional):**
- Purpose: Convert raw nanopore signals to FASTQ
- Location: Lines 3921-3975
- Contains: `basecaller()`, `get_free_gpu_memory()` functions
- Tools: Guppy (FAST5), Dorado (POD5)

**Quality Control Layer:**
- Purpose: Read filtering, trimming, and quality assessment
- Location: Lines 541-749
- Contains: `QC_func()`, `ngs_qc()`, helper functions
- Tools: fastp, seqkit, porechop, rasusa, hostile, fastqc

**Assembly Layer:**
- Purpose: Genome assembly from short/long reads
- Location: Lines 966-1104, 1384-1738
- Contains: `denovo_asb()`, `reassm_fun()`, `asb_func()`, `polish_func()`
- Tools: SPAdes, Flye, Unicycler, MaSuRCA, Canu, Raven, wtdbg2, miniasm, meta (Megahit+BASALT)

**Annotation Layer:**
- Purpose: Gene prediction and functional annotation
- Location: Lines 1742-1957, 3237-3314
- Contains: `Annotate_func()`, `AnnoEle()`
- Tools: Prokka, bakta, PhiSpy, minced, IslandPath-DIMOB, mefinder

**Species Identification Layer:**
- Purpose: Taxonomic classification and species identification
- Location: Lines 2040-2231
- Contains: `kk2()`, `proc_kra()`, `proc_kra1()`, `exreadsID()`
- Tools: Kraken2, Bracken

**Serotyping Layer:**
- Purpose: Specialized serotype/strain typing for pathogenic bacteria
- Location: Lines 2415-2976
- Contains: 15+ serotyping functions for specific pathogens
- Tools: MLST, SISTR, ectyper, kleborate, VPsero, emm_typing, etc.

**AMR/Virulence Detection Layer:**
- Purpose: Antimicrobial resistance and virulence gene identification
- Location: Lines 1959-2040, 2371-2404
- Contains: `VFDR()`, `assem_vfdr()`, `DrugFinder()`, `getinfo()`
- Tools: abricate (CARD, VFDB), RGI, staramr, ResFinder

**Reporting Layer:**
- Purpose: Result compilation and visualization
- Location: Lines 3316-3478
- Contains: `combine_func()`, result aggregation functions
- Output: TSV files, HTML reports, Krona charts, CGView visualizations

## Data Flow

**Main Pipeline Flow:**

1. **Input Processing:**
   - Detect input type (FASTQ directory, FAST5, POD5, FASTA, barcode folder)
   - Optional basecalling (FAST5/POD5 → FASTQ)
   - Create sample working directories

2. **Quality Control:**
   - Adapter/barcode trimming (porechop)
   - Quality/length filtering (fastp, seqkit)
   - Host decontamination (hostile, kneaddata)
   - QC statistics generation (fastqc, seqkit stat)

3. **Taxonomic Classification:**
   - Kraken2/Bracken classification
   - Species abundance calculation
   - Contamination detection

4. **Assembly:**
   - Read subsampling (rasusa) for performance
   - De novo assembly (method-specific)
   - Assembly polishing (medaka)
   - Quality assessment (CheckM2)

5. **Annotation:**
   - Gene prediction (Prokka)
   - Functional annotation
   - Plasmid detection (PlasFlow, plasmidfinder)

6. **Specialized Analysis:**
   - MLST typing
   - Serotype prediction (species-specific)
   - AMR gene detection
   - Virulence factor identification
   - Mobile element prediction

7. **Result Compilation:**
   - File organization into structured directories
   - Summary table generation
   - Visualization creation (Krona, CGView, mummer2circos)
   - R-based report generation

**State Management:**
- File-based state tracking (e.g., `QC_ok`, `kk2_ok`, `SV_ok` files)
- Directory-based isolation of sample processing
- Global variables for configuration sharing

## Key Abstractions

**Sample Processing:**
- Purpose: Encapsulate all operations for a single sample
- Pattern: Each sample gets dedicated subdirectory under `fastq_analysis/`
- Entry Point: `main_process()` function

**Modular Flow Steps:**
- Purpose: Allow selective execution of pipeline stages
- Pattern: Flow list configuration with conditional execution
- Key Function: `run_fastq_flow()`, `run_fasta_flow()`

**Tool Wrappers:**
- Purpose: Standardize external tool execution
- Pattern: Subprocess calls with consistent logging and error handling
- Example: `subprocess.run()` with shell=True throughout

**Data Classes:**
- Purpose: Structured configuration for QC targets
- Implementation: `@dataclass(frozen=True)` for `QCTarget`
- Usage: Standardized QC processing for different read types

## Entry Points

**Command-Line Entry Point:**
- Location: Lines 3978-4249 (main execution block)
- Triggers: Direct script execution
- Responsibilities: Input validation, directory setup, sample iteration

**Primary Processing Function:**
- Function: `main_process()` (Lines 3845-3874)
- Parameters: 15+ arguments covering all pipeline options
- Branches: FastQ pipeline vs. FastA pipeline

**Flow Dispatch Functions:**
- FastQ: `run_fastq_pipeline()` (Lines 3728-3761)
- FastA: `run_fasta_pipeline()` (Lines 3809-3838)
- Fake: `run_fake_pipeline()` (Lines 3583-3609) - for testing

## Error Handling

**Strategy:** Try-except blocks with graceful degradation

**Patterns:**
- Try-except around main sample processing to prevent total pipeline failure
- File existence checks before processing steps
- Empty result handling with default values
- Subprocess error capture via log files

**Examples:**
- Lines 4100-4104: Sample-level exception handling
- Lines 3744-3748: Flow step error isolation
- Lines 2823-2834: Species inference fallback logic

## Cross-Cutting Concerns

**Logging:**
- Approach: Dedicated log files per processing step
- Pattern: `with open('step.log', 'w') as logf:`
- Tools: Standard Python logging for QC functions

**Validation:**
- Approach: File existence and size checks
- Pattern: `os.path.isfile()`, `os.path.getsize()` checks
- Examples: Lines 767-768, 1720-1722

**Resource Management:**
- Approach: CPU thread limiting, GPU memory checking
- Functions: `get_free_gpu_memory()`, automatic thread capping
- Pattern: Lines 72-75, 3950-3953

**Progress Tracking:**
- Approach: Standardized progress printing
- Function: `print_progress()` (Lines 3511-3516)
- Output: Step count, sample count, runtime, status messages

## Notable Architectural Decisions

1. **Monolithic Design:** Single 4,200-line file rather than module separation
2. **Global State:** Heavy use of global variables for configuration
3. **Subprocess-Heavy:** Minimal use of Python libraries; external tool orchestration
4. **Bilingual Comments:** Mixed Chinese and English documentation
5. **Hardcoded Paths:** Absolute paths to databases and tools throughout
6. **File-Based Communication:** State persistence via marker files

---

*Architecture analysis: 2026-03-16*
