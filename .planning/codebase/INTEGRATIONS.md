# External Integrations

**Analysis Date:** 2026-03-16

## Overview

This is a bacterial genome assembly and analysis pipeline that orchestrates 30+ external bioinformatics tools through subprocess calls. It integrates with specialized genomic databases and uses conda environments for tool management.

## External Tools & APIs

**Genome Assembly (Long-read):**
- **Flye** (conda: base) - De novo assembly for ONT/PacBio
  - Usage: `flye --nano-hq/--nano-raw/--pacbio-hifi/--pacbio-raw`
  - Lines: 974-984
- **Canu** (conda: base) - Long-read assembly
  - Usage: `canu -nanopore-raw/-pacbio-raw`
  - Lines: 1001-1004
- **WTDBG2** (conda: base) - De novo assembler
  - Usage: `wtdbg2.pl` (Perl script)
  - Line: 997
- **Miniasm** (conda: base) - Ultrafast assembly
  - Usage: `minimap2 | miniasm | gfatools`
  - Lines: 988-995
- **Unicycler** (conda: base) - Hybrid assembly
  - Usage: `unicycler -l/-1 -2`
  - Lines: 1006-1007, 1108+
- **Raven** (conda: base) - ONT assembly
  - Usage: `raven`
  - Line: 1009

**Genome Assembly (Short-read):**
- **SPAdes** (conda: base) - Short-read assembly
  - Usage: `spades.py --pe1-1 --pe1-2`
  - Lines: 1020-1024
- **MEGAHIT** (conda: BASALT) - Metagenomic assembly
  - Usage: `megahit -1 -2`
  - Lines: 1043-1045

**Metagenomic Binning:**
- **BASALT** (conda: BASALT) - Automated binning
  - Usage: `BASALT --module autobinning/refinement`
  - Lines: 1061-1072
- **dRep** (conda: BASALT) - Dereplication
  - Usage: `dRep dereplicate`
  - Line: 1080

**Polishing:**
- **Medaka** (conda: medaka) - ONT consensus
  - Usage: `medaka_consensus`
  - Lines: 755-767
- **FreeBayes** (conda: base) - Variant calling (polishing)
  - Usage: `freebayes-parallel`
  - Lines: 1232, 1265, 2797
- **Clair3** (conda: clair3) - ONT variant calling
  - Usage: `run_clair3.sh --platform=ont`
  - Lines: 1265

**Quality Control:**
- **Fastp** (conda: base) - Read preprocessing
  - Usage: `fastp --in1 --in2 --out1 --out2`
  - Lines: 549+, 561+, 570+, 652+, 717+
- **FastQC** (conda: base) - Quality reports
  - Usage: `fastqc`
  - Lines: 620, 687-688, 742
- **SeqKit** (conda: base) - Sequence manipulation
  - Usage: `seqkit seq/grep/stat/fx2tab`
  - Lines: 106, 523-534, 561, 596, 638, 650, 715, 755, 765-767
- **Porechop** (conda: base) - Adapter trimming (ONT)
  - Usage: `porechop -i -o`
  - Line: 561
- **Rasusa** (conda: base) - Read subsampling
  - Usage: `rasusa reads --bases`
  - Lines: 546-548, 647-649, 712-714

**Taxonomic Classification:**
- **Kraken2** (conda: base) - Taxonomic assignment
  - Usage: `kraken2 --db`
  - Lines: 626, 694
  - Database: `/data1/shanghai_pip/meta_genome/database/kraken2_custom_202101_24G`
- **Bracken** (conda: base) - Abundance estimation
  - Usage: `bracken -d -i -o`
  - Lines: 627, 695

**Host Decontamination:**
- **Hostile** (conda: hostile) - Host read removal
  - Usage: `hostile clean --aligner minimap2/bowtie2`
  - Lines: 567-568, 677-679, 735-736
- **KneadData** (conda: kneaddata) - Quality control
  - Usage: `kneaddata --bypass-trf --bypass-trim`
  - Lines: 683-684, 740-741

**Genome Quality Assessment:**
- **CheckM2** (conda: cm2/cm210) - Completeness/contamination
  - Usage: `checkm2 predict`
  - Lines: 842, 1701-1705, 1869
- **GTDB-Tk** (conda: gtdbtk) - Taxonomic classification
  - Usage: `gtdbtk classify_wf`
  - Lines: 836, 893, 953

**Genome Annotation:**
- **Prokka** (conda: base) - Bacterial annotation
  - Usage: `prokka --force --addgenes`
  - Lines: 1641-1662, 1731, 1899

**Antimicrobial Resistance:**
- **Abricate** (conda: base) - Resistance gene screening
  - Usage: `abricate --db vfdb/card`
  - Lines: 851, 853
- **RGI** (conda: RGI_new) - CARD resistance genes
  - Usage: `rgi main --include_loose`
  - Lines: 855, 911-914
- **StarAMR** (conda: base) - ResFinder/PlasmidFinder
  - Usage: `staramr search`
  - Lines: 858, 870-871, 915, 923

**Plasmid Detection:**
- **PlasFlow** (conda: plasflow) - Plasmid prediction
  - Usage: `PlasFlow.py`
  - Lines: 869, 1599, 1812, 1869

**MLST Typing:**
- **MLST** (conda: base) - Sequence typing
  - Usage: `mlst --quiet --csv`
  - Lines: 3137-3143, 3159

**Visualization:**
- **GenoVi** (conda: genovi) - Circular genome maps
  - Usage: `genovi -i -o`
  - Lines: 1731, 1899

**Read Mapping:**
- **Minimap2** (conda: base) - Alignment
  - Usage: `minimap2 -x ava-ont/map-pb/map-hifi`
  - Lines: 988-992
- **CoverM** (conda: coverm) - Coverage calculation
  - Usage: `coverm genome`
  - Line: 887

## Data Storage

**Databases (External):**
- **Kraken2 Database**: `/data1/shanghai_pip/meta_genome/database/kraken2_custom_202101_24G` (24GB)
- **GTDB-Tk Database**: Referenced but path not shown
- **CheckM2 Database**: `/data1/shanghai_pip/meta_genome/uniref100.KO.1.dmnd`
- **Hostile Index**: `/data/deploy/meta_new/Database/Host_Ref/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401.mmi`
- **KneadData Databases**:
  - `/data/Ref/human_hg38_refMrna`
  - `/data/Ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA`
- **MLST Database**: `/home/dell/miniconda3/envs/TB_ONT/db/pubmlst/`
- **Species Reference**: `/data1/shanghai_pip/meta_genome/pathotable.tsv`

**File Storage:**
- Local filesystem only (no cloud storage)
- Hardcoded paths: `/data/`, `/data1/`, `/home/dell/`
- Temporary directories: `tmpdir1` (default)

**No Traditional Databases:**
- No SQL (SQLite, PostgreSQL, MySQL)
- No NoSQL (MongoDB, Redis)
- File-based data exchange (TSV, CSV, FASTA, FASTQ, JSON)

## Authentication & Identity

**Not Applicable:**
- No user authentication system
- No API keys or tokens
- No OAuth integration
- Single-user command-line tool

## Monitoring & Observability

**Logging:**
- Python `logging` module with custom formatter
- Pattern: `[%(asctime)s] %(levelname)s %(name)s: %(message)s`
- Functions: `get_logger()` defined at lines 178, 201
- Log files created per tool (e.g., `gtdbtk.log`, `checkm2.log`)

**Error Tracking:**
- File existence checks before operations
- Try-catch blocks for JSON parsing (lines 190-197, 213-220)
- No external error tracking (Sentry, etc.)

**Progress Tracking:**
- `sys.stdout.flush()` after command echo
- Print statements for command execution

## CI/CD & Deployment

**Not Detected:**
- No CI configuration (GitHub Actions, GitLab CI, Jenkins)
- No containerization (Dockerfile, Singularity)
- No package management (pip requirements, environment.yml)
- Manual deployment to server

**Deployment Pattern:**
- Direct server deployment at `/data/deploy/`
- Conda environment-based tool management
- Hardcoded absolute paths throughout

## Webhooks & Callbacks

**None:**
- No incoming webhooks
- No outgoing callbacks
- No REST API integration
- No HTTP/FTP network calls

## Environment Configuration

**Required Environment Variables:**
- None explicitly required
- All configuration hardcoded

**Conda Environments Used:**
- `TB_ONT` - Main Python environment
- `medaka` - Consensus polishing
- `clair3` - Variant calling
- `hostile` - Host decontamination
- `kneaddata` - Quality control
- `BASALT` - Metagenomic binning
- `gtdbtk` - Taxonomic classification
- `cm2` / `cm210` - Genome quality
- `RGI_new` - Resistance genes
- `plasflow` - Plasmid detection
- `genovi` - Visualization
- `coverm` - Read coverage

**Hardcoded Paths:**
- `/home/dell/miniconda3/` - Conda installation
- `/data1/shanghai_pip/meta_genome/` - Pipeline data
- `/data/deploy/` - Deployment directory
- `/data/Ref/` - Reference databases

---

*Integration audit: 2026-03-16*
