# Technology Stack

**Analysis Date:** 2026-03-16

## Languages

**Primary:**
- Python 3 (TB_ONT conda environment) - Complete codebase

**Secondary:**
- Shell/Bash - External tool orchestration via subprocess
- Perl - Used in external tools (wtdbg2)

## Runtime

**Environment:**
- Conda (Miniconda3) with Python 3.x
- Shebang: `#!/home/dell/miniconda3/envs/TB_ONT/bin/python`
- Platform: Linux (Dell server deployment)

**Package Manager:**
- Conda (primary) - Environment and package management
- pip (implied) - Python package installation

## Frameworks

**Core Scientific Computing:**
- **Biopython** (Bio) - Sequence I/O and phylogenetic analysis
- **pandas** - Data manipulation and analysis
- **numpy** - Numerical computing

**Bioinformatics Specific:**
- **pathogenprofiler** - Pathogen genomic profiling
- **pytaxonkit** - Taxonomic classification utilities

**Standard Library (Heavy Usage):**
- `subprocess` - External tool execution (primary orchestration method)
- `argparse` - Command-line interface
- `multiprocessing` - Parallel processing
- `logging` - Structured logging
- `pathlib`/`os.path` - Path manipulation
- `json`/`csv` - Data serialization
- `dataclasses` - Data structures
- `typing` - Type hints (Optional, Dict, List, Tuple, Set, Any)

## Key Dependencies

**Critical:**
- `pandas` - DataFrame operations throughout pipeline
- `numpy` - Numerical operations
- `biopython` - FASTA/FASTQ parsing, phylogenetics
- `pathogenprofiler` - Custom pathogen analysis
- `pytaxonkit` - Taxonomic operations

**Infrastructure:**
- `subprocess` - Shell command execution (100+ calls)
- `glob` - File pattern matching
- `shutil` - File operations

**Data Processing:**
- `re` - Regular expressions
- `itertools` - Iterator operations
- `math` - Mathematical operations
- `ast` - Abstract syntax trees

## Configuration

**Environment:**
- Hardcoded paths to databases and tools
- Conda environment activation via `/home/dell/miniconda3/bin/conda run -n <env>`
- Database paths in `/data/` and `/data1/shanghai_pip/` directories

**Key Configuration Points:**
- Line 77: `Krdb='/home/dell/kraken2_custom_202101_24G'`
- Line 93: `speciesdb = pd.read_table('/data1/shanghai_pip/meta_genome/pathotable.tsv')`
- Line 755: `medaka_cmd = "/home/dell/miniconda3/bin/conda run -n medaka medaka_consensus"`

**No External Config Files:**
- No YAML/JSON/TOML configuration
- All paths hardcoded in source
- Database locations embedded throughout

## Platform Requirements

**Development:**
- Linux environment (CentOS/Ubuntu)
- Miniconda3 installed at `/home/dell/miniconda3`
- 10+ CPU cores recommended (default thread count)
- Large storage for genomic databases (>100GB)

**Production:**
- High-performance computing cluster
- 24GB+ RAM for genome assembly
- Conda environments for tool isolation
- Access to `/data/` and `/data1/` mount points

**Compute Resources:**
- Default 10 threads (`-t` parameter)
- Auto-detects CPU count and caps accordingly
- Memory-intensive operations (Flye, Canu, SPAdes)

---

*Stack analysis: 2026-03-16*
