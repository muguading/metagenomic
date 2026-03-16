# Codebase Concerns

**Analysis Date:** 2026-03-16

## Tech Debt

**Monolithic Architecture:**
- Issue: Single 4,367-line file (`Bac_assemble_260112_newformat.py`) containing 115 functions with no module separation
- Files: `Bac_assemble_260112_newformat.py`
- Impact: Impossible to test individual components, difficult navigation, high cognitive load
- Fix approach: Extract into modules by concern (QC, assembly, annotation, reporting)

**Global State Pollution:**
- Issue: 30+ global variables defined at module level (lines 72-115), including `argv`, `nt`, `Krdb`, `speciesdb`
- Files: `Bac_assemble_260112_newformat.py` (lines 72-115)
- Impact: Functions have hidden dependencies, impossible to reason about execution context, race conditions in parallel execution
- Fix approach: Pass configuration explicitly via dataclasses or dependency injection

**Hardcoded Paths:**
- Issue: 15+ absolute paths to specific servers and databases
- Files: `Bac_assemble_260112_newformat.py`
- Examples:
  - Line 77: `Krdb='/home/dell/kraken2_custom_202101_24G'`
  - Line 93: `speciesdb = pd.read_table('/data1/shanghai_pip/meta_genome/pathotable.tsv')`
  - Line 95: `sc1='/data/deploy/TB_soft/other_soft/3_kreport2krona.py'`
  - Line 567: `'/home/dell/miniconda3/bin/conda run -n hostile hostile...'`
  - Line 1701: `'/data1/shanghai_pip/meta_genome/uniref100.KO.1.dmnd'`
- Impact: Code only runs on specific server, breaks in any other environment
- Fix approach: Externalize to configuration file or environment variables

**Shell Injection Vulnerabilities:**
- Issue: 410+ subprocess calls with `shell=True` using string interpolation with user input
- Files: `Bac_assemble_260112_newformat.py` (throughout)
- Examples:
  - Line 106: `subprocess.run(f'seqkit seq {tref} > tmp_ref.fa',shell=True)`
  - Line 526: `subprocess.run(f'seqkit grep -n -f {Maintax}_fq1ID.txt {fq1} > {Pre}.R1.fastq',shell=True)`
- Impact: Command injection if any input variable contains shell metacharacters
- Fix approach: Use list arguments instead of shell strings, validate inputs

**Bare Except Clauses:**
- Issue: 19 bare `except:` clauses catching all exceptions silently
- Files: `Bac_assemble_260112_newformat.py`
- Lines: 174, 1734, 1902, 2053, 2088, 3160, 3298, 4001, 4030, 4055, 4072, 4100, 4131, 4159, 4188, 4231, 4249, 4310, 4339
- Impact: Swallows KeyboardInterrupt, SystemExit, and programming errors; makes debugging impossible
- Fix approach: Catch specific exceptions, add logging before re-raising

**Duplicate Function Definitions:**
- Issue: Two identical `get_logger` functions (lines 178-186 and 201-209) and two `safe_read_json` functions (lines 190-198 and 213-220)
- Files: `Bac_assemble_260112_newformat.py`
- Impact: Second definition overwrites first, confusing behavior if code assumes different implementations
- Fix approach: Remove duplicates, consolidate utility functions

**Magic Numbers and Strings:**
- Issue: Hardcoded thresholds, file extensions, and database IDs throughout
- Examples:
  - Line 634: `if float(tmpfile.loc[tmpfile['name']==ONTSpe,'fraction_total_reads'].tolist()[0]) < float(tspeabun)`
  - Line 546: `subprocess.run(f'rasusa reads --bases 2gb -o {Pre}.sub.fastq {inf}',shell=True)`
  - Line 548: `subprocess.run(f'rasusa reads --bases 20gb -o {Pre}.sub.fastq {inf}',shell=True)`
- Impact: Business logic scattered in code, difficult to adjust parameters
- Fix approach: Centralize configuration in constants or config objects

## Known Issues

**Race Condition in File Polling:**
- Issue: `wait_for_file()` function (lines 773-783) has incorrect size comparison logic
- Files: `Bac_assemble_260112_newformat.py` (lines 773-783)
- Bug: Line 779 compares `filepath == lastsize` (string vs int) instead of `curruent_size == lastsize`
- Impact: Function returns immediately instead of waiting for file to stabilize
- Fix approach: Correct variable name typo (`curruent_size` -> `current_size`)

**Silent Failures in QC:**
- Issue: QC functions check for file existence but not process success
- Files: `Bac_assemble_260112_newformat.py` (lines 533-750)
- Impact: Failed subprocess calls continue as if successful, producing invalid downstream results
- Fix approach: Check subprocess return codes, raise on failure

**Incorrect Sleep Call:**
- Issue: Line 775 calls `times.sleep(1)` instead of `time.sleep(1)`
- Files: `Bac_assemble_260112_newformat.py` (line 775)
- Impact: `NameError` exception in wait loop (currently swallowed by bare except)
- Fix approach: Fix typo, remove bare except to expose errors

**Resource Leaks:**
- Issue: File handles opened without context managers in some paths
- Files: `Bac_assemble_260112_newformat.py` (lines 748, 795, 802, 803, 817, 826-830)
- Examples: `open('QC_ok','w').write('已跑过')`, `open(f'binning_name.tsv','a').write(...)`
- Impact: File descriptors exhausted in long-running processes
- Fix approach: Use `with` statement for all file operations

**Unreliable Process Detection:**
- Issue: Uses file existence to determine if steps completed, not success
- Files: `Bac_assemble_260112_newformat.py` (throughout)
- Examples: Lines 1390, 1455, 1648, 1746 checking `os.path.isfile()`
- Impact: Partial/corrupted files treated as successful completion
- Fix approach: Use success markers, checksums, or database state tracking

## Security Considerations

**Command Injection:**
- Risk: Critical - user-provided filenames and paths directly interpolated into shell commands
- Files: Throughout `Bac_assemble_260112_newformat.py`
- Example: Line 106 `subprocess.run(f'seqkit seq {tref} > tmp_ref.fa',shell=True)` - if `tref` contains `; rm -rf /`, code executes it
- Current mitigation: None
- Recommendations:
  - Use `subprocess.run(['seqkit', 'seq', tref], ...)` without shell=True
  - Validate all inputs against whitelist patterns
  - Escape shell arguments if shell=True absolutely required

**World-Readable Temporary Files:**
- Risk: Medium - temporary files created in shared directories with default permissions
- Files: Lines 106, 107 creating `tmp_ref.fa` in current working directory
- Impact: Sensitive genomic data potentially exposed to other users
- Current mitigation: None
- Recommendations: Use `tempfile` module with appropriate permissions

**Hardcoded Credentials in Paths:**
- Risk: Low - username in paths reveals system information
- Files: Line 77, 567, 735, 1701 referencing `/home/dell/`
- Impact: Information disclosure, social engineering vector
- Current mitigation: None
- Recommendations: Externalize all paths to config

## Performance Bottlenecks

**Busy-Waiting Pattern:**
- Problem: `wait_for_file()` polls filesystem every 2 seconds indefinitely
- Files: `Bac_assemble_260112_newformat.py` (lines 773-783)
- Impact: Wasted CPU cycles, potential infinite loops
- Improvement path: Use inotify or file system events instead of polling

**Inefficient Data Structures:**
- Problem: Repeated DataFrame filtering and list conversions in loops
- Files: `Bac_assemble_260112_newformat.py` (lines 628-634, 1709-1726)
- Impact: O(n²) operations where O(n) would suffice
- Improvement path: Use vectorized pandas operations, cache lookups

**Repeated File System Operations:**
- Problem: Multiple `os.path.isfile()` checks for same files in tight loops
- Files: `Bac_assemble_260112_newformat.py` (assembly functions)
- Impact: Unnecessary I/O overhead
- Improvement path: Cache existence checks, use pathlib for efficiency

**No Parallelization:**
- Problem: Sequential processing of independent samples
- Files: Main execution block (lines 4260-4367)
- Impact: Underutilizes multi-core systems
- Improvement path: Use multiprocessing.Pool for sample-level parallelism

## Fragile Areas

**Assembly Function Complexity:**
- Files: `asb_func()` (lines 1385-1740) - 355 lines
- Why fragile: 15 nested conditionals, multiple assembly methods, mixed concerns
- Safe modification: Add unit tests before any changes, extract method-specific logic
- Test coverage: None detected

**Serotyping Logic:**
- Files: Lines 2842-2971 with complex primer matching rules
- Why fragile: 47 hardcoded serotype rules, manual primer name normalization
- Safe modification: Add validation tests with known positive/negative samples
- Test coverage: None detected

**QC Summary Statistics:**
- Files: `summaryfastqc_prod()` and helpers (lines 400-480)
- Why fragile: Fragile parsing of FASTQC output formats, assumes column positions
- Safe modification: Use structured formats (JSON) instead of text parsing
- Test coverage: None detected

**Database Schema Dependencies:**
- Files: Throughout (speciesdb, vfmeta, speciesrefdb loading)
- Why fragile: Assumes specific column names and positions in TSV files
- Safe modification: Add schema validation on load, version database files
- Test coverage: None detected

**Input Type Detection:**
- Files: `check_input()` and main block (lines 4200-4367)
- Why fragile: String-based file type detection, relies on extensions
- Safe modification: Use file magic numbers, validate format with parsers
- Test coverage: None detected

## Scaling Limits

**Memory Usage:**
- Current capacity: Loads entire DataFrames into memory (lines 93, 114, 115, 597, 628)
- Limit: Pandas memory overhead for large genomic datasets
- Scaling path: Use chunked reading, Dask for out-of-core processing

**Concurrent Sample Processing:**
- Current capacity: Sequential processing only
- Limit: Single sample processed at a time regardless of available cores
- Scaling path: Implement job queue (Celery, Snakemake, or Nextflow)

**Database Size:**
- Current capacity: Fixed Kraken2/Bracken database paths
- Limit: Cannot switch databases without code modification
- Scaling path: Make databases configurable per-project

## Dependencies at Risk

**Conda Environment Coupling:**
- Risk: Hardcoded conda environment names (`medaka`, `hostile`, `cm2`, `genovi`, etc.)
- Impact: 15+ external tools with specific environment requirements
- Migration plan: Containerization (Docker/Singularity) with environment.yml

**Unpinned Dependencies:**
- Risk: No requirements.txt, setup.py, or environment.yml detected
- Impact: Different versions of pandas/biopython may break API compatibility
- Migration plan: Add requirements.txt with pinned versions

**Deprecated Bioinformatics Tools:**
- Risk: Some tools (prokka) in maintenance mode
- Impact: No bug fixes or updates for new sequence types
- Migration plan: Evaluate bakta as prokka replacement

## Missing Critical Features

**No Test Suite:**
- Problem: Zero unit tests, integration tests, or validation suite
- Blocks: Safe refactoring, CI/CD implementation
- Priority: Critical

**No Logging Configuration:**
- Problem: Mix of print statements, basic logging, and log files
- Blocks: Production monitoring, debugging distributed runs
- Priority: High

**No Input Validation:**
- Problem: No schema validation for input files or parameters
- Blocks: Early error detection, user-friendly error messages
- Priority: High

**No Resume Capability:**
- Problem: Partial runs leave inconsistent state
- Blocks: Long-running pipeline reliability
- Priority: Medium

**No Resource Monitoring:**
- Problem: No tracking of CPU, memory, or disk usage
- Blocks: Optimization, cost accounting, failure diagnosis
- Priority: Medium

## Test Coverage Gaps

**Assembly Algorithms:**
- What's not tested: All 9 assembly methods (flye, canu, unicycler, etc.)
- Files: `denovo_asb()`, `reassm_fun()` functions
- Risk: Silent failures produce incorrect assemblies
- Priority: Critical

**Database Queries:**
- What's not tested: Kraken2/Bracken parsing logic
- Files: `proc_kra()`, `proc_kra1()` functions
- Risk: Incorrect taxonomic classification
- Priority: High

**File Format Parsers:**
- What's not tested: FASTQ, FASTA, GFF, GBK parsing edge cases
- Files: Throughout
- Risk: Format variations cause crashes
- Priority: High

**Error Recovery:**
- What's not tested: Exception handling paths
- Files: All try/except blocks
- Risk: Except clauses never tested, may hide errors
- Priority: High

---

*Concerns audit: 2026-03-16*
