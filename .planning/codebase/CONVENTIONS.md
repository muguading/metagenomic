# Coding Conventions

**Analysis Date:** 2026-03-16

## Overview

This codebase is a Python-based bioinformatics pipeline for bacterial genome assembly and metagenomic analysis. The script `Bac_assemble_260112_newformat.py` is a monolithic workflow script (~4,367 lines) that orchestrates various bioinformatics tools through subprocess calls.

## Naming Patterns

**Files:**
- Main script: `Bac_assemble_260112_newformat.py` (descriptive with date suffix)
- No consistent module separation - single large file architecture

**Functions:**
- Mixed naming conventions:
  - `snake_case` for most functions: `is_fasta()`, `run_cmd()`, `get_logger()`
  - Abbreviated names: `proc_kra()`, `exreadsID()`, `ngs_qc()`
  - Descriptive names: `normalize_fastqc_images()`, `safe_read_json()`
  - Inconsistent abbreviation: `denovo_asb()` vs `reassm_fun()`

**Variables:**
- Short abbreviations common: `inf`, `ofn`, `nt`, `minl`, `minQ`
- Mixed languages (Chinese comments, English variables)
- Global configuration variables at module level

**Classes:**
- Single class detected: `QCTarget` (dataclass-style)
- Uses `@dataclass` decorator from `dataclasses` module

## Code Style

**Formatting:**
- No automated formatter detected (no Black, autopep8, or Ruff configuration)
- Inconsistent indentation spacing
- Line length varies significantly (some very long lines >150 characters)
- Mixed quote usage: single quotes predominantly, some double quotes

**Import Organization:**
```python
# Standard library imports (grouped at top)
import pandas as pd
import os
import subprocess
import argparse
import sys
import numpy as np
import time
import re

# Bioinformatics libraries
from Bio import SeqIO
from Bio import Phylo

# Utility imports
import itertools
import glob
import math
import json
import csv
import ast
from typing import Iterable, Set, Tuple, Dict, List
import shutil
import logging
from dataclasses import dataclass
from typing import Dict, Any, Optional
from pathlib import Path

# External bioinformatics tools
import pathogenprofiler as pp
import pytaxonkit
```

**Import Issues:**
- `glob` imported twice (lines 15 and 27)
- `typing` imports split (lines 21 and 26)
- No `__future__` annotations

## Type Hints

**Usage Pattern:**
- Partial adoption of type hints in newer code sections
- Example: `def safe_read_json(path: str, logger: logging.Logger) -> Optional[Dict[str, Any]]:`
- Mix of typed and untyped functions
- Uses `typing` module imports: `Iterable`, `Set`, `Tuple`, `Dict`, `List`, `Optional`, `Any`

**Type Hint Locations:**
- QC module functions (lines 190-440) have type hints
- Core pipeline functions lack type hints
- Inconsistent application even within same functional area

## Error Handling

**Patterns:**

1. **Bare except clauses (AVOID):**
```python
try:
    shutil.copy(f, dest)
except:
    pass
```

2. **Specific exception handling (PREFERRED):**
```python
try:
    with open(path) as f:
        return json.load(f)
except FileNotFoundError:
    logger.warning(f"JSON not found: {path}")
except json.JSONDecodeError as e:
    logger.error(f"JSON parse error: {path} ({e}")
return None
```

3. **Exception with logging:**
```python
try:
    # operation
except Exception as e:
    logger.error(f"Operation failed: {e}")
```

**Issues:**
- Many bare `except:` clauses throughout
- Silent failure patterns (pass in except blocks)
- Inconsistent error propagation

## Logging

**Framework:** Python standard `logging` module

**Pattern:**
```python
def get_logger(name="qc", level=logging.INFO):
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(level)
        h = logging.StreamHandler()
        fmt = logging.Formatter("[%(asctime)s] %(levelname)s %(name)s: %(message)s")
        h.setFormatter(fmt)
        logger.addHandler(h)
    return logger
```

**Usage:**
- Logger passed between functions explicitly
- Mix of `logger.info()`, `logger.warning()`, `logger.error()`
- Also uses `print()` statements for user-facing output
- Format: `[2024-03-16 10:30:00] INFO qc: message`

## Comments

**Style:**
- Predominantly Chinese comments
- Inline comments explain "what" not "why"
- Section dividers: `# ---------- logging ----------`
- Module docstrings absent

**Examples:**
```python
# 规定了应使用的Python解释器的路径
#!/home/dell/miniconda3/envs/TB_ONT/bin/python

# 创建一个名为'纯菌+宏测组装分析脚本'的解析器（对象）
parser = argparse.ArgumentParser(description='纯菌+宏测组装分析脚本')

# -i,输入文件路径。type=str，输入的值会被解析为'字符串'
parser.add_argument('--input','-i',type=str,default=False,help='输入文件路径')
```

## Function Design

**Size:**
- Very large functions common (500+ lines)
- `QC_func()`: ~200 lines
- `denovo_asb()`: ~170 lines
- `main_process()`: ~600+ lines

**Parameters:**
- Many functions take 10+ parameters
- Example: `main_process(fq1, fq2, threads, Pre, pts, pst, method, asmt, f, outputfa, ...)`
- Uses `**kwargs` pattern implicitly through `argv`

**Return Values:**
- Inconsistent return patterns
- Some functions return `bool`, others return `None`
- Silent failure returns `None` or empty structures

## Module Design

**Structure:**
- Single-file architecture (monolithic)
- No separation of concerns (I/O, logic, orchestration mixed)
- Global state through module-level variables

**Execution Flow:**
1. Global imports
2. Global constants and paths
3. Argument parsing (at module level)
4. Function definitions
5. Global configuration loading
6. Main execution block (at end of file)

**Entry Point:**
- No `if __name__ == "__main__":` guard
- Script executes immediately on import
- Relies on argument parsing side effects

## Subprocess Usage

**Pattern:**
- Heavy reliance on external bioinformatics tools
- 415+ subprocess calls throughout
- Pattern: `subprocess.run(f'command {variable}', shell=True)`

**Security Issues:**
- Uses `shell=True` universally
- No input sanitization detected
- F-string interpolation directly into shell commands

**Example:**
```python
subprocess.run(f'seqkit seq {tref} > tmp_ref.fa', shell=True)
subprocess.run(f'cat {tinf}/*.f*q* |seqkit rmdup -i |seqkit seq > {Sam}.raw.fastq', shell=True)
```

## Data Handling

**Pandas Usage:**
- Heavy use of pandas for TSV/CSV manipulation
- Pattern: `pd.read_table()` for input, `.to_csv()` for output
- DataFrames passed between functions
- Column access by integer position (fragile)

**File I/O:**
- Mix of `pathlib.Path` and `os.path` operations
- Hardcoded paths to external databases
- Temporary file creation without cleanup

## Configuration

**Approach:**
- Command-line arguments via `argparse`
- 25+ command-line parameters
- Hardcoded database paths in global scope
- Environment-dependent paths (e.g., `/data1/shanghai_pip/`)

**Configuration Variables:**
```python
Krdb='/home/dell/kraken2_custom_202101_24G'
speciesdb = pd.read_table('/data1/shanghai_pip/meta_genome/pathotable.tsv')
```

## Version Management

**Pattern:**
```python
__author__='wsh'
__version__='1.0.0'
__date__='20260306'
```

## Recommended Improvements

1. **Add pre-commit hooks** for formatting (Black, Ruff)
2. **Replace bare except clauses** with specific exception handling
3. **Add type hints** consistently across all functions
4. **Break into modules** by functional area (qc, assembly, annotation)
5. **Add `__main__` guard** for script execution
6. **Sanitize subprocess inputs** to prevent injection
7. **Use pathlib exclusively** instead of os.path
8. **Add unit tests** (currently none detected)

---

*Convention analysis: 2026-03-16*
