# Testing Patterns

**Analysis Date:** 2026-03-16

## Test Framework Status

**Current State:** NO TESTS DETECTED

This codebase has no automated testing infrastructure. No test files, test directories, or testing configuration were found.

## Missing Testing Infrastructure

**Test Framework:** None
- No pytest configuration
- No unittest modules
- No test directory structure
- No CI/CD testing pipeline

**Test Commands:** Not applicable
```bash
# No test commands available
# pytest - Not installed
# python -m unittest - No test modules
```

## Code Structure Challenges for Testing

### Monolithic Architecture
- Single file: `Bac_assemble_260112_newformat.py` (4,367 lines)
- No module separation
- Global state dependencies
- Side effects at import time

### Tight Coupling Issues

**1. Subprocess Coupling**
- 415+ external tool dependencies
- Direct shell command execution
- No abstraction layer for external tools
```python
# Example of untestable pattern:
subprocess.run(f'seqkit seq {tref} > tmp_ref.fa', shell=True)
```

**2. File System Dependencies**
- Hardcoded absolute paths
- External database dependencies
- Temporary file manipulation
```python
Krdb='/home/dell/kraken2_custom_202101_24G'
speciesdb = pd.read_table('/data1/shanghai_pip/meta_genome/pathotable.tsv')
```

**3. Global State**
- Arguments parsed at module level
- Configuration loaded on import
- Logger instances created globally

### No Testable Units

Functions are not designed for unit testing:
- Large functions with multiple responsibilities
- Direct I/O operations
- No dependency injection
- Side effects (file system, subprocess, global state)

## Manual Testing Evidence

**Print Statements as Debugging:**
```python
print(command_line)  # 打印命令行参数
print('测试rna建库')  # Debug output
print(f'{Sam}文件夹内没有三代fastq格式文件')
```

**Fake Mode Implementation:**
```python
parser.add_argument('--fake_pip','-f',type=int,default=0,help='是否是假流程')
```
- Suggests manual testing via "fake" execution path
- No automated validation of fake mode

## Recommended Testing Strategy

### Phase 1: Infrastructure Setup

**1. Create Test Structure:**
```
tests/
├── __init__.py
├── conftest.py
├── unit/
│   ├── __init__.py
│   ├── test_qc.py
│   ├── test_assembly.py
│   └── test_utils.py
├── integration/
│   ├── __init__.py
│   └── test_pipeline.py
└── fixtures/
    ├── sample.fastq
    ├── sample.fasta
    └── sample.kraken.txt
```

**2. Add pytest Configuration:**
```ini
# pytest.ini
[pytest]
testpaths = tests
python_files = test_*.py
python_functions = test_*
addopts = -v --tb=short
```

**3. Install Testing Dependencies:**
```bash
pip install pytest pytest-cov pytest-mock freezegun
```

### Phase 2: Extract Testable Units

**Refactor for Testability:**

1. **Extract pure functions:**
```python
# Current (untestable)
def proc_kra(kraken, tax, lel):
    kradb = pd.read_table(kraken, header=None)
    # ... processing ...
    return result

# Refactored (testable)
def proc_kra_data(kradb: pd.DataFrame, tax: str, lel: str) -> pd.DataFrame:
    """Pure function operating on DataFrame."""
    # ... processing ...
    return result
```

2. **Create abstraction layers:**
```python
class SubprocessRunner:
    def run(self, cmd: str) -> subprocess.CompletedProcess:
        return subprocess.run(cmd, shell=True, capture_output=True)

# In tests:
class MockSubprocessRunner:
    def __init__(self, expected_outputs):
        self.expected_outputs = expected_outputs

    def run(self, cmd: str) -> MockCompletedProcess:
        return self.expected_outputs.get(cmd)
```

3. **Dependency injection:**
```python
def denovo_asb(inf, fq1, fq2, runner: SubprocessRunner = None):
    runner = runner or SubprocessRunner()
    # Use runner instead of subprocess.run directly
```

### Phase 3: Critical Test Coverage

**High Priority Tests:**

1. **QC Functions:**
   - `is_fasta()` - File format detection
   - `is_fastq()` - Extension validation
   - `safe_read_json()` - Error handling
   - `QCTarget` class - Data structure validation

2. **Data Processing:**
   - `proc_kra()` - Kraken report processing
   - `exreadsID()` - Read ID extraction logic
   - `normalize_summary_txt()` - File normalization

3. **Utility Functions:**
   - `safe_float()`, `safe_int()` - Type conversion
   - `pick()` - Dictionary utilities
   - `get_logger()` - Logging setup

### Phase 4: Integration Testing

**End-to-End Pipeline Tests:**

Use small test datasets to verify pipeline stages:
```python
@pytest.fixture
def sample_fastq(tmp_path):
    """Create minimal test FASTQ file."""
    fastq = tmp_path / "test.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")
    return fastq

def test_qc_stage(sample_fastq, tmp_path):
    """Test QC stage with real file."""
    result = run_qc(sample_fastq, output_dir=tmp_path)
    assert (tmp_path / "QC.summary.tsv").exists()
```

## Testing Tools Recommendation

**Core Framework:**
- `pytest` - Test runner
- `pytest-cov` - Coverage reporting
- `pytest-mock` - Mocking utilities
- `freezegun` - Time freezing

**Bioinformatics Testing:**
- `pysam` - Manipulate SAM/BAM/CRAM for testing
- `biopython` - Sequence manipulation in tests
- `pandas.testing` - DataFrame comparison

**Coverage Targets:**
- Initial target: 30% coverage of utility functions
- Medium target: 60% coverage including QC logic
- Long-term: 80% coverage with integration tests

## Risk Assessment

**Current Risks (No Testing):**
- Silent failures in data processing
- Undetected regressions in pipeline stages
- No validation of scientific correctness
- Deployment risks without verification

**Critical Paths Needing Tests:**
1. QC calculation logic (q20/q30 rates)
2. Species identification from Kraken reports
3. Assembly success/failure detection
4. File format validation
5. Database path resolution

## Testing Checklist

**Immediate Actions:**
- [ ] Add pytest to project dependencies
- [ ] Create tests/ directory structure
- [ ] Extract first pure function for testing
- [ ] Add CI workflow for automated testing

**Short-term Goals:**
- [ ] Unit tests for all utility functions
- [ ] Mock-based tests for subprocess calls
- [ ] Fixture-based tests for file I/O
- [ ] Coverage reporting in CI

**Long-term Goals:**
- [ ] Integration tests with small datasets
- [ ] Regression test suite for bug fixes
- [ ] Performance benchmarking tests
- [ ] Property-based testing for data processing

---

*Testing analysis: 2026-03-16*

**Summary:** This codebase requires significant refactoring to support automated testing. The monolithic architecture and tight coupling to external tools make unit testing impossible without abstraction layers. Immediate focus should be on extracting pure functions and creating integration tests with small datasets.
