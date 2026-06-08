# Repository Boundaries

This repository should contain source code, lightweight fixtures, scripts,
configuration templates, and human-written documentation.

It should not contain machine-local state, generated build outputs, virtual
environments, large downloaded databases, packaged desktop applications, or
task/output directories produced by running analyses.

## Keep In Git

- Python source under `bac_analysis_portal/`, `metagenomic_refactor/`,
  `pathosource_refactor/`, and `scripts/`
- Flask templates and static assets required by the source application
- Requirements files, build scripts, and `.spec` templates only when they are
  edited as source
- Small test fixtures and curated demo inputs that are required by tests or
  README walkthroughs
- Documentation and architecture notes

## Keep Out Of Git

- `dist/`, `build/`, `target/`, `Windows_app/`, and generated installers
- PyInstaller caches and bundled runtime libraries
- `.venv*`, `envs/`, `__pycache__/`, `.pytest_cache/`, and `.DS_Store`
- SQLite portal/task state such as `bac_analysis_portal.sqlite3` and
  `bac_analysis_portal/*.db`
- Large reference databases under `database/`, `host_database/`, `genome_db/`,
  `pathogen_database/`, `snpEff/`, and `vadr/`
- Downloaded browser/viewer/tool bundles under `public/` and `soft/`
- Analysis outputs, scratch runs, generated batch inputs, and delivery archives

## Asset Boundary

Large databases and third-party tool bundles are deployment assets, not source
files. They should be provisioned by documented download/build scripts or copied
from an approved deployment location. If a lightweight fixture is needed for a
test, place only the minimal fixture in `tests/` or a documented demo directory.

## Review Rule

Before adding a binary, database, archive, generated report, or packaged app to
Git, ask whether a clean clone needs that exact file to understand or verify the
source. If the answer is "no", keep it outside the repository and document how
to recreate or fetch it.
