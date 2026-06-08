# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a pathogen microbial analysis platform (病原微生物分析工作台) consisting of a Flask web portal, refactored Python analysis pipelines, and standalone delivery packages. The target users are laboratory analysts who submit samples, monitor batch tasks, and review reports.

## Repository Structure

- `bac_analysis_portal/` — Flask web application for sample submission, task management, and result viewing. Uses SQLite for persistence.
- `metagenomic_refactor/` — Refactored metagenomic analysis pipeline package (QC, assembly, taxonomy, strain typing, virus analysis, MAG binning, reporting).
- `pathosource_refactor/` — Pathogen source tracking pipeline (phylogeny, MLST, cgMLST typing).
- `genome_db/` — Genome database management mini-app (Flask + SQLite).
- `delivery/neisseria_meningitidis_pipeline/` — Self-contained delivery pipeline for *Neisseria meningitidis* with its own web UI.
- `scripts/` — Utility scripts, build scripts, and data backfill scripts.
- `analysis_tasks/` — Runtime output directory for submitted analysis tasks.
- `database/` — Reference databases and knowledge base assets.

## Common Development Commands

### Run the web portal
```bash
bash run_analysis_portal.sh
# or directly
/Users/wuhhh/Desktop/徐老师/代码/metagenomic/.venv_web/bin/python -m bac_analysis_portal.app
```

### Run the desktop app (local Flask + pywebview)
```bash
python run_bac_analysis_desktop.py
```

### Run the Neisseria meningitidis delivery web UI
```bash
python delivery/neisseria_meningitidis_pipeline/web_app.py
# Default address: http://127.0.0.1:5088
```

### Run tests
```bash
pytest tests/
# Run a single test file
pytest tests/test_workflow_flow_split.py
# Run a single test
pytest tests/test_workflow_flow_split.py -k test_get_flow_list_expands_legacy_mlst_and_serotype_flow
```

### Build the macOS desktop app
```bash
bash scripts/build_mac_desktop_app.sh
# Build and open immediately
bash scripts/build_mac_desktop_app.sh --open
```

## High-Level Architecture

### Web Portal (`bac_analysis_portal/`)

- `app.py` — Large Flask application exposing routes for login, sample library management, task submission, monitoring, result reports, and knowledge base integration.
- `store.py` — `PortalStore` handles all SQLite schema and queries (users, sample library, submissions, version logs, audit logs, host databases, reference panels).
- `task_manager.py` — `AnalysisTaskManager` validates inputs, builds command-line arguments, and writes JSON task descriptors to `analysis_tasks/`.
- `task_runner.py` — Reads task JSON, spawns the actual pipeline as a subprocess, streams stdout/stderr to a log file, and updates task status. Includes auto-pause logic for species review and Neisseria AMR post-processing.
- `sample_library_manager.py` — Bulk import and metadata management for sample libraries.
- `knowledge_base.py` — Loads knowledge base bundles (pathogens, taxonomies, MLST/serotype associations) from `database/knowledge_base/`.

The portal communicates with pipelines through JSON task files and log files on disk, not via HTTP.

### Refactored Pipeline (`metagenomic_refactor/`)

- `context.py` — Global `RuntimeContext` dataclass used to share state (output dir, threads, species, resource paths) across modules. Initialized via `set_runtime_context()`.
- `config.py` — `PipelineConfig` and `ResourcePaths` dataclasses. Resource paths default to hard-coded server paths but can be overridden via `META_*` environment variables.
- `workflow.py` — High-level workflow orchestration. Translates the `runflow` string into ordered steps, handles influenza/SARS-CoV-2 special cases, and drives QC -> taxonomy -> assembly -> annotation/report.
- `runner.py` — Input validation, directory setup, and sample loop that bridges the web task format to `workflow.py`.
- `qc.py`, `assembly.py`, `taxonomy.py`, `strain_typing.py`, `virus_analysis.py`, `annotation.py`, `report.py` — Module-level implementations of each analysis stage.

### Pathogen Source Tracking (`pathosource_refactor/`)

- `workflow.py` — `main_process()` that runs snippy -> consensus -> tree building -> pairwise distance -> MLST/cgMLST.
- `phylogeny.py` — Snippy core SNP extraction, Gubbins recombination filtering, RAxML/IQ-TREE tree building, and pairwise distance matrices.
- `typing.py` — MLST and cgMLST wrappers.

### Neisseria Delivery Pipeline (`delivery/neisseria_meningitidis_pipeline/`)

A standalone copy of the main pipeline customized for *Neisseria meningitidis*. Entry points:
- `run_neisseria_meningitidis_pipeline.py` — CLI pipeline runner (does not depend on `metagenomic_refactor`).
- `web_app.py` — Small Flask app for local task submission and result review.
- `nm_pipeline/` — Copied pipeline modules (QC, assembly, typing, AMR, report) local to this delivery.

## Design Conventions

- UI style is light, clinical, and information-dense. Avoid decorative gradients or dark neon aesthetics.
- The portal supports both bacteria and virus analysis targets, plus a metagenome mode when `method == "meta"`.
- Progress is communicated from pipeline to portal by printing lines matching:
  ```
  task_step：{step}/{total_step}\t样本进度：{sample_index}/{sample_total}\t样本：{sample}\t{message}
  ```
- Task files are JSON documents stored under `analysis_tasks/{timestamp}_{task_id}/task.json`. Logs go to `pipeline.log` in the same directory.

## Environment and Dependencies

- Python 3.10+.
- Web portal virtual environment: `.venv_web/` (Flask, SQLAlchemy, pywebview).
- Test dependencies: `pytest`, `SQLAlchemy` (see `requirements-dev.txt`).
- Many bioinformatics tools are assumed to be installed via conda and available in PATH or invoked via `conda run -n <env>` (e.g., `VFind`, `BASALT`). See `requirements-viral-assembly.txt` for the virus-assembly node requirements.

## Important File Paths

- `bac_analysis_portal.sqlite3` — Production SQLite DB for the portal.
- `genome_db.sqlite3` — Genome database mini-app DB.
- `database/knowledge_base/` — JSON knowledge base files for pathogen taxonomy, MLST, and serotype associations.
- `analysis_tasks/` — All submitted tasks and their outputs.
