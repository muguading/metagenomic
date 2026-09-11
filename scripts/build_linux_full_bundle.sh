#!/usr/bin/env bash
set -Eeuo pipefail

VERSION=""
SOURCE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_ROOT=""
DATABASE_ROOT=""
SOFT_ROOT=""
PUBLIC_ROOT=""
SMOKE_TESTS_ROOT=""
OUTPUT=""
ARCHIVE="0"
REQUIRED_ENVS="web_runtime,meta_main,genomad_aux,amr_aux,host_filter,mag_aux,cm210,sistr_hicap,qiime2,microeco,genovi,plasflow,longread_aux,chewie,ncov,choleraefinder,report_env,Vlib,tNGS"
MIN_FREE_GB="80"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/build_linux_full_bundle.sh \
    --version 2026.06 \
    --source-root /path/to/metagenomic \
    --conda-root /opt/miniconda3 \
    --database-root /data/pathogen-db \
    --output /mnt/usb/pathogen-workbench-full-linux

Options:
  --version VERSION          Release version label. Required.
  --source-root PATH         Source repository root. Default: current repo.
  --conda-root PATH          Prepared Conda root containing bin/conda and envs/. Required.
  --database-root PATH       Full database asset root. Required.
  --soft-root PATH           Tool/model asset root. Default: SOURCE_ROOT/soft.
  --public-root PATH         Public visual asset root. Default: SOURCE_ROOT/public.
  --smoke-tests-root PATH    Optional smoke test directory to include.
  --output PATH              Output release directory. Required.
  --required-envs CSV        Conda env names to pack.
  --min-free-gb N            Minimum target free-space hint in manifest. Default: 80.
  --archive                  Also create OUTPUT-<version>.tar.zst next to output dir.
  -h, --help                 Show this help.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --version) VERSION="${2:?}"; shift 2 ;;
    --source-root) SOURCE_ROOT="${2:?}"; shift 2 ;;
    --conda-root) CONDA_ROOT="${2:?}"; shift 2 ;;
    --database-root) DATABASE_ROOT="${2:?}"; shift 2 ;;
    --soft-root) SOFT_ROOT="${2:?}"; shift 2 ;;
    --public-root) PUBLIC_ROOT="${2:?}"; shift 2 ;;
    --smoke-tests-root) SMOKE_TESTS_ROOT="${2:?}"; shift 2 ;;
    --output) OUTPUT="${2:?}"; shift 2 ;;
    --required-envs) REQUIRED_ENVS="${2:?}"; shift 2 ;;
    --min-free-gb) MIN_FREE_GB="${2:?}"; shift 2 ;;
    --archive) ARCHIVE="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

SOURCE_ROOT="$(cd "$SOURCE_ROOT" && pwd)"
SOFT_ROOT="${SOFT_ROOT:-${SOURCE_ROOT}/soft}"
PUBLIC_ROOT="${PUBLIC_ROOT:-${SOURCE_ROOT}/public}"

if [[ -z "$VERSION" || -z "$CONDA_ROOT" || -z "$DATABASE_ROOT" || -z "$OUTPUT" ]]; then
  usage >&2
  exit 2
fi

CONDA_ROOT="$(cd "$CONDA_ROOT" && pwd)"
DATABASE_ROOT="$(cd "$DATABASE_ROOT" && pwd)"
SOFT_ROOT="$(cd "$SOFT_ROOT" && pwd)"
PUBLIC_ROOT="$(cd "$PUBLIC_ROOT" && pwd)"
OUTPUT_PARENT="$(dirname "$OUTPUT")"
mkdir -p "$OUTPUT_PARENT"
OUTPUT="$(cd "$OUTPUT_PARENT" && pwd)/$(basename "$OUTPUT")"

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: required command not found: $1" >&2
    exit 2
  fi
}

copy_tree() {
  local src="$1"
  local dst="$2"
  mkdir -p "$dst"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a "$src"/ "$dst"/
  else
    cp -a "$src"/. "$dst"/
  fi
}

copy_app_source() {
  local src="$1"
  local dst="$2"
  mkdir -p "$dst"
  local entries=(
    "bac_analysis_portal"
    "metagenomic_refactor"
    "pathosource_refactor"
    "virus_pipeline_skill"
    "scripts"
    "docs"
    "envs"
    "Bac_assemble_260112_newformat.py"
    "Virus_WTSpip.py"
    "CommunityAnalysis.py"
    "PathoSource.py"
    "run_metagenome_analysis.sh"
    "run_bac_analysis_desktop.py"
    "requirements-web.txt"
    "requirements-release.lock"
    "requirements-dev.txt"
    "requirements-viral-assembly.txt"
    "pytest.ini"
    "README.md"
    "README_desktop_app.md"
    "README_refactor.md"
    "README_cgmlst_download.md"
    "deploy_bac_analysis_portal_ubuntu.sh"
    "HP_0601logo.png"
    "logo_transparent_bg.png"
  )
  for entry in "${entries[@]}"; do
    if [[ ! -e "${src}/${entry}" ]]; then
      continue
    fi
    if [[ -d "${src}/${entry}" ]]; then
      mkdir -p "${dst}/${entry}"
      if command -v rsync >/dev/null 2>&1; then
        rsync -a \
          --exclude='__pycache__/' \
          --exclude='.pytest_cache/' \
          --exclude='.ruff_cache/' \
          --exclude='*.pyc' \
          --exclude='*.zip' \
          --exclude='datasets_linux/' \
          --exclude='datasets_macos/' \
          "${src}/${entry}"/ "${dst}/${entry}"/
      else
        copy_tree "${src}/${entry}" "${dst}/${entry}"
      fi
    else
      cp -p "${src}/${entry}" "${dst}/${entry}"
    fi
  done
  mkdir -p "${dst}/deployment"
  if [[ -d "${src}/deployment" ]]; then
    copy_tree "${src}/deployment" "${dst}/deployment"
  fi
}

require_cmd python3
require_cmd tar
require_cmd sha256sum
require_cmd zstd
require_cmd conda-pack

if [[ ! -x "${CONDA_ROOT}/bin/conda" ]]; then
  echo "ERROR: --conda-root must contain bin/conda: ${CONDA_ROOT}" >&2
  exit 2
fi

IFS=',' read -r -a env_names <<< "$REQUIRED_ENVS"
for env_name in "${env_names[@]}"; do
  if [[ ! -d "${CONDA_ROOT}/envs/${env_name}" ]]; then
    echo "ERROR: required Conda environment not found: ${CONDA_ROOT}/envs/${env_name}" >&2
    exit 2
  fi
done

rm -rf "$OUTPUT"
mkdir -p "$OUTPUT"/{app,conda-runtime,conda-packs,database,soft,public,smoke-tests}

echo "Copying application source..."
copy_app_source "$SOURCE_ROOT" "$OUTPUT/app"

echo "Adding full offline installer..."
install -m 0755 "$SOURCE_ROOT/deployment/full-linux/install.sh" "$OUTPUT/install.sh"

echo "Copying Conda runtime..."
if command -v rsync >/dev/null 2>&1; then
  rsync -a \
    --exclude='envs/' \
    --exclude='pkgs/' \
    --exclude='conda-meta/history' \
    "$CONDA_ROOT"/ "$OUTPUT/conda-runtime"/
else
  copy_tree "$CONDA_ROOT" "$OUTPUT/conda-runtime"
  rm -rf "$OUTPUT/conda-runtime/envs" "$OUTPUT/conda-runtime/pkgs"
fi

pack_env() {
  local env_name="$1"
  local env_dir="${CONDA_ROOT}/envs/${env_name}"
  local target="${OUTPUT}/conda-packs/${env_name}.tar.zst"
  local tmp_dir
  tmp_dir="$(mktemp -d)"
  local tmp_gz="${tmp_dir}/${env_name}.tar.gz"

  echo "Packing Conda env: ${env_name}"
  conda-pack -p "$env_dir" -o "$tmp_gz" --force
  gzip -dc "$tmp_gz" | zstd -T0 -19 -q -o "$target"
  rm -rf "$tmp_dir"
}

for env_name in "${env_names[@]}"; do
  pack_env "$env_name"
done

echo "Copying database assets..."
copy_tree "$DATABASE_ROOT" "$OUTPUT/database"

echo "Copying tool and public assets..."
copy_tree "$SOFT_ROOT" "$OUTPUT/soft"
copy_tree "$PUBLIC_ROOT" "$OUTPUT/public"

if [[ -n "$SMOKE_TESTS_ROOT" ]]; then
  SMOKE_TESTS_ROOT="$(cd "$SMOKE_TESTS_ROOT" && pwd)"
  copy_tree "$SMOKE_TESTS_ROOT" "$OUTPUT/smoke-tests"
fi

cat > "$OUTPUT/README_INSTALL.md" <<EOF
# Pathogen Workbench Full Linux Bundle

Version: ${VERSION}

Install on the target Ubuntu server:

\`\`\`bash
sudo bash $(basename "$OUTPUT")/install.sh --yes
\`\`\`

This bundle is designed for offline installation. Conda environments are already packed under \`conda-packs/\`.
EOF

echo "Writing manifest..."
python3 - "$OUTPUT" "$VERSION" "$REQUIRED_ENVS" "$MIN_FREE_GB" <<'PY'
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

root = Path(sys.argv[1]).resolve()
version = sys.argv[2]
required_envs = [item.strip() for item in sys.argv[3].split(",") if item.strip()]
min_free_gb = float(sys.argv[4])

files = []
for path in sorted(root.rglob("*")):
    if not path.is_file() or path.name == "manifest.json":
        continue
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    files.append({
        "path": path.relative_to(root).as_posix(),
        "size": size,
        "sha256": digest.hexdigest(),
    })

manifest = {
    "name": "pathogen-workbench-full-linux",
    "version": version,
    "generated_at": datetime.now(timezone.utc).isoformat(),
    "required_platform": {"os": "ubuntu", "arch": "linux-x86_64"},
    "required_conda_envs": required_envs,
    "min_free_gb": min_free_gb,
    "install_defaults": {
        "prefix": "/opt/pathogen-workbench",
        "data_root": "/data/pathogen-workbench",
        "db_root": "/data/pathogen-db",
        "port": 5055,
    },
    "files": files,
}
(root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY

if [[ "$ARCHIVE" == "1" ]]; then
  archive_path="${OUTPUT_PARENT}/$(basename "$OUTPUT")-${VERSION}.tar.zst"
  echo "Creating archive: ${archive_path}"
  tar -I 'zstd -T0 -19' -cf "$archive_path" -C "$OUTPUT_PARENT" "$(basename "$OUTPUT")"
fi

echo
echo "Full Linux bundle ready: $OUTPUT"
echo "Required envs: $REQUIRED_ENVS"
