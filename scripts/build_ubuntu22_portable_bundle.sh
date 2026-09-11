#!/usr/bin/env bash
set -Eeuo pipefail

SOURCE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT="/Volumes/10gelizi/pathogen-workbench-ubuntu22-portable"
VERSION="$(date +%Y.%m.%d)"
ARCHIVE="0"
INCLUDE_DEMO="0"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/build_ubuntu22_portable_bundle.sh [options]

Options:
  --source-root PATH  Source repository root. Default: current repo.
  --output PATH       Output bundle directory.
                      Default: /Volumes/10gelizi/pathogen-workbench-ubuntu22-portable
  --version VERSION   Version label. Default: current date.
  --archive           Also create OUTPUT-VERSION.tar.gz next to the output dir.
  --include-demo      Include the large demo dataset. Disabled by default.
  -h, --help          Show this help.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-root) SOURCE_ROOT="${2:?}"; shift 2 ;;
    --output) OUTPUT="${2:?}"; shift 2 ;;
    --version) VERSION="${2:?}"; shift 2 ;;
    --archive) ARCHIVE="1"; shift ;;
    --include-demo) INCLUDE_DEMO="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

SOURCE_ROOT="$(cd "$SOURCE_ROOT" && pwd)"
OUTPUT_PARENT="$(dirname "$OUTPUT")"
mkdir -p "$OUTPUT_PARENT"
OUTPUT="$(cd "$OUTPUT_PARENT" && pwd)/$(basename "$OUTPUT")"

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
    "deployment"
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
    "virus_detect.xlsx"
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
}

for required in database soft public requirements-web.txt requirements-release.lock bac_analysis_portal scripts/check_deployment.py; do
  if [[ ! -e "${SOURCE_ROOT}/${required}" ]]; then
    echo "ERROR: missing required project asset: ${SOURCE_ROOT}/${required}" >&2
    exit 2
  fi
done

rm -rf "$OUTPUT"
mkdir -p "$OUTPUT"/{app,database,soft,public}

echo "Copying application source..."
copy_app_source "$SOURCE_ROOT" "$OUTPUT/app"

echo "Adding Ubuntu 22.04 installer..."
install -m 0755 "$SOURCE_ROOT/deployment/ubuntu22-portable/install.sh" "$OUTPUT/install.sh"
cat > "$OUTPUT/deploy.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec sudo bash "${BUNDLE_DIR}/install.sh" --yes "$@"
EOF
chmod 0755 "$OUTPUT/deploy.sh"

echo "Copying database assets..."
copy_tree "$SOURCE_ROOT/database" "$OUTPUT/database"

echo "Copying tool and public assets..."
copy_tree "$SOURCE_ROOT/soft" "$OUTPUT/soft"
copy_tree "$SOURCE_ROOT/public" "$OUTPUT/public"

if [[ "$INCLUDE_DEMO" == "1" && -d "$SOURCE_ROOT/demo_data" ]]; then
  echo "Copying demo data..."
  copy_tree "$SOURCE_ROOT/demo_data" "$OUTPUT/demo_data"
fi

cat > "$OUTPUT/README_UBUNTU22_INSTALL.md" <<EOF
# Pathogen Workbench Ubuntu 22.04 Portable Bundle

Version: ${VERSION}

This bundle contains the application source, database assets, soft assets, and
public assets copied from:

\`\`\`text
${SOURCE_ROOT}
\`\`\`

## One-command install on a fresh Ubuntu 22.04 server

Mount or copy this directory to the server, then run:

\`\`\`bash
cd /path/to/pathogen-workbench-ubuntu22-portable
sudo bash install.sh --yes
\`\`\`

Custom install example:

\`\`\`bash
sudo bash install.sh \\
  --prefix /opt/pathogen-workbench \\
  --data-root /data/pathogen-workbench \\
  --db-root /data/pathogen-db \\
  --port 5055 \\
  --server-name _ \\
  --yes
\`\`\`

After installation, open:

\`\`\`text
http://SERVER_IP/login
\`\`\`

The installer prints the generated initial \`admin\` password and saves it to:

\`\`\`text
/data/pathogen-workbench/state/initial_admin_password.txt
\`\`\`

## Important boundary

This bundle installs the Web Portal and copies database/runtime assets. Full
sequencing analysis still requires Linux Conda environments under
\`/opt/miniconda3/envs/\`. A macOS Conda environment cannot be copied to Ubuntu.

After preparing those Linux environments, validate the full analysis profile:

\`\`\`bash
/opt/pathogen-workbench/.venv_web/bin/python \\
  /opt/pathogen-workbench/app/scripts/check_deployment.py \\
  --project-root /opt/pathogen-workbench/app \\
  --env-file /opt/pathogen-workbench/app/deployment/portal.env \\
  --profile full-analysis
\`\`\`
EOF

echo "Writing checksum manifest..."
(
  cd "$OUTPUT"
  find . -type f \
    ! -name manifest.sha256 \
    ! -name '._*' \
    -print0 | sort -z | xargs -0 shasum -a 256 > manifest.sha256
)

if [[ "$ARCHIVE" == "1" ]]; then
  archive_path="${OUTPUT_PARENT}/$(basename "$OUTPUT")-${VERSION}.tar.gz"
  echo "Creating archive: ${archive_path}"
  tar -czf "$archive_path" -C "$OUTPUT_PARENT" "$(basename "$OUTPUT")"
fi

echo
echo "Ubuntu 22.04 portable bundle ready: $OUTPUT"
du -sh "$OUTPUT"
