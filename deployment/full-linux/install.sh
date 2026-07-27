#!/usr/bin/env bash
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="/opt/pathogen-workbench"
DATA_ROOT="/data/pathogen-workbench"
DB_ROOT="/data/pathogen-db"
PORT="5055"
SERVER_NAME="_"
YES="0"
FORCE_CONDA="0"
SERVICE_NAME="bac-analysis-portal"
SERVICE_USER="pathogen-workbench"
SERVICE_GROUP="pathogen-workbench"
REQUIRED_ENVS="web_runtime,meta_main,genomad_aux,amr_aux,host_filter,mag_aux,cm210,sistr_hicap,qiime2,microeco,genovi,plasflow,longread_aux,chewie,ncov,choleraefinder,report_env,Vlib,tNGS"

usage() {
  cat <<'USAGE'
Usage:
  sudo bash install.sh [options]

Options:
  --prefix PATH       Installation root. Default: /opt/pathogen-workbench
  --data-root PATH    Portal state/task root. Default: /data/pathogen-workbench
  --db-root PATH      Database install root. Default: /data/pathogen-db
  --port PORT         Portal HTTP port behind nginx. Default: 5055
  --server-name NAME  nginx server_name. Default: _
  --yes              Do not prompt before installing
  --force-conda      Replace existing unpacked Conda environments
  -h, --help          Show this help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="${2:?}"; shift 2 ;;
    --data-root) DATA_ROOT="${2:?}"; shift 2 ;;
    --db-root) DB_ROOT="${2:?}"; shift 2 ;;
    --port) PORT="${2:?}"; shift 2 ;;
    --server-name) SERVER_NAME="${2:?}"; shift 2 ;;
    --yes) YES="1"; shift ;;
    --force-conda) FORCE_CONDA="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

APP_DIR="${PREFIX}/app"
CONDA_ROOT="${PREFIX}/conda"
CONDA_ENVS_DIR="${CONDA_ROOT}/envs"
STATE_DIR="${DATA_ROOT}/state"
TASK_ROOT="${DATA_ROOT}/tasks"
OUTPUT_ROOT="${DATA_ROOT}/outputs"
ENV_FILE="${APP_DIR}/deployment/portal.env"
MANIFEST_FILE="${BUNDLE_DIR}/manifest.json"
ADMIN_PASSWORD_FILE="${STATE_DIR}/initial_admin_password.txt"

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

sha256_of() {
  sha256sum "$1" | awk '{print $1}'
}

random_token() {
  python3 - <<'PY'
import secrets
print(secrets.token_urlsafe(48))
PY
}

verify_manifest() {
  echo "Verifying manifest checksums..."
  python3 - "$BUNDLE_DIR" "$MANIFEST_FILE" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
manifest = Path(sys.argv[2])
if not manifest.is_file():
    raise SystemExit(f"manifest not found: {manifest}")
data = json.loads(manifest.read_text(encoding="utf-8"))
errors = []
for item in data.get("files", []):
    rel = item.get("path", "")
    expected = item.get("sha256", "")
    target = (root / rel).resolve()
    if not str(target).startswith(str(root)):
        errors.append(f"unsafe manifest path: {rel}")
        continue
    if not target.is_file():
        errors.append(f"missing file: {rel}")
        continue
    digest = hashlib.sha256()
    with target.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    if actual != expected:
        errors.append(f"sha256 mismatch: {rel}")
if errors:
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    raise SystemExit(1)
print(f"Verified {len(data.get('files', []))} files.")
PY
}

load_manifest_required_envs() {
  python3 - "$MANIFEST_FILE" <<'PY'
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
if not manifest.is_file():
    raise SystemExit(0)
data = json.loads(manifest.read_text(encoding="utf-8"))
envs = [str(item).strip() for item in data.get("required_conda_envs", []) if str(item).strip()]
if envs:
    print(",".join(envs))
PY
}

extract_pack() {
  local pack="$1"
  local dest="$2"
  mkdir -p "$dest"
  case "$pack" in
    *.tar.zst|*.tzst) tar -I zstd -xf "$pack" -C "$dest" ;;
    *.tar.gz|*.tgz) tar -xzf "$pack" -C "$dest" ;;
    *.tar.xz|*.txz) tar -xJf "$pack" -C "$dest" ;;
    *) echo "ERROR: unsupported conda pack format: $pack" >&2; exit 2 ;;
  esac
}

install_conda_env() {
  local env_name="$1"
  local pack="${BUNDLE_DIR}/conda-packs/${env_name}.tar.zst"
  local dest="${CONDA_ENVS_DIR}/${env_name}"
  local stamp="${dest}/.pathogen-pack.sha256"

  if [[ ! -f "$pack" ]]; then
    echo "ERROR: missing Conda pack: $pack" >&2
    exit 2
  fi

  local digest
  digest="$(sha256_of "$pack")"
  if [[ -d "$dest" && "$FORCE_CONDA" != "1" && -f "$stamp" && "$(cat "$stamp")" == "$digest" ]]; then
    echo "Conda env already installed: $env_name"
    return
  fi
  if [[ -d "$dest" ]]; then
    local backup="${dest}.bak.$(date +%Y%m%d_%H%M%S)"
    echo "Backing up existing Conda env: $dest -> $backup"
    mv "$dest" "$backup"
  fi
  mkdir -p "$dest"
  echo "Extracting Conda env: $env_name"
  extract_pack "$pack" "$dest"
  if [[ -x "${dest}/bin/conda-unpack" ]]; then
    "${dest}/bin/conda-unpack"
  else
    echo "WARN: conda-unpack not found in ${env_name}; prefix relocation may be incomplete." >&2
  fi
  echo "$digest" > "$stamp"
}

if [[ "$(id -u)" -ne 0 ]]; then
  echo "ERROR: run as root, for example: sudo bash install.sh --yes" >&2
  exit 2
fi

require_cmd python3
require_cmd tar
require_cmd sha256sum
require_cmd systemctl
require_cmd zstd

if [[ ! -f /etc/os-release ]] || ! grep -qi '^ID=ubuntu' /etc/os-release; then
  echo "ERROR: this installer currently supports Ubuntu only." >&2
  exit 2
fi
if ! command -v nginx >/dev/null 2>&1; then
  echo "ERROR: nginx is not installed. Install nginx before running this offline installer." >&2
  exit 2
fi

for required_dir in app conda-packs database soft public; do
  if [[ ! -e "${BUNDLE_DIR}/${required_dir}" ]]; then
    echo "ERROR: release bundle is missing ${required_dir}/" >&2
    exit 2
  fi
done
if [[ ! -d "${BUNDLE_DIR}/conda-runtime" || ! -x "${BUNDLE_DIR}/conda-runtime/bin/conda" ]]; then
  echo "ERROR: release bundle must include conda-runtime/bin/conda." >&2
  exit 2
fi

if [[ "$YES" != "1" ]]; then
  cat <<EOF
Install Pathogen Workbench full offline bundle:
  app        ${APP_DIR}
  conda      ${CONDA_ROOT}
  state      ${STATE_DIR}
  tasks      ${TASK_ROOT}
  database   ${DB_ROOT}
  port       ${PORT}

EOF
  read -r -p "Continue? [y/N] " answer
  case "$answer" in
    y|Y|yes|YES) ;;
    *) echo "Aborted."; exit 1 ;;
  esac
fi

verify_manifest
manifest_envs="$(load_manifest_required_envs || true)"
if [[ -n "$manifest_envs" ]]; then
  REQUIRED_ENVS="$manifest_envs"
fi

echo "Creating install directories..."
mkdir -p "$PREFIX" "$APP_DIR" "$CONDA_ROOT" "$CONDA_ENVS_DIR" "$STATE_DIR" "$TASK_ROOT" "$OUTPUT_ROOT" "$DB_ROOT"
if ! getent group "$SERVICE_GROUP" >/dev/null; then
  groupadd --system "$SERVICE_GROUP"
fi
if ! id "$SERVICE_USER" >/dev/null 2>&1; then
  useradd --system --gid "$SERVICE_GROUP" --home-dir /nonexistent --shell /usr/sbin/nologin "$SERVICE_USER"
fi
chown "$SERVICE_USER:$SERVICE_GROUP" "$STATE_DIR" "$TASK_ROOT" "$OUTPUT_ROOT"
chmod 0750 "$STATE_DIR" "$TASK_ROOT" "$OUTPUT_ROOT"

echo "Installing application files..."
copy_tree "${BUNDLE_DIR}/app" "$APP_DIR"

echo "Installing Conda runtime..."
copy_tree "${BUNDLE_DIR}/conda-runtime" "$CONDA_ROOT"

IFS=',' read -r -a env_names <<< "$REQUIRED_ENVS"
for env_name in "${env_names[@]}"; do
  install_conda_env "$env_name"
done

echo "Installing database and runtime assets..."
copy_tree "${BUNDLE_DIR}/database" "$DB_ROOT"
copy_tree "${BUNDLE_DIR}/soft" "${APP_DIR}/soft"
copy_tree "${BUNDLE_DIR}/public" "${APP_DIR}/public"
if [[ -d "${BUNDLE_DIR}/smoke-tests" ]]; then
  copy_tree "${BUNDLE_DIR}/smoke-tests" "${APP_DIR}/smoke-tests"
fi

SECRET_KEY="$(random_token)"
INITIAL_ADMIN_PASSWORD="$(random_token)"
if [[ -f "$ENV_FILE" ]]; then
  # Preserve existing production secrets on re-install.
  # shellcheck disable=SC1090
  source "$ENV_FILE" || true
  SECRET_KEY="${PORTAL_SECRET_KEY:-$SECRET_KEY}"
  if [[ -f "$ADMIN_PASSWORD_FILE" ]]; then
    INITIAL_ADMIN_PASSWORD="$(cat "$ADMIN_PASSWORD_FILE")"
  fi
fi
mkdir -p "$(dirname "$ENV_FILE")"
cat > "$ENV_FILE" <<EOF
PORTAL_MODE=production
PORTAL_SECRET_KEY=${SECRET_KEY}
PORTAL_INITIAL_ADMIN_PASSWORD=${INITIAL_ADMIN_PASSWORD}
PORTAL_COOKIE_SECURE=0
PORTAL_STATE_DIR=${STATE_DIR}
PORTAL_TASK_DIR=${TASK_ROOT}
PORTAL_DB_PATH=${STATE_DIR}/bac_analysis_portal.sqlite3
BAC_ANALYSIS_TASK_ROOT=${TASK_ROOT}
BAC_ANALYSIS_OUTPUT_ROOT=${OUTPUT_ROOT}
META_DATABASE_ROOT=${DB_ROOT}
CONDA_ROOT=${CONDA_ROOT}
PORTAL_CONDA_ROOT=${CONDA_ROOT}
CONDA_EXE=${CONDA_ROOT}/bin/conda
PORTAL_REQUIRED_CONDA_ENVS=${REQUIRED_ENVS}
APP_HOST=127.0.0.1
APP_PORT=${PORT}
GUNICORN_WORKERS=2
SERVICE_NAME=${SERVICE_NAME}
SERVER_NAME=${SERVER_NAME}
META_KRAKEN_DB=${DB_ROOT}/kraken2
META_VIRUS_KRAKEN_DB=${DB_ROOT}/virus_kraken2
META_VIRSORTER2_DB=${DB_ROOT}/virsorter2
META_CHECKV_DB=${DB_ROOT}/checkv-db
META_GENOMAD_DB=${DB_ROOT}/genomad
META_MOBILEOG_DB=${DB_ROOT}/mobileOG-db
META_MOBILEOG_META=${DB_ROOT}/mobileOG-db-beatrix.csv
EOF
chmod 600 "$ENV_FILE"
echo "$INITIAL_ADMIN_PASSWORD" > "$ADMIN_PASSWORD_FILE"
chmod 600 "$ADMIN_PASSWORD_FILE"

echo "Seeding portal runtime settings..."
set -a
# shellcheck disable=SC1090
source "$ENV_FILE"
set +a
"${CONDA_ENVS_DIR}/web_runtime/bin/python" - <<PY
from pathlib import Path

from bac_analysis_portal.store import PortalStore

project_root = Path(${APP_DIR@Q})
store = PortalStore.from_project_root(
    project_root,
    initial_admin_password=${INITIAL_ADMIN_PASSWORD@Q},
    rotate_weak_admin=True,
)
settings = {
    "workspace_root": str(project_root),
    "pipeline_script": "Bac_assemble_260112_newformat.py",
    "pipeline_python": "meta_main",
    "conda_root": ${CONDA_ROOT@Q},
    "database_root": ${DB_ROOT@Q},
    "max_concurrent_tasks": "2",
}
conda_envs = {
    "vfind": "genomad_aux",
    "hamronization": "amr_aux",
    "rgi": "amr_aux",
    "rgi_new": "amr_aux",
    "hostile": "host_filter",
    "kneaddata": "host_filter",
    "basalt": "mag_aux",
    "coverm": "mag_aux",
    "cm210": "mag_aux",
    "gtdbtk": "mag_aux",
    "sistr_hicap": "sistr_hicap",
    "qiime2": "qiime2",
    "microeco": "microeco",
    "genovi": "genovi",
    "plasflow": "plasflow",
    "medaka": "longread_aux",
    "clair3": "longread_aux",
    "chewie": "chewie",
    "tb_profiler": "ncov",
    "choleraefinder": "choleraefinder",
    "report_env": "report_env",
    "vlib": "Vlib",
    "tngs": "tNGS",
}
for key, value in settings.items():
    store.set_setting(key, value)
for key, value in conda_envs.items():
    store.set_setting(f"conda_env_{key}", value)
PY

echo "Writing systemd service..."
cat > "/etc/systemd/system/${SERVICE_NAME}.service" <<EOF
[Unit]
Description=Pathogen Workbench Portal
After=network.target

[Service]
User=${SERVICE_USER}
Group=${SERVICE_GROUP}
WorkingDirectory=${APP_DIR}
EnvironmentFile=${ENV_FILE}
Environment=PYTHONUNBUFFERED=1
ExecStart=${CONDA_ENVS_DIR}/web_runtime/bin/gunicorn -w 2 -b 127.0.0.1:${PORT} 'bac_analysis_portal:create_app()'
Restart=always
RestartSec=5
UMask=0077
NoNewPrivileges=true
PrivateTmp=true
ProtectHome=true
ProtectSystem=full
ReadOnlyPaths=${APP_DIR} ${DB_ROOT} ${CONDA_ROOT}
ReadWritePaths=${STATE_DIR} ${TASK_ROOT} ${OUTPUT_ROOT}
CapabilityBoundingSet=

[Install]
WantedBy=multi-user.target
EOF

echo "Writing nginx site..."
cat > "/etc/nginx/sites-available/${SERVICE_NAME}" <<EOF
server {
    listen 80;
    server_name ${SERVER_NAME};
    client_max_body_size 2G;

    location / {
        proxy_pass http://127.0.0.1:${PORT};
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
    }
}
EOF
ln -sf "/etc/nginx/sites-available/${SERVICE_NAME}" "/etc/nginx/sites-enabled/${SERVICE_NAME}"

echo "Checking deployment..."
"${CONDA_ENVS_DIR}/web_runtime/bin/python" "${APP_DIR}/scripts/check_deployment.py" \
  --project-root "$APP_DIR" \
  --env-file "$ENV_FILE" \
  --profile full-analysis

systemctl daemon-reload
systemctl enable --now "$SERVICE_NAME"
nginx -t
systemctl reload nginx

echo
echo "Pathogen Workbench installed."
echo "URL: http://$(hostname -I | awk '{print $1}')/login"
echo "Initial admin username: admin"
echo "Initial admin password: ${INITIAL_ADMIN_PASSWORD}"
echo "Password saved once at: ${ADMIN_PASSWORD_FILE}"
