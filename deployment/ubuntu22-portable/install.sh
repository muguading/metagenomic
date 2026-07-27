#!/usr/bin/env bash
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="/opt/pathogen-workbench"
DATA_ROOT="/data/pathogen-workbench"
DB_ROOT="/data/pathogen-db"
PORT="5055"
SERVER_NAME="_"
SERVICE_NAME="bac-analysis-portal"
YES="0"
INSTALL_APT="1"
DEFAULT_REQUIRED_ENVS="web_runtime,meta_main,genomad_aux,amr_aux,host_filter,mag_aux,cm210,sistr_hicap,qiime2,microeco,genovi,plasflow,longread_aux,chewie,ncov,choleraefinder,report_env,Vlib,tNGS"

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
  --no-apt            Do not install Ubuntu packages
  --yes              Do not prompt before installing
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
    --no-apt) INSTALL_APT="0"; shift ;;
    --yes) YES="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

APP_DIR="${PREFIX}/app"
VENV_DIR="${PREFIX}/.venv_web"
PACKED_CONDA_ROOT="${PREFIX}/conda"
PACKED_CONDA_ENVS_DIR="${PACKED_CONDA_ROOT}/envs"
STATE_DIR="${DATA_ROOT}/state"
TASK_ROOT="${DATA_ROOT}/tasks"
ENV_FILE="${APP_DIR}/deployment/portal.env"
ADMIN_PASSWORD_FILE="${STATE_DIR}/initial_admin_password.txt"
MANIFEST_FILE="${BUNDLE_DIR}/manifest.sha256"
if [[ -d "${BUNDLE_DIR}/conda-runtime" || -f "${BUNDLE_DIR}/conda-runtime.tar.zst" || -d "${BUNDLE_DIR}/conda-packs" ]]; then
  ACTIVE_CONDA_ROOT="$PACKED_CONDA_ROOT"
else
  ACTIVE_CONDA_ROOT="/opt/miniconda3"
fi

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
    rsync -a --delete "$src"/ "$dst"/
  else
    cp -a "$src"/. "$dst"/
  fi
}

random_token() {
  python3 - <<'PY'
import secrets
print(secrets.token_urlsafe(48))
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
    *) echo "ERROR: unsupported Conda pack format: $pack" >&2; exit 2 ;;
  esac
}

if [[ "$(id -u)" -ne 0 ]]; then
  echo "ERROR: run as root, for example: sudo bash install.sh --yes" >&2
  exit 2
fi

if [[ ! -f /etc/os-release ]] || ! grep -qi '^ID=ubuntu' /etc/os-release; then
  echo "ERROR: this installer supports Ubuntu only." >&2
  exit 2
fi
if ! grep -q '^VERSION_ID="22.04"' /etc/os-release && ! grep -q '^VERSION_ID=22.04' /etc/os-release; then
  echo "WARN: this bundle was written for Ubuntu 22.04; continuing anyway." >&2
fi

for required_dir in app database soft public; do
  if [[ ! -e "${BUNDLE_DIR}/${required_dir}" ]]; then
    echo "ERROR: bundle is missing ${required_dir}/" >&2
    exit 2
  fi
done

if [[ "$INSTALL_APT" == "1" ]]; then
  apt-get update
  DEBIAN_FRONTEND=noninteractive apt-get install -y \
    ca-certificates curl nginx python3 python3-pip python3-venv rsync zstd
fi

require_cmd python3
require_cmd nginx
if [[ -d "${BUNDLE_DIR}/conda-packs" || -f "${BUNDLE_DIR}/conda-runtime.tar.zst" ]]; then
  require_cmd zstd
fi

if [[ -f "$MANIFEST_FILE" ]]; then
  echo "Verifying bundle checksums..."
  (cd "$BUNDLE_DIR" && grep -vE '  \./\._' manifest.sha256 | sha256sum -c -)
fi

if [[ "$YES" != "1" ]]; then
  cat <<EOF
Install Pathogen Workbench:
  app        ${APP_DIR}
  venv       ${VENV_DIR}
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

echo "Creating install directories..."
mkdir -p "$PREFIX" "$APP_DIR" "$STATE_DIR" "$TASK_ROOT" "$DB_ROOT"

echo "Installing application files..."
copy_tree "${BUNDLE_DIR}/app" "$APP_DIR"
copy_tree "${BUNDLE_DIR}/database" "$DB_ROOT"
copy_tree "${BUNDLE_DIR}/soft" "${APP_DIR}/soft"
copy_tree "${BUNDLE_DIR}/public" "${APP_DIR}/public"
if [[ -d "${BUNDLE_DIR}/demo_data" ]]; then
  copy_tree "${BUNDLE_DIR}/demo_data" "${APP_DIR}/demo_data"
fi
if [[ -d "${BUNDLE_DIR}/test_data" ]]; then
  copy_tree "${BUNDLE_DIR}/test_data" "${APP_DIR}/test_data"
fi

if [[ -f "${BUNDLE_DIR}/conda-runtime.tar.zst" ]]; then
  echo "Installing packed Conda runtime archive..."
  rm -rf "$PACKED_CONDA_ROOT"
  mkdir -p "$PACKED_CONDA_ROOT"
  extract_pack "${BUNDLE_DIR}/conda-runtime.tar.zst" "$PACKED_CONDA_ROOT"
elif [[ -d "${BUNDLE_DIR}/conda-runtime" ]]; then
  echo "Installing packed Conda runtime..."
  mkdir -p "$PACKED_CONDA_ROOT"
  copy_tree "${BUNDLE_DIR}/conda-runtime" "$PACKED_CONDA_ROOT"
fi
if [[ -d "${BUNDLE_DIR}/conda-packs" ]]; then
  echo "Installing packed Conda environments..."
  mkdir -p "$PACKED_CONDA_ENVS_DIR"
  shopt -s nullglob
  for pack in "${BUNDLE_DIR}/conda-packs"/*.tar.zst "${BUNDLE_DIR}/conda-packs"/*.tgz "${BUNDLE_DIR}/conda-packs"/*.tar.gz; do
    env_name="$(basename "$pack")"
    env_name="${env_name%.tar.zst}"
    env_name="${env_name%.tar.gz}"
    env_name="${env_name%.tgz}"
    env_dir="${PACKED_CONDA_ENVS_DIR}/${env_name}"
    rm -rf "$env_dir"
    echo "Extracting Conda env: ${env_name}"
    extract_pack "$pack" "$env_dir"
    if [[ -x "${env_dir}/bin/conda-unpack" ]]; then
      "${env_dir}/bin/conda-unpack"
    fi
  done
  shopt -u nullglob
fi

echo "Creating Python virtual environment..."
python3 -m venv "$VENV_DIR"
"${VENV_DIR}/bin/python" -m pip install --upgrade pip
"${VENV_DIR}/bin/python" -m pip install -r "${APP_DIR}/requirements-release.lock"

SECRET_KEY="$(random_token)"
INITIAL_ADMIN_PASSWORD="$(random_token)"
if [[ -f "$ENV_FILE" ]]; then
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
META_DATABASE_ROOT=${DB_ROOT}
CONDA_ROOT=${ACTIVE_CONDA_ROOT}
PORTAL_CONDA_ROOT=${ACTIVE_CONDA_ROOT}
CONDA_EXE=${ACTIVE_CONDA_ROOT}/bin/conda
PORTAL_REQUIRED_CONDA_ENVS=${DEFAULT_REQUIRED_ENVS}
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

echo "Seeding portal settings..."
set -a
# shellcheck disable=SC1090
source "$ENV_FILE"
set +a
cd "$APP_DIR"
"${VENV_DIR}/bin/python" - <<PY
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
    "conda_root": ${ACTIVE_CONDA_ROOT@Q},
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
WorkingDirectory=${APP_DIR}
EnvironmentFile=${ENV_FILE}
Environment=PYTHONUNBUFFERED=1
ExecStart=${VENV_DIR}/bin/gunicorn -w 2 -b 127.0.0.1:${PORT} 'bac_analysis_portal:create_app()'
Restart=always
RestartSec=5

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
rm -f /etc/nginx/sites-enabled/default

echo "Checking web deployment..."
"${VENV_DIR}/bin/python" "${APP_DIR}/scripts/check_deployment.py" \
  --project-root "$APP_DIR" \
  --env-file "$ENV_FILE" \
  --profile web

systemctl daemon-reload
systemctl enable --now "$SERVICE_NAME"
nginx -t
systemctl reload nginx

cat > "${APP_DIR}/INSTALL_ANALYSIS_ENVIRONMENT.md" <<'EOF'
# Full Analysis Environment

The portal is installed and can manage tasks, users, reports, and deployment
settings. Real sequencing analysis still requires Linux Conda environments with
the tools named in `PORTAL_REQUIRED_CONDA_ENVS`.

This portable bundle was built from macOS source assets, so it intentionally
does not include relocated Linux Conda environments. Do not copy macOS Conda
environments onto Ubuntu. Build or import Linux environments under:

```bash
/opt/miniconda3/envs/
```

Then run:

```bash
source /opt/pathogen-workbench/app/deployment/portal.env
/opt/pathogen-workbench/.venv_web/bin/python \
  /opt/pathogen-workbench/app/scripts/check_deployment.py \
  --project-root /opt/pathogen-workbench/app \
  --env-file /opt/pathogen-workbench/app/deployment/portal.env \
  --profile full-analysis
```
EOF

echo
echo "Pathogen Workbench installed."
echo "URL: http://$(hostname -I | awk '{print $1}')/login"
echo "Initial admin username: admin"
echo "Initial admin password: ${INITIAL_ADMIN_PASSWORD}"
echo "Password saved at: ${ADMIN_PASSWORD_FILE}"
echo
echo "Real analysis requires Linux Conda environments. See:"
echo "  ${APP_DIR}/INSTALL_ANALYSIS_ENVIRONMENT.md"
