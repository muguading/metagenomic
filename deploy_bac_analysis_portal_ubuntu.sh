#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${ENV_FILE:-${PROJECT_DIR}/deployment/portal.env}"
ENV_EXAMPLE="${PROJECT_DIR}/deployment/portal.env.example"

if [[ ! -f "${ENV_FILE}" ]]; then
  if [[ -f "${ENV_EXAMPLE}" ]]; then
    mkdir -p "$(dirname "${ENV_FILE}")"
    cp "${ENV_EXAMPLE}" "${ENV_FILE}"
    cat >&2 <<WARNEOF
Created ${ENV_FILE} from portal.env.example.
Edit PORTAL_SECRET_KEY and PORTAL_INITIAL_ADMIN_PASSWORD before starting production.
WARNEOF
  else
    echo "ERROR: env file not found and template is missing: ${ENV_EXAMPLE}" >&2
    exit 2
  fi
fi

set -a
# shellcheck disable=SC1090
source "${ENV_FILE}"
set +a

VENV_DIR="${PROJECT_DIR}/.venv_web"
SERVICE_NAME="${SERVICE_NAME:-bac-analysis-portal}"
APP_USER="${APP_USER:-$(id -un)}"
APP_GROUP="${APP_GROUP:-$(id -gn)}"
APP_HOST="${APP_HOST:-127.0.0.1}"
APP_PORT="${APP_PORT:-5055}"
GUNICORN_WORKERS="${GUNICORN_WORKERS:-2}"
SERVER_NAME="${SERVER_NAME:-_}"
INSTALL_SYSTEM_PACKAGES="${INSTALL_SYSTEM_PACKAGES:-1}"
PORTAL_DB_PATH="${PORTAL_DB_PATH:-${PROJECT_DIR}/bac_analysis_portal.sqlite3}"
BAC_ANALYSIS_TASK_ROOT="${BAC_ANALYSIS_TASK_ROOT:-${PROJECT_DIR}/analysis_tasks}"

echo "Project dir: ${PROJECT_DIR}"
echo "Env file: ${ENV_FILE}"
echo "App user: ${APP_USER}:${APP_GROUP}"
echo "Bind: ${APP_HOST}:${APP_PORT}"
echo "Portal DB: ${PORTAL_DB_PATH}"
echo "Task root: ${BAC_ANALYSIS_TASK_ROOT}"

if [[ "${PORTAL_MODE:-}" == "production" ]]; then
  if [[ -z "${PORTAL_SECRET_KEY:-}" || "${PORTAL_SECRET_KEY}" == replace-with-* ]]; then
    echo "ERROR: PORTAL_SECRET_KEY must be set in ${ENV_FILE} for production." >&2
    exit 2
  fi
  if [[ -z "${PORTAL_INITIAL_ADMIN_PASSWORD:-}" || "${PORTAL_INITIAL_ADMIN_PASSWORD}" == replace-with-* || "${PORTAL_INITIAL_ADMIN_PASSWORD}" == "admin123" ]]; then
    echo "ERROR: PORTAL_INITIAL_ADMIN_PASSWORD must be a real strong password in ${ENV_FILE}." >&2
    exit 2
  fi
fi

mkdir -p "$(dirname "${PORTAL_DB_PATH}")" "${BAC_ANALYSIS_TASK_ROOT}"

if [[ "${INSTALL_SYSTEM_PACKAGES}" == "1" ]]; then
  echo "Installing Ubuntu system packages..."
  sudo apt update
  sudo apt install -y python3 python3-venv python3-pip nginx
fi

echo "Creating virtual environment..."
python3 -m venv "${VENV_DIR}"

echo "Installing Python dependencies..."
"${VENV_DIR}/bin/pip" install --upgrade pip
"${VENV_DIR}/bin/pip" install -r "${PROJECT_DIR}/requirements-release.lock"

echo "Writing Linux run script..."
cat > "${PROJECT_DIR}/run_analysis_portal_linux.sh" <<RUNEOF
#!/usr/bin/env bash
set -euo pipefail
cd "${PROJECT_DIR}"
set -a
source "${ENV_FILE}"
set +a
exec "${VENV_DIR}/bin/python" -m bac_analysis_portal.app
RUNEOF
chmod +x "${PROJECT_DIR}/run_analysis_portal_linux.sh"

echo "Writing systemd service template..."
cat > "${PROJECT_DIR}/${SERVICE_NAME}.service" <<SERVICEEOF
[Unit]
Description=Bac Analysis Portal
After=network.target

[Service]
User=${APP_USER}
Group=${APP_GROUP}
WorkingDirectory=${PROJECT_DIR}
ExecStart=${VENV_DIR}/bin/gunicorn -w ${GUNICORN_WORKERS} -b ${APP_HOST}:${APP_PORT} 'bac_analysis_portal:create_app()'
Restart=always
RestartSec=5
Environment=PYTHONUNBUFFERED=1
EnvironmentFile=${ENV_FILE}

[Install]
WantedBy=multi-user.target
SERVICEEOF

echo "Writing nginx site template..."
cat > "${PROJECT_DIR}/${SERVICE_NAME}.nginx.conf" <<NGINXEOF
server {
    listen 80;
    server_name ${SERVER_NAME};

    client_max_body_size 2G;

    location / {
        proxy_pass http://${APP_HOST}:${APP_PORT};
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
    }
}
NGINXEOF

cat <<INFOEOF

Deployment bootstrap complete.

Files created:
  ${PROJECT_DIR}/run_analysis_portal_linux.sh
  ${PROJECT_DIR}/${SERVICE_NAME}.service
  ${PROJECT_DIR}/${SERVICE_NAME}.nginx.conf

Deployment check:
  "${VENV_DIR}/bin/python" "${PROJECT_DIR}/scripts/check_deployment.py" --env-file "${ENV_FILE}"

Quick test:
  cd "${PROJECT_DIR}"
  "${VENV_DIR}/bin/python" -m bac_analysis_portal.app

Recommended production start:
  "${VENV_DIR}/bin/gunicorn" -w ${GUNICORN_WORKERS} -b ${APP_HOST}:${APP_PORT} 'bac_analysis_portal:create_app()'

To install the systemd service:
  sudo cp "${PROJECT_DIR}/${SERVICE_NAME}.service" /etc/systemd/system/${SERVICE_NAME}.service
  sudo systemctl daemon-reload
  sudo systemctl enable --now ${SERVICE_NAME}
  sudo systemctl status ${SERVICE_NAME}

To install the nginx config:
  sudo cp "${PROJECT_DIR}/${SERVICE_NAME}.nginx.conf" /etc/nginx/sites-available/${SERVICE_NAME}
  sudo ln -sf /etc/nginx/sites-available/${SERVICE_NAME} /etc/nginx/sites-enabled/${SERVICE_NAME}
  sudo nginx -t
  sudo systemctl reload nginx

If you need a public bind instead of nginx reverse proxy:
  APP_HOST=0.0.0.0 APP_PORT=5055 bash "${PROJECT_DIR}/deploy_bac_analysis_portal_ubuntu.sh"

INFOEOF
