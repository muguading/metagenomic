#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
APP_NAME="病原微生物分析工作台"
SPEC_FILE="${PROJECT_ROOT}/${APP_NAME}.spec"
DIST_APP="${PROJECT_ROOT}/dist/${APP_NAME}.app"
APPLICATIONS_APP="/Applications/${APP_NAME}.app"
PYTHON_BIN="${PYTHON_BIN:-}"

resolve_python_bin() {
  if [[ -n "${PYTHON_BIN}" ]]; then
    echo "${PYTHON_BIN}"
    return
  fi

  local candidates=(
    "${PROJECT_ROOT}/.venv_web/bin/python"
    "${PROJECT_ROOT}/.venv_web/bin/python3"
    "${PROJECT_ROOT}/.venv/bin/python"
    "${PROJECT_ROOT}/.venv/bin/python3"
  )

  local candidate
  for candidate in "${candidates[@]}"; do
    if [[ -x "${candidate}" ]]; then
      echo "${candidate}"
      return
    fi
  done

  command -v python3
}

PYTHON_BIN="$(resolve_python_bin)"

if [[ ! -f "${SPEC_FILE}" ]]; then
  echo "未找到 spec 文件: ${SPEC_FILE}" >&2
  exit 1
fi

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "未找到可用的 Python 解释器，请先安装 Python 3 或创建 .venv_web。" >&2
  exit 1
fi

if ! "${PYTHON_BIN}" -c "import PyInstaller" >/dev/null 2>&1; then
  echo "当前解释器缺少 PyInstaller: ${PYTHON_BIN}" >&2
  echo "请先执行: ${PYTHON_BIN} -m pip install -r ${PROJECT_ROOT}/requirements-web.txt" >&2
  exit 1
fi

if ! "${PYTHON_BIN}" -c "import flask, webview" >/dev/null 2>&1; then
  echo "当前解释器缺少桌面运行依赖: ${PYTHON_BIN}" >&2
  echo "请先执行: ${PYTHON_BIN} -m pip install -r ${PROJECT_ROOT}/requirements-web.txt" >&2
  exit 1
fi

echo "==> 项目目录: ${PROJECT_ROOT}"
echo "==> Python 解释器: ${PYTHON_BIN}"
echo "==> 使用 spec: ${SPEC_FILE}"
echo "==> 开始构建 macOS App..."
"${PYTHON_BIN}" -m PyInstaller "${SPEC_FILE}" --noconfirm

if [[ ! -d "${DIST_APP}" ]]; then
  echo "构建完成但未找到产物: ${DIST_APP}" >&2
  exit 1
fi

echo "==> 覆盖安装到应用程序目录..."
ditto "${DIST_APP}" "${APPLICATIONS_APP}"

echo "==> 构建完成"
echo "产物目录: ${DIST_APP}"
echo "安装位置: ${APPLICATIONS_APP}"

if [[ "${1:-}" == "--open" ]]; then
  echo "==> 正在打开 App..."
  open "${APPLICATIONS_APP}"
fi
