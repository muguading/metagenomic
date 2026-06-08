#!/usr/bin/env bash

set -euo pipefail

ENV_NAME="${1:-mag_binning}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_YAML="${PROJECT_ROOT}/envs/mag_binning_env.yml"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda 未安装或不在 PATH 中。" >&2
    exit 1
fi

if [ ! -f "${ENV_YAML}" ]; then
    echo "ERROR: 环境文件不存在: ${ENV_YAML}" >&2
    exit 1
fi

if command -v mamba >/dev/null 2>&1; then
    CONDA_BIN="mamba"
else
    CONDA_BIN="conda"
fi

TMP_YAML="$(mktemp "${TMPDIR:-/tmp}/mag_binning_env.XXXXXX.yml")"
trap 'rm -f "${TMP_YAML}"' EXIT

sed "s/^name: .*/name: ${ENV_NAME}/" "${ENV_YAML}" > "${TMP_YAML}"

echo "[1/3] 使用 ${CONDA_BIN} 创建环境: ${ENV_NAME}"
"${CONDA_BIN}" env create -f "${TMP_YAML}" --force

echo "[2/3] 校验关键软件"
conda run -n "${ENV_NAME}" python --version
conda run -n "${ENV_NAME}" megahit --version || true
conda run -n "${ENV_NAME}" bwa 2>&1 | head -n 2 || true
conda run -n "${ENV_NAME}" samtools --version | head -n 2 || true
conda run -n "${ENV_NAME}" metabat2 --help >/dev/null
conda run -n "${ENV_NAME}" SemiBin2 --help >/dev/null
conda run -n "${ENV_NAME}" vamb --help >/dev/null
conda run -n "${ENV_NAME}" DAS_Tool --help >/dev/null

echo "[3/3] 完成"
echo "激活命令:"
echo "  conda activate ${ENV_NAME}"
echo
echo "建议测试:"
echo "  python -m metagenomic_refactor.mag_binning --help"
