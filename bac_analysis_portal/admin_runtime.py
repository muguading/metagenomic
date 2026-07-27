from __future__ import annotations

import copy
import os
import re
import shutil
from datetime import datetime
from pathlib import Path

from metagenomic_refactor.common import resolve_conda_env_name

from .store import PortalStore
from .task_manager import ValidationError

WORKSTATION_MODULES = {
    "bacteria": {
        "key": "bacteria",
        "label": "细菌分析",
        "subtitle": "面向单菌细菌样本的组装、注释、耐药毒力与分型分析工作台。",
        "pipeline_script": "Bac_assemble_260112_newformat.py",
    },
    "virus": {
        "key": "virus",
        "label": "病毒分析",
        "subtitle": "面向病毒样本的参考驱动分析、鉴定与病毒分型工作台。",
        "pipeline_script": "Bac_assemble_260112_newformat.py",
    },
    "metagenome": {
        "key": "metagenome",
        "label": "宏基因组分析",
        "subtitle": "面向混合样本的宏基因组组装、物种识别与结果汇总工作台。",
        "pipeline_script": "Bac_assemble_260112_newformat.py",
    },
    "community": {
        "key": "community",
        "label": "群落分析",
        "subtitle": "面向多样本群落项目的 alpha / beta 多样性与 microeco Biomarker（LEfSe / RF）后分析工作台。",
        "pipeline_script": "CommunityAnalysis.py",
    },
    "pathosource": {
        "key": "pathosource",
        "label": "分子溯源",
        "subtitle": "围绕分子分型、谱系判读与溯源分析组织任务和结果视图。",
        "pipeline_script": "PathoSource.py",
    },
}

ADMIN_REALTIME_MONITOR_PRESETS = {
    "bacteria": {
        "standard": {
            "asm_type": "shortasm",
            "method": "spades",
        },
        "fast-track": {
            "asm_type": "shortasm",
            "method": "spades",
            "runflow": "fastp质控,耐药基因鉴定,毒力基因鉴定",
        },
    },
    "virus": {
        "reference": {
            "asm_type": "shortref",
            "method": "bwa",
        },
    },
    "metagenome": {
        "screening": {
            "asm_type": "shortasm",
            "method": "meta",
        },
        "deep": {
            "asm_type": "shortasm",
            "method": "meta",
            "thread": 16,
        },
    },
}

CONDA_ENV_SETTINGS = [
    {"key": "vfind", "label": "genomad_aux", "default": "genomad_aux", "group": "病毒发现 / cgMLST", "description": "VirSorter2、CheckV、chewBBACA 默认使用这个环境。"},
    {"key": "hamronization", "label": "amr_aux", "default": "amr_aux", "group": "耐药整合", "description": "hamronize、ResFinder、AMRFinder 等整合工具链。"},
    {"key": "rgi", "label": "amr_aux", "default": "amr_aux", "group": "耐药整合", "description": "RGI、pmga、salty、lissero 等分型补充工具。"},
    {"key": "rgi_new", "label": "amr_aux", "default": "amr_aux", "group": "耐药整合", "description": "宏基因组 bin 的新版 RGI 流程。"},
    {"key": "hostile", "label": "host_filter", "default": "host_filter", "group": "去宿主 / 清洗", "description": "hostile 去宿主流程。"},
    {"key": "kneaddata", "label": "host_filter", "default": "host_filter", "group": "去宿主 / 清洗", "description": "kneaddata 去宿主与去污染流程。"},
    {"key": "basalt", "label": "mag_aux", "default": "mag_aux", "group": "宏基因组", "description": "宏基因组 megahit / BASALT 分箱主流程。"},
    {"key": "coverm", "label": "mag_aux", "default": "mag_aux", "group": "宏基因组", "description": "宏基因组覆盖度与 TPM 定量。"},
    {"key": "cm210", "label": "cm210", "default": "cm210", "group": "宏基因组", "description": "checkm2、ectyper 等质量评估与补充分型。"},
    {"key": "gtdbtk", "label": "mag_aux", "default": "mag_aux", "group": "宏基因组", "description": "宏基因组 bin 的 GTDB-Tk 分类。"},
    {"key": "sistr_hicap", "label": "sistr_hicap", "default": "sistr_hicap", "group": "分型补充", "description": "SISTR 与 HICAP 补充分型工具。"},
    {"key": "qiime2", "label": "qiime2", "default": "qiime2", "group": "群落分析", "description": "Alpha / Beta 多样性、样本分类器与多样本生态统计。"},
    {"key": "microeco", "label": "microeco", "default": "microeco", "group": "群落分析", "description": "LEfSe、生态统计与可视化分析。"},
    {"key": "genovi", "label": "genovi", "default": "genovi", "group": "图谱与注释", "description": "genovi 图谱与 COG 统计。"},
    {"key": "plasflow", "label": "plasflow", "default": "plasflow", "group": "图谱与注释", "description": "质粒判定相关流程。"},
    {"key": "medaka", "label": "longread_aux", "default": "longread_aux", "group": "三代纠错 / 变异", "description": "三代数据抛光。"},
    {"key": "clair3", "label": "longread_aux", "default": "longread_aux", "group": "三代纠错 / 变异", "description": "三代数据变异检测。"},
    {"key": "chewie", "label": "chewie", "default": "chewie", "group": "分型补充", "description": "chewie / chewBBACA 相关分型分析。"},
    {"key": "tb_profiler", "label": "ncov", "default": "ncov", "group": "分型补充", "description": "结核分枝杆菌 tb-profiler 家系分析环境。"},
    {"key": "choleraefinder", "label": "choleraefinder", "default": "choleraefinder", "group": "分型补充", "description": "霍乱弧菌专用分型工具。"},
    {"key": "report_env", "label": "report_env", "default": "report_env", "group": "报告生成", "description": "R 报告输出与结果整合。"},
    {"key": "vlib", "label": "Vlib", "default": "Vlib", "group": "病毒旧流程兼容", "description": "旧版病毒分析兼容环境。"},
    {"key": "tngs", "label": "tNGS", "default": "tNGS", "group": "病毒旧流程兼容", "description": "旧版病毒流程补充环境。"},
]

CONDA_ENV_SETTINGS_BY_KEY = {item["key"]: item for item in CONDA_ENV_SETTINGS}


def _detect_conda_root() -> Path | None:
    candidates: list[Path] = []
    conda_exe = str(os.environ.get("CONDA_EXE", "") or "").strip()
    if conda_exe:
        candidates.append(Path(conda_exe).expanduser())
    which_conda = shutil.which("conda")
    if which_conda:
        candidates.append(Path(which_conda).expanduser())
    for candidate in candidates:
        try:
            resolved = candidate.resolve()
        except OSError:
            continue
        if resolved.name == "conda" and resolved.parent.name in {"bin", "Scripts"}:
            root = resolved.parent.parent
            if root.is_dir():
                return root
    return None

def _list_conda_env_names(conda_root: Path | None) -> list[str]:
    if conda_root is None or not conda_root.is_dir():
        return []
    env_names: set[str] = set()
    if (conda_root / "bin" / "python").is_file() or (conda_root / "python.exe").is_file():
        env_names.add("base")
    envs_dir = conda_root / "envs"
    if envs_dir.is_dir():
        for child in envs_dir.iterdir():
            if child.is_dir():
                env_names.add(child.name)
    return sorted(env_names, key=lambda item: item.lower())

def _resolve_workstation_module(raw_value: object) -> dict[str, str]:
    key = str(raw_value or "bacteria").strip().lower()
    if key in {"pathogen", "single"}:
        key = "bacteria"
    return WORKSTATION_MODULES.get(key, WORKSTATION_MODULES["bacteria"])

def _load_conda_env_settings(store: PortalStore) -> dict[str, str]:
    return {
        item["key"]: resolve_conda_env_name(
            str(
                store.get_setting(
                    f"conda_env_{item['key']}",
                    store.get_setting("conda_env_gtdbtk_caps", item["default"]) if item["key"] == "sistr_hicap" else item["default"],
                )
                or item["default"]
            ).strip() or item["default"]
        )
        for item in CONDA_ENV_SETTINGS
    }

def _load_conda_root_setting(store: PortalStore) -> str:
    raw = str(store.get_setting("conda_root", "") or "").strip()
    if raw:
        return raw
    detected = _detect_conda_root()
    return str(detected) if detected else ""

def _resolve_conda_exe_from_root(conda_root: str | Path | None) -> str:
    raw = str(conda_root or "").strip()
    if not raw:
        return "conda"
    root = Path(raw).expanduser()
    candidates = (
        root / "bin" / "conda",
        root / "condabin" / "conda",
        root / "Scripts" / "conda.exe",
        root / "conda.exe",
    )
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return str(candidates[0])

def _build_admin_realtime_monitor_payload(
    *,
    module_key: str,
    preset_key: str,
    preset_overrides: dict[str, object] | None,
    monitor_name: str,
    input_path: Path,
    output_root: Path,
    watch_stable_minutes: int,
    watch_poll_minutes: int,
    watch_max_samples: int,
) -> dict[str, object]:
    workstation = _resolve_workstation_module(module_key)
    resolved_key = workstation["key"]
    preset_group = ADMIN_REALTIME_MONITOR_PRESETS.get(resolved_key)
    if not preset_group:
        raise ValidationError(f"实时监控暂不支持模块：{resolved_key}")
    preset = copy.deepcopy(preset_group.get(preset_key))
    if not preset:
        raise ValidationError(f"{workstation['label']} 不支持预设参数：{preset_key}")
    if preset_overrides:
        if not isinstance(preset_overrides, dict):
            raise ValidationError(f"{workstation['label']} 的预设参数覆盖必须是对象")
        preset.update(preset_overrides)
    if resolved_key == "virus" and str(preset.get("species") or "").strip() in {"", "False", "false", "-", "none", "None"}:
        raise ValidationError("病毒分析请先选择目标病毒类型")
    timestamp_label = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_name = re.sub(r"[^0-9A-Za-z\u4e00-\u9fff._-]+", "_", monitor_name.strip() or input_path.name or resolved_key).strip("._-") or resolved_key
    task_name = f"{base_name}_{resolved_key}_{timestamp_label}"
    output_dir = output_root / resolved_key / task_name
    payload: dict[str, object] = {
        "task_name": task_name,
        "input_path": str(input_path),
        "output_dir": str(output_dir),
        "inputtype": "fastq",
        "watch_mode": "1",
        "watch_stable_minutes": watch_stable_minutes,
        "watch_poll_minutes": watch_poll_minutes,
        "watch_max_samples": watch_max_samples,
        "workstation_key": resolved_key,
        "thread": preset.get("thread", 10),
        "species": "False",
        "rna": "0",
        "fake_pip": 0,
    }
    payload.update(preset)
    return payload

def _first_allowed_module(user: dict[str, object] | None) -> str:
    allowed = user.get("allowed_modules") if isinstance(user, dict) else None
    if isinstance(allowed, list):
        for item in allowed:
            key = str(item or "").strip().lower()
            if key in {"pathogen", "single"}:
                key = "bacteria"
            if key in WORKSTATION_MODULES:
                return key
    return "bacteria"

def _resolve_home_and_desktop_paths() -> tuple[Path, Path]:
    home_path = Path.home().resolve()

    env_desktop = str(os.environ.get("XDG_DESKTOP_DIR", "")).strip()
    if env_desktop:
        desktop_candidate = Path(env_desktop.replace("$HOME", str(home_path))).expanduser()
        if desktop_candidate.is_dir():
            return home_path, desktop_candidate.resolve()

    user_dirs_file = home_path / ".config" / "user-dirs.dirs"
    if user_dirs_file.is_file():
        try:
            for raw_line in user_dirs_file.read_text(encoding="utf-8", errors="ignore").splitlines():
                line = raw_line.strip()
                if not line or line.startswith("#") or "XDG_DESKTOP_DIR" not in line:
                    continue
                _, value = line.split("=", 1)
                value = value.strip().strip('"').replace("$HOME", str(home_path))
                desktop_candidate = Path(value).expanduser()
                if desktop_candidate.is_dir():
                    return home_path, desktop_candidate.resolve()
        except OSError:
            pass

    for folder_name in ("Desktop", "桌面", "desktop"):
        desktop_candidate = home_path / folder_name
        if desktop_candidate.is_dir():
            return home_path, desktop_candidate.resolve()

    return home_path, (home_path / "Desktop").resolve()
