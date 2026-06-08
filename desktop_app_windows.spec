# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path


APP_NAME = "PathogenWorkbench"
PROJECT_ROOT = Path(SPECPATH).resolve()
WINDOWS_ICON_PATH = PROJECT_ROOT / "bac_analysis_portal" / "static" / "app_icon.ico"


def collect_tree(relative_root: str) -> list[tuple[str, str]]:
    source_root = PROJECT_ROOT / relative_root
    if not source_root.exists():
        return []

    collected: list[tuple[str, str]] = []
    for path in source_root.rglob("*"):
        if not path.is_file():
            continue
        if path.name == ".DS_Store" or path.suffix == ".pyc" or "__pycache__" in path.parts:
            continue
        collected.append((str(path), str(path.relative_to(PROJECT_ROOT).parent)))
    return collected


datas = []
for folder in (
    "bac_analysis_portal/templates",
    "bac_analysis_portal/static",
    "public",
    "scripts",
):
    datas.extend(collect_tree(folder))

for filename in (
    "PathoSource.py",
    "Bac_assemble_260112_newformat.py",
    "CommunityAnalysis.py",
    "Virus_WTSpip.py",
):
    file_path = PROJECT_ROOT / filename
    if file_path.is_file():
        datas.append((str(file_path), "."))


icon_path = str(WINDOWS_ICON_PATH) if WINDOWS_ICON_PATH.is_file() else None

a = Analysis(
    ["run_bac_analysis_desktop.py"],
    pathex=[str(PROJECT_ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=["webview", "webview.platforms.winforms"],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name=APP_NAME,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    icon=icon_path,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name=APP_NAME,
)
