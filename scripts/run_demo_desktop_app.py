from pathlib import Path
import os
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _ensure_runtime_python() -> None:
    venv_python = ROOT / ".venv_web" / "bin" / "python"
    if os.environ.get("METAGENOMIC_DEMO_VENV_READY") == "1":
        return
    try:
        import flask  # noqa: F401
        return
    except ModuleNotFoundError:
        if venv_python.exists():
            os.environ["METAGENOMIC_DEMO_VENV_READY"] = "1"
            os.execv(str(venv_python), [str(venv_python), __file__])
        raise


_ensure_runtime_python()


from bac_analysis_portal.demo_desktop_app import launch_demo_desktop_app


if __name__ == "__main__":
    launch_demo_desktop_app()
