from __future__ import annotations

import subprocess

import pytest

from metagenomic_refactor.mag_binning import _vamb_cuda_available


def test_vamb_cuda_available_returns_false_when_nvidia_smi_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    def fake_run(*args, **kwargs):
        raise FileNotFoundError("nvidia-smi not found")

    monkeypatch.setattr(subprocess, "run", fake_run)

    assert _vamb_cuda_available(tmp_path) is False
    assert "nvidia-smi not found" in (
        tmp_path / "vamb_cuda_check" / "host.stderr.log"
    ).read_text(encoding="utf-8")


def test_vamb_cuda_available_returns_false_when_env_probe_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    calls = []

    def fake_run(args, *posargs, **kwargs):
        calls.append(args)
        if args[:2] == ["nvidia-smi", "-L"]:
            return subprocess.CompletedProcess(args=args, returncode=0)
        raise FileNotFoundError("conda env not available")

    monkeypatch.setattr(subprocess, "run", fake_run)

    assert _vamb_cuda_available(tmp_path) is False
    assert len(calls) == 2
    assert "conda env not available" in (
        tmp_path / "vamb_cuda_check" / "env.stderr.log"
    ).read_text(encoding="utf-8")
