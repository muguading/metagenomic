#!/usr/bin/env python3
from __future__ import annotations

import argparse
import http.client
import json
import os
import platform
import subprocess
import shutil
import sys
import time
import zipfile
from pathlib import Path
from typing import Iterable
from urllib.parse import quote
from urllib.request import Request, urlopen


USER_AGENT = "metagenomic-download-test/1.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download a genome FASTA from NCBI using a GCF/GCA accession."
    )
    parser.add_argument("accession", help="NCBI assembly accession, e.g. GCF_000001405.40")
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory. Default: current working directory.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Maximum retry count for interrupted downloads. Default: 5.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed stage logs for debugging.",
    )
    parser.add_argument(
        "--mode",
        choices=("auto", "datasets", "api"),
        default="auto",
        help="Download backend: auto (prefer ncbi datasets), datasets, or api. Default: auto.",
    )
    return parser.parse_args()


def sanitize_name(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return "genome"
    normalized = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in text)
    while "__" in normalized:
        normalized = normalized.replace("__", "_")
    return normalized.strip("._") or "genome"


def normalize_accession(value: str) -> str:
    text = str(value or "").strip()
    if text.startswith(("GCF_", "GCA_")) and "." not in text:
        parts = text.rsplit("_", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return f"{parts[0]}.{parts[1]}"
    return text


def fasta_suffixes() -> tuple[str, ...]:
    return (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def find_fasta_members(names: Iterable[str]) -> list[str]:
    return [name for name in names if name.lower().endswith((".fa", ".fasta", ".fna"))]


def build_download_url(accession: str) -> str:
    return (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{quote(accession)}/download?include_annotation_type=GENOME_FASTA"
    )


def datasets_binary() -> str | None:
    script_dir = Path(__file__).resolve().parent
    system = platform.system().lower()
    candidates: list[Path] = []
    if system == "linux":
        candidates.append(script_dir / "datasets_linux")
    elif system == "darwin":
        candidates.append(script_dir / "datasets_macos")
    candidates.extend(
        [
            script_dir / "datasets",
            Path(shutil.which("datasets") or ""),
        ]
    )
    for candidate in candidates:
        if not candidate:
            continue
        if isinstance(candidate, Path) and str(candidate) == ".":
            continue
        path = Path(candidate)
        if not path.exists() or not path.is_file():
            continue
        try:
            mode = path.stat().st_mode
            if not os.access(path, os.X_OK):
                path.chmod(mode | 0o111)
        except OSError:
            pass
        if os.access(path, os.X_OK):
            return str(path)
    return None


def log(message: str, *, verbose: bool = True) -> None:
    if verbose:
        sys.stderr.write(message + "\n")
        sys.stderr.flush()


def download_with_datasets(accession: str, outdir: Path, *, verbose: bool, max_attempts: int = 3) -> Path:
    datasets = datasets_binary()
    if not datasets:
        raise RuntimeError("`datasets` command not found in PATH.")
    zip_path = outdir / f"{sanitize_name(accession)}.ncbi_dataset.zip"
    last_error: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        if zip_path.exists():
            zip_path.unlink(missing_ok=True)
        command = [
            datasets,
            "download",
            "genome",
            "accession",
            accession,
            "--filename",
            str(zip_path),
            "--include",
            "genome",
        ]
        log(f"[datasets] attempt {attempt}/{max_attempts}: {' '.join(command)}", verbose=verbose)
        process = subprocess.Popen(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        last_size = -1
        idle_rounds = 0
        while True:
            return_code = process.poll()
            current_size = zip_path.stat().st_size if zip_path.exists() else 0
            if current_size != last_size:
                idle_rounds = 0
                last_size = current_size
                if current_size > 0:
                    log(f"[datasets] zip growing: {current_size} bytes", verbose=verbose)
            else:
                idle_rounds += 1
                if verbose and idle_rounds % 5 == 0:
                    log("[datasets] waiting for download progress...", verbose=verbose)
            if return_code is not None:
                if return_code != 0:
                    last_error = RuntimeError(f"datasets download failed with exit code {return_code}")
                    log(f"[datasets] attempt {attempt}/{max_attempts} failed: exit code {return_code}", verbose=verbose)
                    break
                if not zip_path.exists():
                    last_error = RuntimeError(f"`datasets` finished but zip not found: {zip_path}")
                    log(f"[datasets] attempt {attempt}/{max_attempts} failed: zip missing", verbose=verbose)
                    break
                log(f"[datasets] Download complete: {zip_path}", verbose=verbose)
                return zip_path
            time.sleep(1)
        if attempt < max_attempts:
            time.sleep(min(2 * attempt, 6))
    if last_error is not None:
        raise last_error
    if not zip_path.exists():
        raise RuntimeError(f"`datasets` finished but zip not found: {zip_path}")
    return zip_path


def download_with_resume(url: str, destination: Path, max_attempts: int, *, verbose: bool) -> None:
    downloaded_bytes = destination.stat().st_size if destination.exists() else 0
    total_bytes = 0
    chunk_size = 1024 * 1024
    last_error: Exception | None = None

    for attempt in range(1, max_attempts + 1):
        log(f"[attempt {attempt}/{max_attempts}] Preparing request", verbose=verbose)
        request = Request(url, headers={"User-Agent": USER_AGENT})
        if downloaded_bytes > 0:
            request.add_header("Range", f"bytes={downloaded_bytes}-")
            log(f"[attempt {attempt}/{max_attempts}] Resuming from byte {downloaded_bytes}", verbose=verbose)
        try:
            log(f"[attempt {attempt}/{max_attempts}] Opening connection to NCBI", verbose=verbose)
            with urlopen(request, timeout=180) as response:
                log(f"[attempt {attempt}/{max_attempts}] Connected, HTTP status: {getattr(response, 'status', 'unknown')}", verbose=verbose)
                content_range = str(response.headers.get("Content-Range") or "").strip()
                content_length = str(response.headers.get("Content-Length") or "").strip()
                log(
                    f"[attempt {attempt}/{max_attempts}] Headers: Content-Length={content_length or '-'}, Content-Range={content_range or '-'}",
                    verbose=verbose,
                )
                if downloaded_bytes > 0 and not content_range:
                    downloaded_bytes = 0
                    if destination.exists():
                        destination.unlink(missing_ok=True)
                    log(f"[attempt {attempt}/{max_attempts}] Server does not support resume, restarting from byte 0", verbose=verbose)
                with destination.open("ab" if downloaded_bytes > 0 else "wb") as handle:
                    if content_range and "/" in content_range:
                        try:
                            total_bytes = int(content_range.rsplit("/", 1)[1])
                        except (TypeError, ValueError):
                            total_bytes = total_bytes
                    elif response.headers.get("Content-Length"):
                        try:
                            length = int(response.headers.get("Content-Length") or 0)
                        except (TypeError, ValueError):
                            length = 0
                        total_bytes = downloaded_bytes + length if downloaded_bytes > 0 else length

                    saw_first_chunk = False
                    while True:
                        try:
                            chunk = response.read(chunk_size)
                        except http.client.IncompleteRead as exc:
                            chunk = exc.partial or b""
                            if chunk:
                                handle.write(chunk)
                                downloaded_bytes += len(chunk)
                                report_progress(downloaded_bytes, total_bytes, attempt, retrying=True)
                                log(
                                    f"[attempt {attempt}/{max_attempts}] IncompleteRead after partial chunk, downloaded={downloaded_bytes}",
                                    verbose=verbose,
                                )
                            raise
                        if not chunk:
                            log(f"[attempt {attempt}/{max_attempts}] Response stream ended normally", verbose=verbose)
                            break
                        if not saw_first_chunk:
                            saw_first_chunk = True
                            log(
                                f"[attempt {attempt}/{max_attempts}] First data chunk received ({len(chunk)} bytes)",
                                verbose=verbose,
                            )
                        handle.write(chunk)
                        downloaded_bytes += len(chunk)
                        report_progress(downloaded_bytes, total_bytes, attempt, retrying=False)
                    log(
                        f"[attempt {attempt}/{max_attempts}] Download pass finished, file size now {downloaded_bytes} bytes",
                        verbose=verbose,
                    )
            sys.stderr.write("\n")
            log("[done] Download completed", verbose=verbose)
            return
        except http.client.IncompleteRead as exc:
            last_error = exc
            log(f"[attempt {attempt}/{max_attempts}] IncompleteRead caught: {exc}", verbose=verbose)
            time.sleep(min(2 * attempt, 8))
            continue
        except Exception as exc:
            last_error = exc
            log(f"[attempt {attempt}/{max_attempts}] Download error: {type(exc).__name__}: {exc}", verbose=verbose)
            if attempt >= max_attempts:
                break
            time.sleep(min(2 * attempt, 8))
            continue

    raise RuntimeError(f"Download failed after {max_attempts} attempts: {last_error}")


def report_progress(downloaded_bytes: int, total_bytes: int, attempt: int, retrying: bool) -> None:
    if total_bytes > 0:
        percent = downloaded_bytes * 100.0 / total_bytes
        message = f"\rDownloading: {percent:6.2f}% ({downloaded_bytes}/{total_bytes} bytes)"
    else:
        message = f"\rDownloading: {downloaded_bytes} bytes"
    if retrying:
        message += f" | connection interrupted, retry {attempt}"
    sys.stderr.write(message)
    sys.stderr.flush()


def raise_if_error_payload(downloaded_file: Path) -> None:
    head = downloaded_file.read_bytes()[:4096]
    text = head.decode("utf-8", errors="ignore").strip()
    if not text:
        raise RuntimeError("Downloaded file is empty.")
    if text.startswith(">"):
        return
    if text.startswith("{") or text.startswith("["):
        try:
            payload = json.loads(text)
        except json.JSONDecodeError:
            raise RuntimeError(f"Downloaded response is not a FASTA file: {text[:200]}")
        raise RuntimeError(f"NCBI returned an error payload: {json.dumps(payload, ensure_ascii=False)[:400]}")
    raise RuntimeError(f"Downloaded response is not a FASTA file: {text[:200]}")


def extract_fasta(downloaded_file: Path, accession: str, outdir: Path, *, verbose: bool) -> Path:
    log(f"[extract] Checking downloaded file: {downloaded_file}", verbose=verbose)
    if zipfile.is_zipfile(downloaded_file):
        log("[extract] Downloaded file is a zip archive, scanning members", verbose=verbose)
        with zipfile.ZipFile(downloaded_file) as archive:
            fasta_members = find_fasta_members(archive.namelist())
            if not fasta_members:
                raise RuntimeError("No FASTA file found inside the downloaded archive.")
            chosen = fasta_members[0]
            log(f"[extract] Selected FASTA member: {chosen}", verbose=verbose)
            suffix = Path(chosen).suffix or ".fa"
            target = outdir / f"{sanitize_name(accession)}{suffix}"
            with archive.open(chosen) as src, target.open("wb") as dest:
                shutil.copyfileobj(src, dest)
            log(f"[extract] FASTA extracted to: {target}", verbose=verbose)
            return target

    raise_if_error_payload(downloaded_file)
    suffix = downloaded_file.suffix or ".fa"
    target = outdir / f"{sanitize_name(accession)}{suffix}"
    if downloaded_file.resolve() != target.resolve():
        downloaded_file.replace(target)
    log(f"[extract] Non-zip response normalized to: {target}", verbose=verbose)
    return target


def main() -> int:
    args = parse_args()
    accession = normalize_accession(args.accession)
    if not accession.startswith(("GCF_", "GCA_")):
        sys.stderr.write("Accession must start with GCF_ or GCA_.\n")
        return 2

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    url = build_download_url(accession)

    sys.stderr.write(f"Downloading {accession} from NCBI...\n")
    sys.stderr.write(f"Mode: {args.mode}\n")
    if args.mode != "datasets":
        sys.stderr.write(f"URL: {url}\n")
    sys.stderr.write(f"Output directory: {outdir}\n")
    try:
        download_path: Path
        if args.mode == "datasets":
            download_path = download_with_datasets(accession, outdir, verbose=args.verbose, max_attempts=3)
        elif args.mode == "api":
            download_path = outdir / f"{sanitize_name(accession)}.download"
            download_with_resume(url, download_path, max_attempts=max(1, int(args.retries)), verbose=args.verbose)
        else:
            datasets = datasets_binary()
            if datasets:
                log(f"[auto] Found datasets binary: {datasets}", verbose=args.verbose)
                download_path = download_with_datasets(accession, outdir, verbose=args.verbose, max_attempts=3)
            else:
                log("[auto] `datasets` not found, falling back to API", verbose=args.verbose)
                download_path = outdir / f"{sanitize_name(accession)}.download"
                download_with_resume(url, download_path, max_attempts=max(1, int(args.retries)), verbose=args.verbose)
        fasta_path = extract_fasta(download_path, accession, outdir, verbose=args.verbose)
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        return 1

    sys.stderr.write(f"[done] Final FASTA path: {fasta_path}\n")
    sys.stdout.write(str(fasta_path) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
