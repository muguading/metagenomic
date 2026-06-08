from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_ROOT = PROJECT_ROOT / "public" / "zika-tutorial"
LAT_LONGS_TEMPLATE = TEMPLATE_ROOT / "config" / "lat_longs.tsv"
DEFAULT_NEXTSTRAIN = Path("/Users/wuhhh/.nextstrain/cli-standalone/nextstrain")
COLOR_PALETTE = [
    "#3F51CC",
    "#4475CE",
    "#5AA793",
    "#7EB973",
    "#AABD55",
    "#CFB441",
    "#E39A39",
    "#E1622D",
    "#C7473B",
    "#7C62C9",
]
REGION_COLOR_MAP = {
    "Asia": "#447CCD",
    "Southeast Asia": "#4274CE",
    "Oceania": "#88BB6C",
    "Africa": "#CEB541",
    "Europe": "#E39B39",
    "South America": "#E56C2F",
    "North America": "#DC2F24",
    "Unknown": "#7893B3",
}
COUNTRY_REGION_MAP = {
    "china": "Asia",
    "thailand": "Southeast Asia",
    "singapore": "Southeast Asia",
    "vietnam": "Southeast Asia",
    "japan": "Asia",
    "usa": "North America",
    "united states": "North America",
    "american samoa": "Oceania",
    "brazil": "South America",
    "colombia": "South America",
    "ecuador": "South America",
    "venezuela": "South America",
    "guatemala": "North America",
    "honduras": "North America",
    "nicaragua": "North America",
    "panama": "North America",
    "dominican republic": "North America",
    "puerto rico": "North America",
    "french polynesia": "Oceania",
    "fiji": "Oceania",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a local Nextstrain workspace from sample-library selections.")
    parser.add_argument("--input", required=True, help="JSON manifest generated from selected sample-library rows.")
    parser.add_argument("--output", required=True, help="Workspace directory to create.")
    parser.add_argument("--ref", required=True, help="Reference GenBank file.")
    parser.add_argument("--thread", default="1")
    parser.add_argument("--nextstrain-bin", default=str(DEFAULT_NEXTSTRAIN))
    args, _unknown = parser.parse_known_args()
    return args


def load_manifest(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Manifest must be a JSON object")
    rows = payload.get("samples")
    if not isinstance(rows, list) or not rows:
        raise ValueError("Manifest does not contain any selected samples")
    return payload


def read_first_fasta_record(path: Path) -> tuple[str, str]:
    name = ""
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name and chunks:
                break
            name = line[1:].strip() or path.stem
            continue
        chunks.append(line)
    sequence = "".join(chunks).upper()
    if not sequence:
        raise ValueError(f"No FASTA sequence found in {path}")
    return name or path.stem, sequence


def normalize_country(value: object) -> str:
    return str(value or "").strip()


def infer_region(country: str) -> str:
    normalized = country.strip().lower()
    return COUNTRY_REGION_MAP.get(normalized, "Unknown")


def normalize_date(value: object, fallback: str) -> str:
    text = str(value or "").strip()
    if not text:
        return fallback
    return text.replace("XX", "01")


def load_lat_longs() -> dict[tuple[str, str], tuple[str, str]]:
    mapping: dict[tuple[str, str], tuple[str, str]] = {}
    for line in LAT_LONGS_TEMPLATE.read_text(encoding="utf-8").splitlines():
        parts = line.strip().split("\t")
        if len(parts) != 4:
            continue
        scope, name, latitude, longitude = parts
        mapping[(scope.strip().lower(), name.strip().lower())] = (latitude.strip(), longitude.strip())
    return mapping


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def apply_non_timetree_fallback(workspace_dir: Path) -> bool:
    snakefile_path = workspace_dir / "Snakefile"
    if not snakefile_path.is_file():
        return False
    original = snakefile_path.read_text(encoding="utf-8")
    updated = original
    for token in [
        "            --timetree \\\n",
        "            --coalescent {params.coalescent} \\\n",
        "            --date-confidence \\\n",
        "            --date-inference {params.date_inference} \\\n",
        "            --clock-filter-iqd {params.clock_filter_iqd}\n",
    ]:
        updated = updated.replace(token, "")
    if updated == original:
        return False
    snakefile_path.write_text(updated, encoding="utf-8")
    return True


def prepare_workspace(workspace_dir: Path, reference_gb: Path, manifest: dict) -> tuple[str, list[dict]]:
    if workspace_dir.exists():
        shutil.rmtree(workspace_dir)
    (workspace_dir / "data").mkdir(parents=True, exist_ok=True)
    (workspace_dir / "config").mkdir(parents=True, exist_ok=True)
    (workspace_dir / "envs").mkdir(parents=True, exist_ok=True)

    shutil.copy2(TEMPLATE_ROOT / "Snakefile", workspace_dir / "Snakefile")
    shutil.copy2(TEMPLATE_ROOT / "envs" / "nextstrain.yaml", workspace_dir / "envs" / "nextstrain.yaml")
    shutil.copy2(TEMPLATE_ROOT / "config" / "dropped_strains.txt", workspace_dir / "config" / "dropped_strains.txt")
    shutil.copy2(reference_gb, workspace_dir / "config" / "zika_outgroup.gb")

    build_name = str(manifest.get("build_name") or workspace_dir.name).strip() or workspace_dir.name
    today = datetime.now().strftime("%Y-%m-%d")
    config = {
        "title": build_name,
        "maintainers": [{"name": "本地样本数据库", "url": ""}],
        "build_url": "",
        "colorings": [
            {"key": "gt", "title": "Genotype", "type": "categorical"},
            {"key": "num_date", "title": "Date", "type": "continuous"},
            {"key": "author", "title": "Author", "type": "categorical"},
            {"key": "country", "title": "Country", "type": "categorical"},
            {"key": "region", "title": "Region", "type": "categorical"},
        ],
        "geo_resolutions": ["country", "region"],
        "panels": ["tree", "map", "entropy"],
        "display_defaults": {"map_triplicate": True},
        "filters": ["country", "region", "author"],
    }
    write_text(workspace_dir / "config" / "auspice_config.json", json.dumps(config, ensure_ascii=False, indent=2))

    rows = manifest["samples"]
    fallback_date = today
    fasta_lines: list[str] = []
    metadata_rows: list[dict[str, str]] = []
    countries_in_use: list[str] = []
    regions_in_use: list[str] = []

    for row in rows:
        fasta_path = Path(str(row.get("final_fasta_path") or "")).expanduser().resolve()
        if not fasta_path.is_file():
            raise FileNotFoundError(f"Missing FASTA for sample: {fasta_path}")
        sample_name = str(row.get("sample_name") or row.get("sample_key") or fasta_path.stem).strip()
        _header, sequence = read_first_fasta_record(fasta_path)
        fasta_lines.append(f">{sample_name}")
        fasta_lines.append(sequence)

        country = normalize_country(row.get("country")) or "Unknown"
        region = infer_region(country)
        location = row.get("location") if isinstance(row.get("location"), dict) else {}
        division = str(location.get("province") or location.get("city") or "").strip()
        city = str(location.get("detail") or location.get("district") or "").strip()
        collection_date = normalize_date(row.get("collection_date"), fallback_date)
        accession = str(row.get("genome_id") or row.get("sample_alias") or sample_name).strip()
        author = str(row.get("owner") or row.get("host_info") or "local").strip() or "local"
        metadata_rows.append(
            {
                "strain": sample_name,
                "virus": str(row.get("species_name") or manifest.get("species_name") or "Unknown virus").strip(),
                "accession": accession,
                "date": collection_date,
                "region": region,
                "country": country,
                "division": division,
                "city": city,
                "db": "sample_library",
                "segment": "",
                "authors": author,
                "url": "",
                "title": str(row.get("note") or build_name).strip(),
                "journal": "local build",
                "paper_url": "",
            }
        )
        if country not in countries_in_use:
            countries_in_use.append(country)
        if region not in regions_in_use:
            regions_in_use.append(region)

    write_text(workspace_dir / "data" / "sequences.fasta", "\n".join(fasta_lines) + "\n")
    metadata_columns = [
        "strain",
        "virus",
        "accession",
        "date",
        "region",
        "country",
        "division",
        "city",
        "db",
        "segment",
        "authors",
        "url",
        "title",
        "journal",
        "paper_url",
    ]
    metadata_lines = ["\t".join(metadata_columns)]
    for row in metadata_rows:
        metadata_lines.append("\t".join(str(row.get(column, "")).replace("\t", " ").strip() for column in metadata_columns))
    write_text(workspace_dir / "data" / "metadata.tsv", "\n".join(metadata_lines) + "\n")

    lat_longs_index = load_lat_longs()
    lat_long_lines: list[str] = []
    for country in countries_in_use:
        coords = lat_longs_index.get(("country", country.lower()))
        if coords:
            lat_long_lines.append(f"country\t{country.lower()}\t{coords[0]}\t{coords[1]}")
    for region in regions_in_use:
        coords = lat_longs_index.get(("region", region.lower()))
        if coords:
            lat_long_lines.append(f"region\t{region.lower()}\t{coords[0]}\t{coords[1]}")
    if not lat_long_lines:
        fallback = lat_longs_index.get(("country", "china"))
        if fallback:
            lat_long_lines.append(f"country\tchina\t{fallback[0]}\t{fallback[1]}")
    write_text(workspace_dir / "config" / "lat_longs.tsv", "\n".join(lat_long_lines) + "\n")

    color_lines: list[str] = []
    for index, country in enumerate(countries_in_use):
        color_lines.append(f"country\t{country.lower()}\t{COLOR_PALETTE[index % len(COLOR_PALETTE)]}")
    for region in regions_in_use:
        color_lines.append(f"region\t{region.lower()}\t{REGION_COLOR_MAP.get(region, REGION_COLOR_MAP['Unknown'])}")
    write_text(workspace_dir / "config" / "colors.tsv", "\n".join(color_lines) + "\n")

    manifest_snapshot = {
        "build_name": build_name,
        "species_name": str(manifest.get("species_name") or "").strip(),
        "sample_count": len(rows),
        "samples": metadata_rows,
        "reference_genbank": str(reference_gb),
    }
    write_text(workspace_dir / "build_manifest.json", json.dumps(manifest_snapshot, ensure_ascii=False, indent=2))
    return build_name, metadata_rows


def run_nextstrain(nextstrain_bin: Path, workspace_dir: Path, cpus: int) -> None:
    command = [
        str(nextstrain_bin),
        "build",
        "--cpus",
        str(max(1, cpus)),
        str(workspace_dir),
    ]
    result = subprocess.run(
        command,
        cwd=str(PROJECT_ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if result.returncode == 0:
        return

    failure_text = result.stdout or ""
    if "rate estimate is negative" in failure_text and apply_non_timetree_fallback(workspace_dir):
        print("[nextstrain] TreeTime 时间树失败，自动降级为普通系统发育树重新构建。")
        retry = subprocess.run(
            command,
            cwd=str(PROJECT_ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if retry.stdout:
            print(retry.stdout, end="" if retry.stdout.endswith("\n") else "\n")
        if retry.returncode == 0:
            return
        raise subprocess.CalledProcessError(retry.returncode, command)

    raise subprocess.CalledProcessError(result.returncode, command)


def main() -> int:
    args = parse_args()
    manifest_path = Path(args.input).expanduser().resolve()
    workspace_dir = Path(args.output).expanduser().resolve()
    reference_gb = Path(args.ref).expanduser().resolve()
    nextstrain_bin = Path(args.nextstrain_bin).expanduser().resolve()

    if not manifest_path.is_file():
        raise FileNotFoundError(f"Manifest file not found: {manifest_path}")
    if not reference_gb.is_file():
        raise FileNotFoundError(f"Reference GenBank not found: {reference_gb}")
    if not nextstrain_bin.is_file():
        raise FileNotFoundError(f"Nextstrain executable not found: {nextstrain_bin}")

    manifest = load_manifest(manifest_path)
    build_name, metadata_rows = prepare_workspace(workspace_dir, reference_gb, manifest)
    print(f"[nextstrain] Workspace prepared: {workspace_dir}")
    print(f"[nextstrain] Build name: {build_name}")
    print(f"[nextstrain] Samples: {len(metadata_rows)}")
    run_nextstrain(nextstrain_bin, workspace_dir, int(args.thread or 1))

    dataset_path = workspace_dir / "auspice" / "zika.json"
    if not dataset_path.is_file():
        raise FileNotFoundError(f"Build finished but dataset JSON is missing: {dataset_path}")
    print(f"[nextstrain] Dataset ready: {dataset_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
