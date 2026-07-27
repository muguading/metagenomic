from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from urllib.parse import unquote

from .parse_utils import _safe_int

_NCOV_NGDC_KB_DIR = Path(__file__).resolve().parent.parent / "database" / "virus" / "ncov" / "ngdc_mutation_kb"
_NCOV_NGDC_KB_INDEX_PATH = _NCOV_NGDC_KB_DIR / "knowledge_index.json"
_NCOV_ORF1A_AA_LENGTH = 4401


def _resolve_report_artifact_path(report_dir: Path, preferred_names: list[str], fallback_patterns: list[str] | None = None) -> Path:
    for name in preferred_names:
        candidate = report_dir / name
        if candidate.is_file():
            return candidate
    for pattern in fallback_patterns or []:
        matches = sorted(
            [item for item in report_dir.glob(pattern) if item.is_file() and not item.name.startswith(".")],
            key=lambda item: item.name.lower(),
        )
        if matches:
            return matches[0]
    return Path("")

def _parse_gff_attributes(raw_text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in str(raw_text or "").strip().split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        key = key.strip()
        value = unquote(value.strip())
        if key:
            attributes[key] = value
    return attributes

def _read_gff_genome_features(gff_path: Path) -> list[dict[str, object]]:
    if not gff_path.is_file():
        return []
    features: list[dict[str, object]] = []
    sequence_regions: list[tuple[str, int, int]] = []
    try:
        with gff_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                text = line.strip()
                if not text:
                    continue
                if text.startswith("##sequence-region"):
                    parts = text.split()
                    if len(parts) >= 4:
                        try:
                            sequence_regions.append((parts[1], int(parts[2]), int(parts[3])))
                        except (TypeError, ValueError):
                            pass
                    continue
                if text.startswith("#"):
                    continue
                parts = text.split("\t")
                if len(parts) < 9:
                    continue
                seqid, _source, feature_type, start_text, end_text, _score, strand, _phase, attributes_text = parts[:9]
                if feature_type not in {"gene", "CDS", "five_prime_UTR", "three_prime_UTR"}:
                    continue
                try:
                    start = int(start_text)
                    end = int(end_text)
                except (TypeError, ValueError):
                    continue
                if end < start:
                    continue
                attributes = _parse_gff_attributes(attributes_text)
                if feature_type in {"gene", "CDS"}:
                    primary_product = attributes.get("product") if feature_type == "CDS" else None
                    label = (
                        attributes.get("gene_name")
                        or attributes.get("gene")
                        or primary_product
                        or attributes.get("Name")
                        or attributes.get("locus_tag")
                        or attributes.get("ID")
                        or "Gene"
                    )
                    category = "gene"
                elif feature_type == "five_prime_UTR":
                    label = "5'UTR"
                    category = "utr"
                else:
                    label = "3'UTR"
                    category = "utr"
                features.append(
                    {
                        "seqid": seqid,
                        "feature_type": feature_type,
                        "category": category,
                        "label": label,
                        "start": start,
                        "end": end,
                        "length": end - start + 1,
                        "strand": strand or "+",
                    }
                )
    except OSError:
        return []
    if len(sequence_regions) > 1:
        offsets: dict[str, int] = {}
        cursor = 0
        for seqid, region_start, region_end in sequence_regions:
            offsets[seqid] = cursor - region_start + 1
            cursor += max(0, region_end - region_start + 1)
        for item in features:
            offset = offsets.get(str(item.get("seqid") or ""))
            if offset is None:
                continue
            item["segment_start"] = item.get("start")
            item["segment_end"] = item.get("end")
            item["start"] = int(item.get("start") or 0) + offset
            item["end"] = int(item.get("end") or 0) + offset
            item["length"] = int(item["end"]) - int(item["start"]) + 1
    gene_features = [item for item in features if str(item.get("feature_type") or "") == "gene"]
    utr_features = [item for item in features if str(item.get("category") or "") == "utr"]
    cds_features = [item for item in features if str(item.get("feature_type") or "") == "CDS"]
    display_features = utr_features + (gene_features if gene_features else cds_features)
    return sorted(display_features, key=lambda item: (int(item.get("start") or 0), int(item.get("end") or 0)))

def _simplify_virus_coverage_features(features: list[dict[str, object]], species_hint: str) -> list[dict[str, object]]:
    if not features:
        return []
    normalized_species = str(species_hint or "").strip().lower()
    if any(token in normalized_species for token in ("hantavirus", "orthohantavirus", "汉坦", "汉他")):
        return _simplify_orthohantavirus_coverage_features(features)
    if not any(token in normalized_species for token in ("ebola", "ebolavirus", "orthoebolavirus", "ebov", "埃博拉")):
        return features

    canonical_order = ["5'UTR", "NP", "VP35", "VP40", "GP", "VP30", "VP24", "L", "3'UTR"]
    canonical_rank = {label: index for index, label in enumerate(canonical_order)}

    def normalize_ebola_label(value: object, feature_type: object) -> str:
        text = str(value or "").strip()
        compact = re.sub(r"[^a-z0-9]+", "", text.lower())
        if str(feature_type or "").strip() == "five_prime_UTR":
            return "5'UTR"
        if str(feature_type or "").strip() == "three_prime_UTR":
            return "3'UTR"
        mapping = {
            "np": "NP",
            "nucleoprotein": "NP",
            "vp35": "VP35",
            "vp40": "VP40",
            "matrixprotein": "VP40",
            "gp": "GP",
            "sgp": "GP",
            "ssgp": "GP",
            "glycoprotein": "GP",
            "virionspikeglycoproteinprecursor": "GP",
            "vp30": "VP30",
            "polymerasecomplexprotein": "VP30",
            "vp24": "VP24",
            "l": "L",
            "polymerase": "L",
            "rnadependentrnapolymerase": "L",
        }
        return mapping.get(compact, text)

    merged: dict[str, dict[str, object]] = {}
    for feature in features:
        label = normalize_ebola_label(feature.get("label"), feature.get("feature_type"))
        if label not in canonical_rank:
            continue
        start = _safe_int(feature.get("start"))
        end = _safe_int(feature.get("end"))
        if start is None or end is None:
            continue
        if label not in merged:
            merged[label] = {**feature, "label": label, "start": start, "end": end}
            if label not in {"5'UTR", "3'UTR"}:
                merged[label]["feature_type"] = "gene"
                merged[label]["category"] = "gene"
            continue
        merged[label]["start"] = min(int(merged[label].get("start") or start), start)
        merged[label]["end"] = max(int(merged[label].get("end") or end), end)
        merged[label]["length"] = int(merged[label]["end"]) - int(merged[label]["start"]) + 1

    return sorted(
        merged.values(),
        key=lambda item: (
            canonical_rank.get(str(item.get("label") or ""), 999),
            int(item.get("start") or 0),
            int(item.get("end") or 0),
        ),
    )

def _simplify_orthohantavirus_coverage_features(features: list[dict[str, object]]) -> list[dict[str, object]]:
    canonical_order = ["5'UTR", "N", "GPC", "Gn", "Gc", "L", "3'UTR"]
    canonical_rank = {label: index for index, label in enumerate(canonical_order)}

    def normalize_hantavirus_label(value: object, feature_type: object) -> str:
        text = str(value or "").strip()
        compact = re.sub(r"[^a-z0-9]+", "", text.lower())
        if str(feature_type or "").strip() == "five_prime_UTR":
            return "5'UTR"
        if str(feature_type or "").strip() == "three_prime_UTR":
            return "3'UTR"
        if compact in {"n", "np", "nucleocapsid", "nucleocapsidprotein"}:
            return "N"
        if compact in {"gpc", "glycoprotein", "glycoproteinprecursor", "glycoproteinprecusor", "envelopeglycoproteinprecursor"}:
            return "GPC"
        if compact in {"gn", "glycoproteinn"}:
            return "Gn"
        if compact in {"gc", "glycoproteinc"}:
            return "Gc"
        if compact in {"l", "polymerase", "rdrp", "rnadependentrnapolymerase"}:
            return "L"
        return text

    merged: dict[str, dict[str, object]] = {}
    fallback: list[dict[str, object]] = []
    for feature in features:
        label = normalize_hantavirus_label(feature.get("label"), feature.get("feature_type"))
        start = _safe_int(feature.get("start"))
        end = _safe_int(feature.get("end"))
        if start is None or end is None:
            continue
        if label not in canonical_rank:
            fallback.append({**feature, "label": label, "start": start, "end": end})
            continue
        if label not in merged:
            merged[label] = {**feature, "label": label, "start": start, "end": end}
            if label not in {"5'UTR", "3'UTR"}:
                merged[label]["feature_type"] = "gene"
                merged[label]["category"] = "gene"
            continue
        merged[label]["start"] = min(int(merged[label].get("start") or start), start)
        merged[label]["end"] = max(int(merged[label].get("end") or end), end)
        merged[label]["length"] = int(merged[label]["end"]) - int(merged[label]["start"]) + 1

    display_features = list(merged.values()) or fallback
    return sorted(
        display_features,
        key=lambda item: (
            int(item.get("start") or 0),
            int(item.get("end") or 0),
            canonical_rank.get(str(item.get("label") or ""), 999),
        ),
    )

def _read_ncov_genome_features() -> list[dict[str, object]]:
    gff_path = Path(__file__).resolve().parent.parent / "database" / "virus" / "ncov" / "genomic.gff"
    return _read_gff_genome_features(gff_path)

def _normalize_sars_cov_2_gene_name(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    normalized = text.lower().replace(" ", "")
    mapping = {
        "orf1ab": "ORF1ab",
        "orf1a": "ORF1a",
        "orf1b": "ORF1b",
        "orf3a": "ORF3a",
        "orf6": "ORF6",
        "orf7a": "ORF7a",
        "orf7b": "ORF7b",
        "orf8": "ORF8",
        "orf9b": "ORF9b",
        "orf10": "ORF10",
        "spike": "S",
    }
    return mapping.get(normalized, text)

def _replace_first_numeric_token(text: str, replacement: int) -> str:
    return re.sub(r"\d+", str(replacement), text, count=1)

def _build_ncov_aa_alias_keys(gene: object, mutation: object) -> list[str]:
    gene_name = _normalize_sars_cov_2_gene_name(gene)
    change = str(mutation or "").strip()
    if not gene_name or not change:
        return []
    aliases = {f"{gene_name}:{change}"}
    match = re.search(r"(\d+)", change)
    if not match:
        return sorted(aliases)
    position = int(match.group(1))
    if gene_name == "ORF1ab":
        if position <= _NCOV_ORF1A_AA_LENGTH:
            aliases.add(f"ORF1a:{_replace_first_numeric_token(change, position)}")
        else:
            aliases.add(f"ORF1b:{_replace_first_numeric_token(change, position - _NCOV_ORF1A_AA_LENGTH)}")
    elif gene_name == "ORF1a":
        aliases.add(f"ORF1ab:{_replace_first_numeric_token(change, position)}")
    elif gene_name == "ORF1b":
        aliases.add(f"ORF1ab:{_replace_first_numeric_token(change, position + _NCOV_ORF1A_AA_LENGTH)}")
    return sorted(aliases)

def _split_comma_tokens(value: object) -> list[str]:
    if value is None:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]

def _load_ncov_ngdc_mutation_knowledge() -> dict:
    if not _NCOV_NGDC_KB_INDEX_PATH.is_file():
        return {}
    try:
        payload = json.loads(_NCOV_NGDC_KB_INDEX_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}

def _collect_ncov_ngdc_matches(columns: list[str], row: list[object]) -> dict:
    knowledge_payload = _load_ncov_ngdc_mutation_knowledge()
    by_aa_key = knowledge_payload.get("by_aa_key") if isinstance(knowledge_payload, dict) else {}
    if not isinstance(by_aa_key, dict) or not by_aa_key:
        return {
            "status": "missing",
            "source": "NGDC",
            "source_page": "https://ngdc.cncb.ac.cn/ncov/knowledge/mutation",
            "downloaded_at": "",
            "record_count": 0,
            "matched_record_count": 0,
            "aa_matches": {},
        }

    row_map = {
        str(column): row[index] if index < len(row) else ""
        for index, column in enumerate(columns)
    }
    aa_tokens: list[tuple[str, str, str]] = []
    for field_name, change_type in [
        ("aaSubstitutions", "substitution"),
        ("aaDeletions", "deletion"),
        ("aaInsertions", "insertion"),
    ]:
        for token in _split_comma_tokens(row_map.get(field_name)):
            gene, _, change = token.partition(":")
            if not gene or not change:
                continue
            aa_tokens.append((gene.strip(), change.strip(), change_type))

    aa_matches: dict[str, list[dict]] = {}
    matched_total = 0
    for gene, change, change_type in aa_tokens:
        label = f"{_normalize_sars_cov_2_gene_name(gene)}:{change}"
        record_pool: list[dict] = []
        seen: set[tuple[str, str, str, str, str]] = set()
        for alias in _build_ncov_aa_alias_keys(gene, change):
            for item in by_aa_key.get(alias, []):
                if not isinstance(item, dict):
                    continue
                dedupe_key = (
                    str(item.get("section_key") or ""),
                    str(item.get("gene") or ""),
                    str(item.get("mutation") or ""),
                    str(item.get("effect") or ""),
                    str(item.get("detail") or ""),
                )
                if dedupe_key in seen:
                    continue
                seen.add(dedupe_key)
                record_pool.append({
                    "source": str(item.get("source") or "NGDC"),
                    "section_key": str(item.get("section_key") or "").strip(),
                    "section_label": str(item.get("section_label") or "").strip(),
                    "gene": str(item.get("gene") or "").strip(),
                    "mutation": str(item.get("mutation") or "").strip(),
                    "display_label": f"{str(item.get('gene') or '').strip()}:{str(item.get('mutation') or '').strip()}".strip(":"),
                    "genomic_position": str(item.get("genomic_position") or "").strip(),
                    "nuc_change": str(item.get("nuc_change") or "").strip(),
                    "effect": str(item.get("effect") or "").strip(),
                    "effect_zh": str(item.get("effect_zh") or item.get("effect") or "").strip(),
                    "detail": str(item.get("detail") or "").strip(),
                    "detail_zh": str(item.get("detail_zh") or item.get("detail") or "").strip(),
                    "method": str(item.get("method") or "").strip(),
                    "method_zh": str(item.get("method_zh") or item.get("method") or "").strip(),
                    "pmid": str(item.get("pmid") or "").strip(),
                    "extra": item.get("extra") if isinstance(item.get("extra"), dict) else {},
                })
        if record_pool:
            aa_matches[label] = record_pool
            matched_total += len(record_pool)

    return {
        "status": "ready",
        "source": str(knowledge_payload.get("source") or "NGDC"),
        "source_page": str(knowledge_payload.get("source_page") or "https://ngdc.cncb.ac.cn/ncov/knowledge/mutation"),
        "downloaded_at": str(knowledge_payload.get("downloaded_at") or ""),
        "record_count": int(knowledge_payload.get("record_count") or 0),
        "matched_record_count": matched_total,
        "aa_matches": aa_matches,
    }
