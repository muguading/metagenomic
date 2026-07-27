from __future__ import annotations

import importlib.util
import re
import sys
import tempfile
from functools import lru_cache
from pathlib import Path

from .table_io import _read_tsv_rows

def _read_pathosource_tree(path: Path, label: str) -> dict:
    if not path.is_file():
        return {"status": "empty", "label": label, "newick": "", "leaf_count": 0, "file_name": path.name}
    try:
        newick = path.read_text(encoding="utf-8", errors="ignore").strip()
    except OSError:
        return {"status": "empty", "label": label, "newick": "", "leaf_count": 0, "file_name": path.name}
    if not newick:
        return {"status": "empty", "label": label, "newick": "", "leaf_count": 0, "file_name": path.name}
    leaf_count = len(re.findall(r"(?:(?<=\()|(?<=,))\s*([^():;,]+)\s*:", newick))
    return {
        "status": "ready",
        "label": label,
        "newick": newick,
        "leaf_count": leaf_count,
        "char_count": len(newick),
        "file_name": path.name,
    }

def _build_norovirus_gene_phylogeny(report_dir: Path) -> dict:
    phylogeny_root = report_dir / f"{report_dir.name}_norovirus_reference_selection" / "phylogeny"
    summary_path = phylogeny_root / "summary.tsv"
    summary_table = _read_tsv_rows(summary_path)
    summary_columns = summary_table.get("columns", [])
    summary_rows = summary_table.get("rows", [])
    if not summary_columns or not summary_rows:
        return {"status": "empty", "trees": []}
    trees: list[dict[str, object]] = []
    for row in summary_rows:
        if not isinstance(row, list):
            continue
        record = {
            summary_columns[index]: row[index] if index < len(row) else ""
            for index in range(len(summary_columns))
        }
        status = str(record.get("status") or "").strip().lower()
        tree_path_text = str(record.get("tree_path") or "").strip()
        if status != "ready" or not tree_path_text:
            continue
        tree_path = Path(tree_path_text)
        gene_label = str(record.get("gene") or "").strip().upper()
        genogroup = str(record.get("genogroup") or "").strip().upper()
        subtype = str(record.get("subtype") or "").strip()
        label = f"Norovirus {gene_label} {genogroup or '-'} 片段系统发育树"
        tree_section = _read_pathosource_tree(tree_path, label)
        tree_section["gene"] = gene_label
        tree_section["genogroup"] = genogroup
        tree_section["subtype"] = subtype
        tree_section["member_count"] = str(record.get("member_count") or "")
        trees.append(_attach_itol_payload(tree_section, tree_path))
    return {
        "status": "ready" if trees else "empty",
        "trees": trees,
        "summary": summary_table,
    }

def _build_rhinovirus_gene_phylogeny(report_dir: Path) -> dict:
    phylogeny_root = report_dir / f"{report_dir.name}_rhinovirus_reference_selection" / "phylogeny"
    summary_path = phylogeny_root / "summary.tsv"
    summary_table = _read_tsv_rows(summary_path)
    summary_columns = summary_table.get("columns", [])
    summary_rows = summary_table.get("rows", [])
    if not summary_columns or not summary_rows:
        return {"status": "empty", "trees": []}
    trees: list[dict[str, object]] = []
    for row in summary_rows:
        if not isinstance(row, list):
            continue
        record = {
            summary_columns[index]: row[index] if index < len(row) else ""
            for index in range(len(summary_columns))
        }
        status = str(record.get("status") or "").strip().lower()
        tree_path_text = str(record.get("tree_path") or "").strip()
        if status != "ready" or not tree_path_text:
            continue
        tree_path = Path(tree_path_text)
        gene_label = str(record.get("gene") or "").strip().upper()
        genogroup = str(record.get("genogroup") or "").strip().upper()
        subtype = str(record.get("subtype") or "").strip()
        label = f"Rhinovirus {gene_label} {genogroup or '-'} 片段系统发育树"
        tree_section = _read_pathosource_tree(tree_path, label)
        tree_section["gene"] = gene_label
        tree_section["genogroup"] = genogroup
        tree_section["subtype"] = subtype
        tree_section["member_count"] = str(record.get("member_count") or "")
        trees.append(_attach_itol_payload(tree_section, tree_path))
    return {
        "status": "ready" if trees else "empty",
        "trees": trees,
        "summary": summary_table,
    }

def _build_enterovirus_gene_phylogeny(report_dir: Path) -> dict:
    phylogeny_root = report_dir / f"{report_dir.name}_enterovirus_reference_selection" / "phylogeny"
    summary_path = phylogeny_root / "summary.tsv"
    summary_table = _read_tsv_rows(summary_path)
    summary_columns = summary_table.get("columns", [])
    summary_rows = summary_table.get("rows", [])
    if not summary_columns or not summary_rows:
        return {"status": "empty", "trees": []}
    trees: list[dict[str, object]] = []
    for row in summary_rows:
        if not isinstance(row, list):
            continue
        record = {
            summary_columns[index]: row[index] if index < len(row) else ""
            for index in range(len(summary_columns))
        }
        status = str(record.get("status") or "").strip().lower()
        tree_path_text = str(record.get("tree_path") or "").strip()
        if status != "ready" or not tree_path_text:
            continue
        tree_path = Path(tree_path_text)
        gene_label = str(record.get("gene") or "").strip().upper()
        genogroup = str(record.get("genogroup") or "").strip().upper()
        subtype = str(record.get("subtype") or "").strip()
        label = f"Enterovirus {gene_label or 'VP1'} {genogroup or '-'} 片段系统发育树"
        tree_section = _read_pathosource_tree(tree_path, label)
        tree_section["gene"] = gene_label
        tree_section["genogroup"] = genogroup
        tree_section["subtype"] = subtype
        tree_section["member_count"] = str(record.get("member_count") or "")
        trees.append(_attach_itol_payload(tree_section, tree_path))
    return {
        "status": "ready" if trees else "empty",
        "trees": trees,
        "summary": summary_table,
    }

def _build_astroviridae_gene_phylogeny(report_dir: Path) -> dict:
    phylogeny_root = report_dir / f"{report_dir.name}_astroviridae_reference_selection" / "phylogeny"
    summary_path = phylogeny_root / "summary.tsv"
    summary_table = _read_tsv_rows(summary_path)
    summary_columns = summary_table.get("columns", [])
    summary_rows = summary_table.get("rows", [])
    if not summary_columns or not summary_rows:
        return {"status": "empty", "trees": []}
    trees: list[dict[str, object]] = []
    for row in summary_rows:
        if not isinstance(row, list):
            continue
        record = {
            summary_columns[index]: row[index] if index < len(row) else ""
            for index in range(len(summary_columns))
        }
        status = str(record.get("status") or "").strip().lower()
        tree_path_text = str(record.get("tree_path") or "").strip()
        if status != "ready" or not tree_path_text:
            continue
        tree_path = Path(tree_path_text)
        gene_label = str(record.get("gene") or "").strip().upper()
        genus = str(record.get("genogroup") or "").strip()
        subtype = str(record.get("subtype") or "").strip()
        label = f"Astrovirus {gene_label or 'ORF2'} {genus or '-'} 片段系统发育树"
        tree_section = _read_pathosource_tree(tree_path, label)
        tree_section["gene"] = gene_label
        tree_section["genogroup"] = genus
        tree_section["subtype"] = subtype
        tree_section["member_count"] = str(record.get("member_count") or "")
        trees.append(_attach_itol_payload(tree_section, tree_path))
    return {
        "status": "ready" if trees else "empty",
        "trees": trees,
        "summary": summary_table,
    }

def _build_seasonal_hcov_spike_phylogeny(report_dir: Path) -> dict:
    phylogeny_root = report_dir / f"{report_dir.name}_seasonal_hcov_reference_selection" / "phylogeny"
    summary_path = phylogeny_root / "summary.tsv"
    summary_table = _read_tsv_rows(summary_path)
    summary_columns = summary_table.get("columns", [])
    summary_rows = summary_table.get("rows", [])
    if not summary_columns or not summary_rows:
        return {"status": "empty", "trees": []}
    trees: list[dict[str, object]] = []
    for row in summary_rows:
        if not isinstance(row, list):
            continue
        record = {
            summary_columns[index]: row[index] if index < len(row) else ""
            for index in range(len(summary_columns))
        }
        status = str(record.get("status") or "").strip().lower()
        tree_path_text = str(record.get("tree_path") or "").strip()
        if status != "ready" or not tree_path_text:
            continue
        tree_path = Path(tree_path_text)
        hcov_type = str(record.get("hcov_type") or "").strip()
        subtype = str(record.get("subtype") or "").strip()
        label = f"{hcov_type or 'Seasonal HCoV'} Spike 系统发育树"
        tree_section = _read_pathosource_tree(tree_path, label)
        tree_section["gene"] = "S"
        tree_section["subtype"] = subtype
        tree_section["genogroup"] = hcov_type
        tree_section["member_count"] = str(record.get("member_count") or "")
        trees.append(_attach_itol_payload(tree_section, tree_path))
    return {
        "status": "ready" if trees else "empty",
        "trees": trees,
        "summary": summary_table,
    }

def _load_itol_helper(helper_name: str, helper_path: str):
    spec = importlib.util.spec_from_file_location(helper_name, helper_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load helper: {helper_name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[helper_name] = module
    spec.loader.exec_module(module)
    return module

def _attach_itol_payload(tree_section: dict, tree_path: Path, tree_text: str | None = None) -> dict:
    if tree_section.get("status") != "ready" or not tree_path.is_file():
        return tree_section
    project_root = Path(__file__).resolve().parents[1]
    nwk2json_path = project_root / "public" / "itol" / "nwk2json.py"
    treeinfo_path = project_root / "public" / "itol" / "newick2treeinfo.py"
    convert_path = tree_path
    tmp_path: Path | None = None
    try:
        if tree_text is not None:
            with tempfile.NamedTemporaryFile("w", suffix=".nwk", encoding="utf-8", delete=False) as handle:
                handle.write(tree_text)
                tmp_path = Path(handle.name)
            convert_path = tmp_path
        nwk2json_module = _load_itol_helper("portal_itol_nwk2json", str(nwk2json_path))
        treeinfo_module = _load_itol_helper("portal_itol_treeinfo", str(treeinfo_path))
        tree_json = nwk2json_module.convert_newick_to_itol_json(convert_path, tree_name=tree_path.stem)
        treeinfo_rows = treeinfo_module.newick_to_treeinfo(convert_path)
    except Exception as exc:
        tree_section["itol_error"] = str(exc)
        return tree_section
    finally:
        if tmp_path is not None:
            try:
                tmp_path.unlink()
            except OSError:
                pass
    tree_section["itol_tree_json"] = tree_json
    tree_section["itol_treeinfo"] = {
        "columns": ["parent", "node", "branch.length", "label", "isTip", "x", "y", "branch", "angle"],
        "rows": [
            [
                row.get("parent", ""),
                row.get("node", ""),
                row.get("branch.length", ""),
                row.get("label", ""),
                row.get("isTip", ""),
                row.get("x", ""),
                row.get("y", ""),
                row.get("branch", ""),
                row.get("angle", ""),
            ]
            for row in treeinfo_rows
        ],
    }
    return tree_section

def _discover_nextclade_tree_path(report_dir: Path) -> Path | None:
    def _has_content(path: Path) -> bool:
        try:
            return path.is_file() and path.stat().st_size > 0
        except OSError:
            return False

    preferred = report_dir / "nextclade_output" / "nextclade.nwk"
    if _has_content(preferred):
        return preferred
    candidates = sorted(
        [
            path
            for path in report_dir.rglob("*.nwk")
            if not path.name.startswith(".") and _has_content(path)
        ],
        key=lambda item: (
            0 if "nextclade" in str(item.relative_to(report_dir)).lower() else 1,
            len(item.parts),
            item.name.lower(),
        ),
    )
    return candidates[0] if candidates else None

def _load_newick_parser_module():
    project_root = Path(__file__).resolve().parents[1]
    parser_path = project_root / "public" / "itol" / "newick_parser.py"
    return _load_itol_helper("portal_itol_newick_parser", str(parser_path))

def _annotate_newick_tree(root) -> tuple[list[object], dict[int, object], dict[int, float]]:
    parent_map: dict[int, object] = {}
    distance_map: dict[int, float] = {}
    leaves: list[object] = []
    stack: list[tuple[object, object | None, float]] = [(root, None, 0.0)]
    while stack:
        node, parent, distance = stack.pop()
        parent_map[id(node)] = parent
        distance_map[id(node)] = distance
        children = list(getattr(node, "children", []) or [])
        if not children:
            leaves.append(node)
            continue
        for child in reversed(children):
            child_length = getattr(child, "branch_length", None)
            stack.append((child, node, distance + (float(child_length) if child_length is not None else 0.0)))
    return leaves, parent_map, distance_map

def _compute_leaf_distance(node_a: object, node_b: object, parent_map: dict[int, object], distance_map: dict[int, float]) -> float:
    ancestors: dict[int, float] = {}
    cursor = node_a
    while cursor is not None:
        ancestors[id(cursor)] = distance_map.get(id(cursor), 0.0)
        cursor = parent_map.get(id(cursor))
    cursor = node_b
    while cursor is not None:
        cursor_id = id(cursor)
        if cursor_id in ancestors:
            lca_distance = ancestors[cursor_id]
            return distance_map.get(id(node_a), 0.0) + distance_map.get(id(node_b), 0.0) - (2.0 * lca_distance)
        cursor = parent_map.get(cursor_id)
    return float("inf")

def _clone_newick_node(parser_module, source: object, *, children: list[object] | None = None, branch_length: float | None = None):
    clone = parser_module.NewickNode(
        name=getattr(source, "name", None),
        branch_length=getattr(source, "branch_length", None) if branch_length is None else branch_length,
        confidence=getattr(source, "confidence", None),
        children=list(children or []),
    )
    return clone

def _prune_newick_tree_to_labels(root: object, keep_labels: set[str], parser_module):
    children = list(getattr(root, "children", []) or [])
    if not children:
        label = str(getattr(root, "name", "") or "")
        if label in keep_labels:
            return _clone_newick_node(parser_module, root, children=[])
        return None

    pruned_children: list[object] = []
    for child in children:
        pruned = _prune_newick_tree_to_labels(child, keep_labels, parser_module)
        if pruned is not None:
            pruned_children.append(pruned)
    if not pruned_children:
        return None

    if len(pruned_children) == 1:
        only_child = pruned_children[0]
        merged_length = (getattr(only_child, "branch_length", None) or 0.0) + (getattr(root, "branch_length", None) or 0.0)
        only_child.branch_length = merged_length if merged_length else None
        return only_child

    return _clone_newick_node(parser_module, root, children=pruned_children)

def _serialize_newick_node(node: object, *, is_root: bool = False) -> str:
    children = list(getattr(node, "children", []) or [])
    label = str(getattr(node, "name", "") or "")
    confidence = getattr(node, "confidence", None)
    parts: list[str] = []
    if children:
        parts.append("(" + ",".join(_serialize_newick_node(child, is_root=False) for child in children) + ")")
        if label:
            parts.append(label)
        elif confidence is not None:
            if float(confidence).is_integer():
                parts.append(str(int(confidence)))
            else:
                parts.append(f"{float(confidence):g}")
    elif label:
        parts.append(label)
    branch_length = getattr(node, "branch_length", None)
    if branch_length is not None and not is_root:
        parts.append(f":{float(branch_length):g}")
    return "".join(parts)

def _build_nextclade_nearest_subtree(newick_text: str, query_labels: list[str], neighbor_limit: int = 50) -> tuple[str, list[str]]:
    parser_module = _load_newick_parser_module()
    root = parser_module.parse_newick_text(newick_text)
    leaves, parent_map, distance_map = _annotate_newick_tree(root)
    label_to_leaf: dict[str, object] = {}
    for leaf in leaves:
        label = str(getattr(leaf, "name", "") or "")
        if label and label not in label_to_leaf:
            label_to_leaf[label] = leaf
    target_labels = [label for label in query_labels if label in label_to_leaf]
    if not target_labels:
        return newick_text, []

    target_leaves = [label_to_leaf[label] for label in target_labels]
    ranked: list[tuple[float, str]] = []
    for label, leaf in label_to_leaf.items():
        if label in target_labels:
            continue
        distance = min(_compute_leaf_distance(target_leaf, leaf, parent_map, distance_map) for target_leaf in target_leaves)
        ranked.append((distance, label))
    ranked.sort(key=lambda item: (item[0], item[1]))
    selected_labels = set(target_labels)
    selected_labels.update(label for _distance, label in ranked[: max(0, neighbor_limit)])
    pruned_root = _prune_newick_tree_to_labels(root, selected_labels, parser_module)
    if pruned_root is None:
        return newick_text, target_labels
    subtree_text = _serialize_newick_node(pruned_root, is_root=True) + ";"
    ordered_neighbors = [label for _distance, label in ranked[: max(0, neighbor_limit)]]
    return subtree_text, ordered_neighbors

def _build_nextclade_phylogeny_tree(report_dir: Path, sequence_names: object = "") -> dict:
    tree_path = _discover_nextclade_tree_path(report_dir)
    if tree_path is None:
        return {"status": "empty", "label": "Nextclade 系统发育树", "newick": "", "leaf_count": 0, "file_name": "nextclade.nwk"}
    tree_section = _read_pathosource_tree(tree_path, "Nextclade 系统发育树")
    newick_text = str(tree_section.get("newick") or "")
    if isinstance(sequence_names, str):
        raw_names = [sequence_names]
    elif isinstance(sequence_names, list):
        raw_names = [str(item or "") for item in sequence_names]
    else:
        raw_names = []
    sanitized_names: list[dict[str, str]] = []
    used_names: set[str] = set()
    for display_name in sorted({item.strip() for item in raw_names if item.strip()}, key=len, reverse=True):
        if display_name not in newick_text:
            continue
        base_name = re.split(r"\s+|\|", display_name, maxsplit=1)[0].strip() or "query_sequence"
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", base_name).strip("_") or "query_sequence"
        original_safe_name = safe_name
        suffix = 2
        while safe_name in used_names:
            safe_name = f"{original_safe_name}_{suffix}"
            suffix += 1
        used_names.add(safe_name)
        newick_text = newick_text.replace(display_name, safe_name)
        sanitized_names.append({"original": display_name, "label": safe_name})
    if sanitized_names:
        tree_section["sanitized_labels"] = sanitized_names
        tree_section["display_sequence_name"] = sanitized_names[0]["original"]
        tree_section["newick"] = newick_text
        tree_section["sanitized_label"] = sanitized_names[0]["label"]
    query_labels = [item["label"] for item in sanitized_names]
    if query_labels:
        try:
            subtree_text, neighbor_labels = _build_nextclade_nearest_subtree(str(tree_section.get("newick") or newick_text), query_labels, neighbor_limit=50)
            tree_section["newick"] = subtree_text
            tree_section["nearest_neighbor_labels"] = neighbor_labels
            tree_section["query_labels"] = query_labels
            tree_section["subset_neighbor_count"] = len(neighbor_labels)
            tree_section["subset_leaf_count"] = len(query_labels) + len(neighbor_labels)
        except Exception as exc:
            tree_section["itol_error"] = f"精简子树失败: {exc}"
    tree_section["leaf_count"] = len(re.findall(r"(?:(?<=\()|(?<=,))\s*([^():;,]+)\s*:", str(tree_section.get("newick") or "")))
    tree_section["char_count"] = len(str(tree_section.get("newick") or ""))
    return _attach_itol_payload(tree_section, tree_path, str(tree_section.get("newick") or "") or None)
