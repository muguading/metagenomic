#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


SEMANTIC_HINTS = (
    "icon", "graphic", "chart", "graph", "network", "pictogram", "symbol",
    "screen", "monitor", "dashboard", "device", "people", "user", "shield",
)


def load_items(path: Path) -> list[dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if isinstance(data, list):
        return data
    for key in ("items", "anchors", "visuals"):
        if isinstance(data.get(key), list):
            return data[key]
    raise ValueError("inventory must be a list or contain items/anchors/visuals")


def is_semantic(item: dict[str, Any]) -> bool:
    route = str(item.get("route") or item.get("classification") or "").lower()
    if route == "imagegen_asset":
        return True
    text = " ".join(str(item.get(k, "")) for k in ("id", "anchor_id", "label", "kind", "visual_kind", "notes")).lower()
    return any(hint in text for hint in SEMANTIC_HINTS)


def prompt_for_group(reference: str, items: list[dict[str, Any]], rows: int, cols: int, margin: int, gap: int, bg: str) -> str:
    lines = []
    for idx, item in enumerate(items, 1):
        label = item.get("label") or item.get("semantic_label") or item.get("id") or item.get("anchor_id")
        lines.append(f"{idx}. {label}")
    return (
        "Create an isolated asset grid for a PowerPoint visual replica.\n"
        f"Use the full reference image ({reference}) for global style and the provided crops/residuals for object identity.\n"
        "Objects in order:\n" + "\n".join(lines) + "\n"
        f"Grid: {rows} rows x {cols} columns, margin {margin}px, gap {gap}px; each object centered in its own cell.\n"
        f"Background: perfectly uniform {bg}; no texture, no shadows, no panels, no labels, no arrows, no surrounding slide context.\n"
        "Text: no readable text, no letters, no numbers, no watermark.\n"
        "Output: one clean grid image for deterministic cell cutting."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Create image generation prompt rows for semantic assets.")
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--rows", type=int, default=2)
    parser.add_argument("--cols", type=int, default=4)
    parser.add_argument("--margin", type=int, default=64)
    parser.add_argument("--gap", type=int, default=48)
    parser.add_argument("--background", default="#00FF00")
    parser.add_argument("--provider", default="imagegen", choices=["imagegen", "openai", "gemini"])
    args = parser.parse_args()

    items = [item for item in load_items(args.inventory) if is_semantic(item)]
    rows = []
    capacity = args.rows * args.cols
    for batch_index in range(0, len(items), capacity):
        batch = items[batch_index:batch_index + capacity]
        prompt_id = f"asset_grid_{batch_index // capacity + 1:03d}"
        rows.append({
            "prompt_id": prompt_id,
            "provider": args.provider,
            "prompt_mode": "asset_grid",
            "source_anchor_ids": [str(item.get("anchor_id") or item.get("id")) for item in batch],
            "reference_inputs": (
                [{"role": "full_reference", "path": args.reference}] +
                [{"role": "object_crop", "anchor_id": str(item.get("anchor_id") or item.get("id")), "path": item.get("source_crop", "")} for item in batch]
            ),
            "grid": {"rows": args.rows, "cols": args.cols, "margin": args.margin, "gap": args.gap, "background": args.background},
            "prompt": prompt_for_group(args.reference, batch, args.rows, args.cols, args.margin, args.gap, args.background),
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(json.dumps(row, ensure_ascii=False) for row in rows) + ("\n" if rows else ""), encoding="utf-8")
    print(json.dumps({"prompt_rows": len(rows), "semantic_items": len(items), "out": str(args.out)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
