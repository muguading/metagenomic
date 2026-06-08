#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from PIL import Image


def parse_color(value: str) -> tuple[int, int, int]:
    value = value.strip()
    if value.startswith("#"):
        value = value[1:]
        return int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16)
    parts = [int(p) for p in value.split(",")]
    if len(parts) != 3:
        raise ValueError("background color must be #RRGGBB or r,g,b")
    return parts[0], parts[1], parts[2]


def remove_background(img: Image.Image, key: tuple[int, int, int], tolerance: int, dominance: int) -> Image.Image:
    out = img.convert("RGBA")
    px = out.load()
    for y in range(out.height):
        for x in range(out.width):
            r, g, b, a = px[x, y]
            dist = math.sqrt((r - key[0]) ** 2 + (g - key[1]) ** 2 + (b - key[2]) ** 2)
            keyed = dist <= tolerance
            if key[1] > key[0] and key[1] > key[2]:
                keyed = keyed or (g > 120 and g - max(r, b) >= dominance)
            if keyed:
                px[x, y] = (r, g, b, 0)
    bbox = out.getchannel("A").getbbox()
    return out.crop(bbox) if bbox else out


def load_spec(path: Path | None, rows: int, cols: int, margin: int, gap: int, ids: list[str]) -> dict[str, Any]:
    if path:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    return {"rows": rows, "cols": cols, "margin": margin, "gap": gap, "ids": ids}


def cells_from_spec(spec: dict[str, Any], width: int, height: int) -> list[dict[str, Any]]:
    if isinstance(spec.get("cells"), list):
        return spec["cells"]
    rows = int(spec["rows"])
    cols = int(spec["cols"])
    margin = int(spec.get("margin", 0))
    gap = int(spec.get("gap", 0))
    ids = list(spec.get("ids") or [])
    cell_w = (width - margin * 2 - gap * (cols - 1)) / cols
    cell_h = (height - margin * 2 - gap * (rows - 1)) / rows
    cells = []
    for r in range(rows):
        for c in range(cols):
            idx = r * cols + c
            cells.append({
                "id": ids[idx] if idx < len(ids) else f"cell_{idx + 1:03d}",
                "bbox": [round(margin + c * (cell_w + gap)), round(margin + r * (cell_h + gap)), round(cell_w), round(cell_h)],
            })
    return cells


def main() -> None:
    parser = argparse.ArgumentParser(description="Cut generated asset grids and transparentize backgrounds.")
    parser.add_argument("--image", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--manifest-out", type=Path)
    parser.add_argument("--spec", type=Path)
    parser.add_argument("--rows", type=int, default=1)
    parser.add_argument("--cols", type=int, default=1)
    parser.add_argument("--margin", type=int, default=0)
    parser.add_argument("--gap", type=int, default=0)
    parser.add_argument("--ids", nargs="*", default=[])
    parser.add_argument("--background", default="#00FF00")
    parser.add_argument("--tolerance", type=int, default=45)
    parser.add_argument("--dominance", type=int, default=35)
    args = parser.parse_args()

    src = Image.open(args.image).convert("RGBA")
    spec = load_spec(args.spec, args.rows, args.cols, args.margin, args.gap, args.ids)
    key = parse_color(str(spec.get("background") or args.background))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for cell in cells_from_spec(spec, src.width, src.height):
        x, y, w, h = [int(round(float(v))) for v in cell["bbox"]]
        asset_id = str(cell.get("id") or f"asset_{len(records) + 1:03d}")
        crop = src.crop((x, y, x + w, y + h))
        final = remove_background(crop, key, args.tolerance, args.dominance)
        out = args.out_dir / f"{asset_id}.png"
        final.save(out)
        records.append({"id": asset_id, "source_grid": str(args.image), "cell_bbox": [x, y, w, h], "path": str(out), "size": [final.width, final.height]})
    if args.manifest_out:
        args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
        args.manifest_out.write_text(json.dumps({"assets": records}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"assets": len(records), "out_dir": str(args.out_dir)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
