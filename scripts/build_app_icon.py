from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image, ImageChops, ImageFilter


ROOT = Path(__file__).resolve().parents[1]
STATIC_DIR = ROOT / "bac_analysis_portal" / "static"
SOURCE_MARK_SVG = STATIC_DIR / "image_v1_mark.svg"
PNG_PATH = STATIC_DIR / "app_icon.png"
ICNS_PATH = STATIC_DIR / "app_icon.icns"
ICONSET_DIR = STATIC_DIR / "app_icon.iconset"


ICONSET_SIZES = {
    "icon_16x16.png": 16,
    "icon_16x16@2x.png": 32,
    "icon_32x32.png": 32,
    "icon_32x32@2x.png": 64,
    "icon_128x128.png": 128,
    "icon_128x128@2x.png": 256,
    "icon_256x256.png": 256,
    "icon_256x256@2x.png": 512,
    "icon_512x512.png": 512,
    "icon_512x512@2x.png": 1024,
}


def _run(command: list[str], cwd: Path | None = None) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def _render_svg_to_png(svg_path: Path, output_png: Path, size: int = 1024) -> None:
    with tempfile.TemporaryDirectory(prefix="app_icon_render_") as temp_dir:
        temp_path = Path(temp_dir)
        _run(["qlmanage", "-t", "-s", str(size), "-o", str(temp_path), str(svg_path)])
        rendered = temp_path / f"{svg_path.name}.png"
        if not rendered.exists():
            raise FileNotFoundError(f"Quick Look did not produce {rendered}")
        shutil.copy2(rendered, output_png)


def _remove_white_background(image: Image.Image, threshold: int = 247) -> Image.Image:
    rgba = image.convert("RGBA")
    result = Image.new("RGBA", rgba.size, (0, 0, 0, 0))
    src = rgba.load()
    dst = result.load()
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = src[x, y]
            if a == 0:
                continue
            brightness = min(r, g, b)
            if brightness >= threshold:
                alpha = max(0, int(a * (255 - brightness) / max(1, 255 - threshold)))
                if alpha == 0:
                    continue
                dst[x, y] = (r, g, b, alpha)
            else:
                dst[x, y] = (r, g, b, a)
    return result


def _content_bbox(image: Image.Image) -> tuple[int, int, int, int]:
    alpha = image.getchannel("A")
    bbox = alpha.getbbox()
    if not bbox:
        raise ValueError("No visible content found in source mark")
    return bbox


def build_icon(size: int = 1024, mark_coverage: float = 0.84) -> Image.Image:
    with tempfile.TemporaryDirectory(prefix="app_icon_source_") as temp_dir:
        temp_png = Path(temp_dir) / "source.png"
        _render_svg_to_png(SOURCE_MARK_SVG, temp_png, size=size)
        source = Image.open(temp_png).convert("RGBA")

    source = _remove_white_background(source)
    source = source.crop(_content_bbox(source))

    shadow = Image.new("RGBA", source.size, (0, 0, 0, 0))
    shadow.alpha_composite(source, (0, 0))
    shadow = shadow.filter(ImageFilter.GaussianBlur(radius=max(6, size // 80)))
    shadow_mask = shadow.getchannel("A").point(lambda value: int(value * 0.18))
    shadow.putalpha(shadow_mask)

    target_span = int(size * mark_coverage)
    scale = min(target_span / source.width, target_span / source.height)
    target_size = (
        max(1, int(round(source.width * scale))),
        max(1, int(round(source.height * scale))),
    )

    resized_mark = source.resize(target_size, Image.Resampling.LANCZOS)
    resized_shadow = shadow.resize(target_size, Image.Resampling.LANCZOS)

    canvas = Image.new("RGBA", (size, size), (255, 255, 255, 255))
    left = (size - target_size[0]) // 2
    top = (size - target_size[1]) // 2 - int(size * 0.015)

    shadow_left = left
    shadow_top = top + int(size * 0.022)
    canvas.alpha_composite(resized_shadow, (shadow_left, shadow_top))
    canvas.alpha_composite(resized_mark, (left, top))
    return canvas


def save_icon_variants() -> None:
    STATIC_DIR.mkdir(parents=True, exist_ok=True)
    icon = build_icon(1024)
    icon.save(PNG_PATH)

    if ICONSET_DIR.exists():
        shutil.rmtree(ICONSET_DIR)
    ICONSET_DIR.mkdir(parents=True, exist_ok=True)

    for filename, icon_size in ICONSET_SIZES.items():
        icon.resize((icon_size, icon_size), Image.Resampling.LANCZOS).save(ICONSET_DIR / filename)

    if ICNS_PATH.exists():
        ICNS_PATH.unlink()
    _run(["iconutil", "-c", "icns", str(ICONSET_DIR), "-o", str(ICNS_PATH)])


if __name__ == "__main__":
    save_icon_variants()
