from __future__ import annotations

from datetime import datetime
from pathlib import Path

from flask import jsonify, request, send_file

from .app_services import get_app_services
from .filesystem_helpers import _assert_within_root, _resolve_browser_path, _resolve_optional_existing_path, _to_browser_path
from .task_manager import ValidationError


def register_batch_input_routes(app) -> None:
    services = get_app_services(app)
    project_root = services.project_root
    login_required = services.access.login_required
    scan_fastq_directory_for_batch_rows = services.batch_inputs.scan_fastq_directory_for_batch_rows

    @app.post("/api/batch-inputs")
    @login_required
    def create_batch_input():
        payload = request.get_json(force=True)
        rows = payload.get("rows", [])
        if not isinstance(rows, list) or not rows:
            raise ValidationError("批量输入不能为空")

        normalized_rows: list[list[str]] = []
        for index, row in enumerate(rows, start=1):
            if not isinstance(row, dict):
                raise ValidationError(f"第 {index} 行格式不正确")
            sample_name = str(row.get("sample_name", "")).strip()
            species = str(row.get("species", "")).strip()
            third_gen = _resolve_optional_existing_path(project_root, row.get("third_gen"))
            short_left = _resolve_optional_existing_path(project_root, row.get("short_left"))
            short_right = _resolve_optional_existing_path(project_root, row.get("short_right"))
            if not sample_name:
                raise ValidationError(f"第 {index} 行样本名称不能为空")
            if not any([third_gen, short_left, short_right]):
                raise ValidationError(f"第 {index} 行至少需要填写一项测序数据")
            normalized_rows.append([sample_name, third_gen, short_left, short_right, species])

        target = services.batch_inputs.write_batch_input(project_root, normalized_rows)

        return jsonify(
            {
                "path": _to_browser_path(project_root.resolve(), target.resolve()),
                "absolute_path": str(target.resolve()),
                "rows": len(normalized_rows),
            }
        ), 201

    @app.post("/api/batch-inputs/fastq-single")
    @login_required
    def create_fastq_single_batch_input():
        payload = request.get_json(force=True)
        task_name = str(payload.get("task_name", "") or "").strip()
        species = "nolevel"
        short_left = _resolve_optional_existing_path(project_root, payload.get("short_left"))
        short_right = _resolve_optional_existing_path(project_root, payload.get("short_right"))
        if not short_left:
            raise ValidationError("请先填写二代左端输入路径")

        left_path = Path(short_left)
        normalized_rows = (
            scan_fastq_directory_for_batch_rows(left_path, species)
            if left_path.is_dir()
            else [[task_name or f"sample_{datetime.now().strftime('%Y%m%d%H%M%S')}", "", short_left, short_right, species]]
        )
        if not normalized_rows:
            raise ValidationError(f"目录中未识别到可组成样本的 fastq 文件: {left_path}")

        target = services.batch_inputs.write_batch_input(project_root, normalized_rows)

        return jsonify(
            {
                "path": _to_browser_path(project_root.resolve(), target.resolve()),
                "absolute_path": str(target.resolve()),
                "rows": len(normalized_rows),
                "scanned_directory": str(left_path) if left_path.is_dir() else "",
            }
        ), 201

    @app.get("/api/batch-inputs/download")
    @login_required
    def download_batch_input():
        batch_path = _resolve_browser_path(project_root.resolve(), request.args.get("path", default="", type=str))
        generated_root = (project_root / "generated_batch_inputs").resolve()
        _assert_within_root(generated_root, batch_path)
        if not batch_path.is_file():
            raise ValidationError(f"批量输入表不存在: {batch_path}")
        return send_file(str(batch_path), mimetype="text/tab-separated-values; charset=utf-8", as_attachment=False)
