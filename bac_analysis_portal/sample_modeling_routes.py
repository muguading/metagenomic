from __future__ import annotations

from flask import jsonify, request, session

from .app_services import get_app_services
from .sample_modeling_service import SampleModelingService
from .task_manager import ValidationError


def register_sample_modeling_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required

    def identity() -> dict[str, str]:
        return {
            "role": str(session.get("role") or ""),
            "username": str(session.get("username") or ""),
            "group_name": str(session.get("group_name") or ""),
        }

    def modeling() -> SampleModelingService:
        return SampleModelingService(store=services.sample_library.store, sample_library=services.sample_library)

    def scope() -> str:
        return str(request.args.get("scope") or "main").strip() or "main"

    def modeling_payload() -> dict:
        payload = request.get_json(silent=True) or {}
        if not isinstance(payload, dict):
            return {}
        return payload

    def modeling_kwargs(payload: dict) -> dict:
        sample_keys = payload.get("sample_keys")
        return {
            "label_field": str(payload.get("label_field") or "").strip(),
            "sample_filter": str(payload.get("sample_filter") or "all").strip() or "all",
            "sample_keys": [str(item or "").strip() for item in sample_keys if str(item or "").strip()] if isinstance(sample_keys, list) else [],
            "algorithm": str(payload.get("algorithm") or "auto_baseline").strip() or "auto_baseline",
        }

    @app.get("/api/database/modeling/readiness")
    @login_required
    def get_database_modeling_readiness():
        return jsonify(modeling().readiness(scope=scope(), **identity()))

    @app.post("/api/database/modeling/readiness")
    @login_required
    def preview_database_modeling_readiness():
        payload = modeling_payload()
        return jsonify(modeling().readiness(scope=scope(), **identity(), **modeling_kwargs(payload)))

    @app.post("/api/database/modeling/label-templates")
    @login_required
    def ensure_database_modeling_label_templates():
        return jsonify(modeling().ensure_label_templates())

    @app.get("/api/database/modeling/dataset")
    @login_required
    def get_database_modeling_dataset():
        return jsonify(modeling().build_dataset(scope=scope(), **identity()))

    @app.post("/api/database/modeling/train-baseline")
    @login_required
    def train_database_modeling_baseline():
        payload = modeling_payload()
        try:
            return jsonify(modeling().train_baseline(scope=scope(), **identity(), **modeling_kwargs(payload)))
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
