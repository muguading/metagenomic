from __future__ import annotations

from flask import Response, jsonify, request, session

from .app_services import get_app_services
from .modeling_service import ModelingPlatformService
from .task_manager import ValidationError


def register_modeling_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required

    def identity() -> dict[str, str]:
        return {
            "role": str(session.get("role") or ""),
            "username": str(session.get("username") or ""),
            "group_name": str(session.get("group_name") or ""),
        }

    def payload() -> dict:
        data = request.get_json(silent=True) or {}
        return data if isinstance(data, dict) else {}

    def scope() -> str:
        return str(request.args.get("scope") or payload().get("scope") or "main").strip() or "main"

    def service() -> ModelingPlatformService:
        return ModelingPlatformService(
            store=services.sample_library.store,
            sample_library=services.sample_library,
            project_root=services.project_root,
        )

    @app.get("/api/modeling/options")
    @login_required
    def get_modeling_options():
        return jsonify(service().options())

    @app.get("/api/modeling/samples/search")
    @login_required
    def search_modeling_samples():
        filters = {key: value for key, value in request.args.items() if key != "scope"}
        return jsonify(service().search_samples(scope=scope(), filters=filters, **identity()))

    @app.post("/api/modeling/samples/search")
    @login_required
    def post_search_modeling_samples():
        data = payload()
        filters = data.get("filters") if isinstance(data.get("filters"), dict) else data
        return jsonify(service().search_samples(scope=scope(), filters=filters, **identity()))

    @app.post("/api/modeling/datasets")
    @login_required
    def create_modeling_dataset():
        return jsonify(service().create_dataset(scope=scope(), payload=payload(), **identity()))

    @app.get("/api/modeling/datasets")
    @login_required
    def list_modeling_datasets():
        return jsonify(service().list_datasets())

    @app.get("/api/modeling/datasets/<path:dataset_id>")
    @login_required
    def get_modeling_dataset(dataset_id: str):
        return jsonify(service().get_dataset(dataset_id))

    @app.post("/api/modeling/feature-sets")
    @login_required
    def create_modeling_feature_set():
        return jsonify(service().create_feature_set(username=identity()["username"], payload=payload()))

    @app.get("/api/modeling/feature-sets")
    @login_required
    def list_modeling_feature_sets():
        return jsonify(service().list_feature_sets())

    @app.post("/api/modeling/train")
    @login_required
    def train_modeling_model():
        try:
            result = service().train(scope=scope(), payload=payload(), **identity())
        except PermissionError as exc:
            raise ValidationError(str(exc)) from exc
        if result.get("error"):
            raise ValidationError(str(result["error"]))
        return jsonify(result)

    @app.get("/api/modeling/train-jobs/<path:job_id>")
    @login_required
    def get_modeling_train_job(job_id: str):
        return jsonify(service().get_train_job(job_id))

    @app.get("/api/modeling/models")
    @login_required
    def list_modeling_models():
        return jsonify(service().list_models())

    @app.get("/api/modeling/models/<path:model_id>")
    @login_required
    def get_modeling_model(model_id: str):
        return jsonify(service().get_model(model_id))

    @app.get("/api/modeling/models/<path:model_id>/report")
    @login_required
    def export_modeling_model_report(model_id: str):
        html = service().export_model_report(model_id)
        return Response(
            html,
            mimetype="text/html; charset=utf-8",
            headers={"Content-Disposition": f'attachment; filename="model-evaluation-{model_id.replace(":", "-")}.html"'},
        )

    @app.post("/api/modeling/models/<path:model_id>/activate")
    @login_required
    def activate_modeling_model(model_id: str):
        try:
            return jsonify(service().set_model_status(model_id, "active", role=identity()["role"]))
        except PermissionError as exc:
            raise ValidationError(str(exc)) from exc

    @app.post("/api/modeling/models/<path:model_id>/deactivate")
    @login_required
    def deactivate_modeling_model(model_id: str):
        try:
            return jsonify(service().set_model_status(model_id, "inactive", role=identity()["role"]))
        except PermissionError as exc:
            raise ValidationError(str(exc)) from exc

    @app.delete("/api/modeling/models/<path:model_id>")
    @login_required
    def delete_modeling_model(model_id: str):
        try:
            return jsonify(service().delete_model(model_id, role=identity()["role"]))
        except PermissionError as exc:
            raise ValidationError(str(exc)) from exc

    @app.post("/api/modeling/models/<path:model_id>/predict")
    @login_required
    def predict_modeling_samples(model_id: str):
        try:
            return jsonify(service().predict(model_id=model_id, scope=scope(), payload=payload(), **identity()))
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc

    @app.get("/api/modeling/predictions")
    @login_required
    def list_modeling_predictions():
        return jsonify(service().list_predictions())

    @app.get("/api/modeling/predictions/<path:prediction_id>")
    @login_required
    def get_modeling_prediction(prediction_id: str):
        return jsonify(service().get_prediction(prediction_id))
