from __future__ import annotations

import json
import os
import sys
import zipfile
from dataclasses import asdict
from datetime import datetime
from functools import wraps
from pathlib import Path
from uuid import uuid4

from flask import Flask, abort, jsonify, redirect, render_template, request, session, url_for
from sqlalchemy import select
from werkzeug.utils import secure_filename

if __package__ in (None, ""):
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from genome_db.auth_manager import AuthenticationError, AuthorizationError, UserManager
    from genome_db.database import DEFAULT_DB_URL
    from genome_db.bulk_import import EXPECTED_COLUMNS, BulkImportFormatError, parse_bulk_import_rows
    from genome_db.genome_manager import DuplicateGenomeError, GenomeManager
    from genome_db.models import Genome
    from genome_db.validators import ValidationError, validate_fasta_file
else:
    from .auth_manager import AuthenticationError, AuthorizationError, UserManager
    from .database import DEFAULT_DB_URL
    from .bulk_import import EXPECTED_COLUMNS, BulkImportFormatError, parse_bulk_import_rows
    from .genome_manager import DuplicateGenomeError, GenomeManager
    from .models import Genome
    from .validators import ValidationError, validate_fasta_file


ALLOWED_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas"}
ALLOWED_BULK_SUFFIXES = {".xlsx", ".csv", ".tsv"}


def create_app(database_url: str = DEFAULT_DB_URL) -> Flask:
    base_dir = Path(__file__).resolve().parent
    upload_dir = base_dir / "uploads"
    metadata_upload_dir = upload_dir / "metadata"
    demo_dir = base_dir / "static" / "demo"
    upload_dir.mkdir(exist_ok=True)
    metadata_upload_dir.mkdir(exist_ok=True)
    demo_dir.mkdir(parents=True, exist_ok=True)
    app = Flask(
        __name__,
        template_folder=str(base_dir / "templates"),
        static_folder=str(base_dir / "static"),
    )
    app.secret_key = os.environ.get("GENOME_DB_SECRET_KEY", "genome-db-dev-secret-key")
    app.config["UPLOAD_FOLDER"] = str(upload_dir)
    app.config["SESSION_COOKIE_HTTPONLY"] = True
    app.config["SESSION_COOKIE_SAMESITE"] = "Lax"

    manager = GenomeManager(database_url=database_url)
    user_manager = UserManager(manager.session_factory)
    _ensure_bulk_import_demo_files(demo_dir)

    def cleanup_unreferenced_upload(genome_file_path: str) -> None:
        try:
            path = Path(genome_file_path).expanduser().resolve()
        except Exception:
            return
        if upload_dir not in path.parents or not path.is_file():
            return
        with manager.session_factory() as db_session:
            existing = db_session.scalar(select(Genome).where(Genome.genome_file_path == str(path)))
            if existing is not None:
                return
        path.unlink(missing_ok=True)

    def cleanup_unreferenced_metadata_files(custom_metadata: object) -> None:
        if not isinstance(custom_metadata, list):
            return
        for item in custom_metadata:
            if not isinstance(item, dict) or item.get("type") != "file":
                continue
            value = str(item.get("value", "")).strip()
            if not value:
                continue
            try:
                path = Path(value).expanduser().resolve()
            except Exception:
                continue
            if metadata_upload_dir not in path.parents or not path.is_file():
                continue
            if manager.is_metadata_file_referenced(str(path)):
                continue
            path.unlink(missing_ok=True)

    @app.context_processor
    def inject_session_user():
        return {
            "current_user": {
                "username": session.get("username"),
                "role": session.get("role"),
            }
            if session.get("username")
            else None
        }

    @app.before_request
    def protect_app():
        if request.endpoint in {"login", "login_post", "register", "register_post", "static"}:
            return None
        if not _is_logged_in():
            if request.path.startswith("/api/"):
                return jsonify({"error": "Authentication required"}), 401
            return redirect(url_for("login"))
        return None

    @app.get("/login")
    def login():
        if _is_logged_in():
            return redirect(url_for("index"))
        return render_template("login.html")

    @app.post("/login")
    def login_post():
        payload = request.get_json(force=True) if request.is_json else request.form
        user = user_manager.authenticate(payload.get("username", ""), payload.get("password", ""))
        session["username"] = user.username
        session["role"] = user.role
        return jsonify({"status": "ok", "user": {"username": user.username, "role": user.role}})

    @app.post("/register")
    def register_post():
        payload = request.get_json(force=True) if request.is_json else request.form
        if payload.get("password", "") != payload.get("confirm_password", ""):
            raise ValidationError("Passwords do not match")
        created = user_manager.register_user(
            username=payload.get("username", ""),
            password=payload.get("password", ""),
            display_name=payload.get("display_name"),
            email=payload.get("email"),
        )
        return jsonify(created), 201

    @app.post("/logout")
    @login_required
    def logout():
        session.clear()
        return jsonify({"status": "logged_out"})

    @app.get("/")
    @login_required
    def index():
        return render_template("index.html", database_url=database_url)

    @app.get("/api/session")
    @login_required
    def current_session():
        return jsonify(user_manager.get_user_profile(session["username"]))

    @app.get("/api/profile")
    @login_required
    def get_profile():
        return jsonify(user_manager.get_user_profile(session["username"]))

    @app.put("/api/profile")
    @login_required
    def update_profile():
        payload = request.get_json(force=True)
        updated = user_manager.update_user_profile(
            username=session["username"],
            new_username=payload.get("username", session["username"]),
            display_name=payload.get("display_name"),
            email=payload.get("email"),
        )
        session["username"] = updated["username"]
        return jsonify(updated)

    @app.post("/api/profile/password")
    @login_required
    def change_profile_password():
        payload = request.get_json(force=True)
        if payload.get("new_password", "") != payload.get("confirm_password", ""):
            raise ValidationError("Passwords do not match")
        user_manager.change_password(
            username=session["username"],
            current_password=payload.get("current_password", ""),
            new_password=payload.get("new_password", ""),
        )
        return jsonify({"status": "password_updated"})

    @app.get("/api/health")
    @login_required
    def health():
        return jsonify({"status": "ok", "database_url": database_url})

    @app.get("/api/users")
    @admin_required
    def list_users():
        return jsonify({"items": user_manager.list_users()})

    @app.post("/api/users")
    @admin_required
    def create_user():
        payload = request.get_json(force=True)
        created = user_manager.create_user(
            username=payload["username"],
            password=payload["password"],
            role=payload["role"],
        )
        return jsonify(created), 201

    @app.put("/api/users/<username>/password")
    @admin_required
    def reset_user_password(username: str):
        payload = request.get_json(force=True)
        updated = user_manager.admin_reset_password(
            target_username=username,
            new_password=payload.get("new_password", ""),
        )
        return jsonify(updated)

    @app.post("/api/upload-fasta")
    @login_required
    def upload_fasta():
        if "file" not in request.files:
            raise ValidationError("No file part found in upload request")

        file_storage = request.files["file"]
        if not file_storage.filename:
            raise ValidationError("No FASTA file selected")

        original_name = secure_filename(file_storage.filename)
        suffix = Path(original_name).suffix.lower()
        if suffix not in ALLOWED_FASTA_SUFFIXES:
            raise ValidationError("Only .fa, .fasta, .fna, or .fas files are supported")

        unique_prefix = datetime.utcnow().strftime("%Y%m%d%H%M%S")
        target_name = f"{unique_prefix}_{uuid4().hex[:8]}_{original_name}"
        target_path = upload_dir / target_name
        file_storage.save(target_path)
        try:
            validated_path, genome_length = validate_fasta_file(str(target_path))
        except Exception:
            if target_path.exists():
                target_path.unlink()
            raise

        return jsonify(
            {
                "filename": original_name,
                "stored_path": validated_path,
                "genome_length": genome_length,
            }
        ), 201

    @app.post("/api/upload-metadata-file")
    @login_required
    def upload_metadata_file():
        if "file" not in request.files:
            raise ValidationError("No file part found in upload request")

        file_storage = request.files["file"]
        if not file_storage.filename:
            raise ValidationError("No file selected")

        original_name = secure_filename(file_storage.filename)
        unique_prefix = datetime.utcnow().strftime("%Y%m%d%H%M%S")
        target_name = f"{unique_prefix}_{uuid4().hex[:8]}_{original_name}"
        target_path = metadata_upload_dir / target_name
        file_storage.save(target_path)
        return jsonify({"filename": original_name, "stored_path": str(target_path.resolve())}), 201

    @app.post("/api/bulk-import-genomes")
    @login_required
    def bulk_import_genomes():
        if "file" not in request.files:
            raise ValidationError("No bulk import file provided")

        file_storage = request.files["file"]
        if not file_storage.filename:
            raise ValidationError("No import file selected")

        original_name = secure_filename(file_storage.filename)
        suffix = Path(original_name).suffix.lower()
        if suffix not in ALLOWED_BULK_SUFFIXES:
            raise ValidationError("Only .xlsx, .csv, and .tsv files are supported")

        rows = parse_bulk_import_rows(original_name, file_storage.read())
        results: list[dict[str, object]] = []
        success_count = 0
        failure_count = 0
        for index, row in enumerate(rows, start=2):
            location = {
                "province": row.get("location_province", ""),
                "city": row.get("location_city", ""),
                "district": row.get("location_district", ""),
                "detail": row.get("location_detail", ""),
            }
            province = str(location["province"] or "").strip()
            city = str(location["city"] or "").strip()
            district = str(location["district"] or "").strip()
            normalized_location = location if all([province, city, district]) else None
            try:
                created = manager.add_genome(
                    genome_id=row.get("genome_id", ""),
                    sample_name=row.get("sample_name", ""),
                    species_name=row.get("species_name", ""),
                    taxid=int(row.get("taxid", "0") or 0),
                    genome_file_path=row.get("genome_file_path", ""),
                    submitter=session["username"],
                    gender=row.get("gender"),
                    country=row.get("country"),
                    location=normalized_location,
                    collection_time=row.get("collection_time"),
                    sample_type=row.get("sample_type"),
                    sequencing_method=row.get("sequencing_method"),
                    description=row.get("description"),
                    custom_metadata=[],
                )
                results.append({"row_number": index, "genome_id": created.get("genome_id"), "status": "SUCCESS"})
                success_count += 1
            except Exception as exc:
                failure_count += 1
                cleanup_unreferenced_upload(row.get("genome_file_path", ""))
                results.append(
                    {
                        "row_number": index,
                        "genome_id": row.get("genome_id", ""),
                        "status": "FAILED",
                        "error": str(exc),
                    }
                )
        return jsonify(
            {
                "total": len(rows),
                "success_count": success_count,
                "failure_count": failure_count,
                "results": results,
            }
        )

    @app.get("/api/genomes")
    @login_required
    def search_genomes():
        submitter = request.args.get("submitter")
        if session.get("role") != "admin" and submitter in (None, ""):
            submitter = session["username"]
        raw_custom_filters = request.args.get("custom_filters")
        try:
            custom_filters = json.loads(raw_custom_filters) if raw_custom_filters else None
        except json.JSONDecodeError as exc:
            raise ValidationError("custom_filters must be valid JSON") from exc
        result = manager.search_genomes(
            species_name=request.args.get("species_name"),
            taxid=_parse_optional_int(request.args.get("taxid")),
            submitter=submitter,
            custom_logic=request.args.get("custom_logic", "and"),
            custom_filters=custom_filters,
            page=_parse_int(request.args.get("page"), default=1),
            page_size=_parse_int(request.args.get("page_size"), default=10),
            operator=session["username"],
        )
        return jsonify(asdict(result))

    @app.get("/api/metadata-templates")
    @login_required
    def list_metadata_templates():
        requested_submitter = request.args.get("submitter")
        submitter = requested_submitter.strip() if requested_submitter else None
        if session.get("role") != "admin":
            submitter = session["username"]
        elif not submitter:
            submitter = session["username"]
        return jsonify({"items": manager.list_metadata_templates(submitter=submitter)})

    @app.get("/api/dashboard-data")
    @login_required
    def get_dashboard_data():
        requested_submitter = request.args.get("submitter")
        submitter = requested_submitter.strip() if requested_submitter else None
        if session.get("role") != "admin":
            submitter = session["username"]
        result = manager.get_dashboard_data(submitter=submitter, operator=session["username"])
        return jsonify({"items": result.items, "templates": result.templates})

    @app.get("/api/genomes/<int:record_id>")
    @login_required
    def get_genome(record_id: int):
        genome = manager.get_genome(record_id, operator=session["username"])
        _ensure_genome_access(genome)
        return jsonify(genome)

    @app.post("/api/genomes")
    @login_required
    def add_genome():
        payload = request.get_json(force=True)
        genome_file_path = payload["genome_file_path"]
        custom_metadata = payload.get("custom_metadata", [])
        try:
            created = manager.add_genome(
                genome_id=payload["genome_id"],
                sample_name=payload["sample_name"],
                species_name=payload["species_name"],
                taxid=int(payload["taxid"]),
                genome_file_path=genome_file_path,
                submitter=session["username"],
                description=payload.get("description"),
                gender=payload.get("gender"),
                country=payload.get("country"),
                location=payload.get("location"),
                collection_time=payload.get("collection_time"),
                sample_type=payload.get("sample_type"),
                sequencing_method=payload.get("sequencing_method"),
                custom_metadata=custom_metadata,
            )
        except Exception:
            cleanup_unreferenced_upload(genome_file_path)
            cleanup_unreferenced_metadata_files(custom_metadata)
            raise
        return jsonify(created), 201

    @app.put("/api/genomes/<int:record_id>")
    @login_required
    def update_genome(record_id: int):
        current = manager.get_genome(record_id, operator=session["username"])
        _ensure_genome_access(current, write=True)
        payload = request.get_json(force=True)
        updates = {"operator": session["username"], "submitter": current["submitter"]}
        uploaded_path = payload.get("genome_file_path")
        custom_metadata = payload.get("custom_metadata", current.get("custom_metadata", []))
        for field in (
            "sample_name",
            "species_name",
            "taxid",
            "genome_file_path",
            "description",
            "gender",
            "country",
            "location",
            "collection_time",
            "sample_type",
            "sequencing_method",
        ):
            value = payload.get(field)
            if value is not None:
                updates[field] = value
        updates["custom_metadata"] = custom_metadata
        try:
            updated = manager.update_genome(record_id, **updates)
        except Exception:
            if uploaded_path and uploaded_path != current.get("genome_file_path"):
                cleanup_unreferenced_upload(uploaded_path)
            cleanup_unreferenced_metadata_files(custom_metadata)
            raise
        return jsonify({k: v for k, v in updated.items() if v is not None})

    @app.delete("/api/genomes/<int:record_id>")
    @login_required
    def delete_genome(record_id: int):
        current = manager.get_genome(record_id, operator=session["username"])
        _ensure_genome_access(current, write=True)
        manager.delete_genome(record_id, operator=session["username"])
        return jsonify({"status": "deleted", "genome_id": current["genome_id"], "id": record_id})

    @app.get("/api/audit-logs")
    @login_required
    def list_audit_logs():
        genome_id = request.args.get("genome_id")
        operation = request.args.get("operation")
        result = manager.list_audit_logs(
            genome_id=genome_id,
            operation=operation,
            page=_parse_int(request.args.get("page"), default=1),
            page_size=_parse_int(request.args.get("page_size"), default=10),
            operator=session["username"],
        )
        payload = asdict(result)
        if session.get("role") != "admin":
            payload["items"] = [item for item in payload["items"] if item.get("operator") == session["username"]]
            payload["total"] = len(payload["items"])
            payload["pages"] = 1 if payload["items"] else 0
            payload["page"] = 1
        return jsonify(payload)

    @app.errorhandler(ValidationError)
    @app.errorhandler(DuplicateGenomeError)
    @app.errorhandler(AuthenticationError)
    @app.errorhandler(BulkImportFormatError)
    def handle_validation_error(error: Exception):
        return jsonify({"error": str(error)}), 400

    @app.errorhandler(AuthorizationError)
    def handle_authorization_error(error: AuthorizationError):
        return jsonify({"error": str(error)}), 403

    @app.errorhandler(KeyError)
    def handle_not_found(error: KeyError):
        return jsonify({"error": str(error)}), 404

    @app.errorhandler(Exception)
    def handle_unknown_error(error: Exception):
        if request.path.startswith("/api/"):
            app.logger.exception("Unhandled genome database UI error")
            return jsonify({"error": str(error)}), 500
        raise error

    def _ensure_genome_access(genome: dict[str, object], write: bool = False) -> None:
        if session.get("role") == "admin":
            return
        if genome.get("submitter") != session.get("username"):
            raise AuthorizationError("You do not have permission to access this genome record")
        if write and genome.get("submitter") != session.get("username"):
            raise AuthorizationError("You do not have permission to modify this genome record")

    return app


def login_required(view_func):
    @wraps(view_func)
    def wrapped(*args, **kwargs):
        if not _is_logged_in():
            if request.path.startswith("/api/"):
                return jsonify({"error": "Authentication required"}), 401
            return redirect(url_for("login"))
        return view_func(*args, **kwargs)

    return wrapped


def admin_required(view_func):
    @wraps(view_func)
    def wrapped(*args, **kwargs):
        if not _is_logged_in():
            return jsonify({"error": "Authentication required"}), 401
        if session.get("role") != "admin":
            raise AuthorizationError("Administrator permission required")
        return view_func(*args, **kwargs)

    return wrapped


def _is_logged_in() -> bool:
    return bool(session.get("username") and session.get("role"))


def _parse_int(raw: str | None, *, default: int) -> int:
    if raw in (None, ""):
        return default
    return int(raw)


def _parse_optional_int(raw: str | None) -> int | None:
    if raw in (None, ""):
        return None
    return int(raw)


def _ensure_bulk_import_demo_files(demo_dir: Path) -> None:
    header = EXPECTED_COLUMNS
    row = [
        "demo_chikv_001",
        "CHIKV-Demo-001",
        "Chikungunya virus",
        "37124",
        "/path/to/your/genome_demo_001.fna",
        "Female",
        "中国",
        "上海市",
        "上海市",
        "黄浦区",
        "People's Square sentinel site",
        "2026-03-18T09:30",
        "Serum",
        "Illumina",
        "Replace genome_file_path with a valid server-side FASTA path before import.",
    ]
    csv_content = ",".join(_csv_escape(item) for item in header) + "\n" + ",".join(_csv_escape(item) for item in row) + "\n"
    tsv_content = "\t".join(header) + "\n" + "\t".join(row) + "\n"
    (demo_dir / "genome_bulk_import_demo.csv").write_text(csv_content, encoding="utf-8")
    (demo_dir / "genome_bulk_import_demo.tsv").write_text(tsv_content, encoding="utf-8")
    _write_simple_xlsx(demo_dir / "genome_bulk_import_demo.xlsx", [header, row])


def _csv_escape(value: str) -> str:
    text = str(value)
    if any(token in text for token in [",", "\"", "\n", "\r", "\t"]):
        escaped = text.replace('"', '""')
        return f'"{escaped}"'
    return text


def _write_simple_xlsx(path: Path, rows: list[list[str]]) -> None:
    shared_strings: list[str] = []
    shared_index: dict[str, int] = {}

    def shared_id(value: str) -> int:
        if value not in shared_index:
            shared_index[value] = len(shared_strings)
            shared_strings.append(value)
        return shared_index[value]

    worksheet_rows: list[str] = []
    for row_index, row in enumerate(rows, start=1):
        cells: list[str] = []
        for col_index, value in enumerate(row, start=1):
            reference = f"{_excel_column_name(col_index)}{row_index}"
            cells.append(f'<c r="{reference}" t="s"><v>{shared_id(str(value))}</v></c>')
        worksheet_rows.append(f'<row r="{row_index}">{"".join(cells)}</row>')

    shared_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        f'count="{len(shared_strings)}" uniqueCount="{len(shared_strings)}">'
        + "".join(f"<si><t>{_xml_escape(value)}</t></si>" for value in shared_strings)
        + "</sst>"
    )
    sheet_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(worksheet_rows)}</sheetData>'
        "</worksheet>"
    )
    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="ImportDemo" sheetId="1" r:id="rId1"/></sheets>'
        "</workbook>"
    )
    workbook_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
        '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/sharedStrings" Target="sharedStrings.xml"/>'
        "</Relationships>"
    )
    root_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>'
        "</Types>"
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types_xml)
        archive.writestr("_rels/.rels", root_rels_xml)
        archive.writestr("xl/workbook.xml", workbook_xml)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        archive.writestr("xl/worksheets/sheet1.xml", sheet_xml)
        archive.writestr("xl/sharedStrings.xml", shared_xml)


def _excel_column_name(index: int) -> str:
    result = ""
    current = index
    while current:
        current, remainder = divmod(current - 1, 26)
        result = chr(65 + remainder) + result
    return result


def _xml_escape(value: str) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


if __name__ == "__main__":
    app = create_app()
    app.run(host="127.0.0.1", port=5000, debug=True)
