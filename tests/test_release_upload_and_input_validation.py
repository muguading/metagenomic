from __future__ import annotations

import io

import pytest

from conftest import login_as


pytestmark = pytest.mark.release


def test_release_fastq_directory_batch_generation_supports_multi_sample_detection(release_app_client, release_data) -> None:
    _app, client = release_app_client
    headers = login_as(client)

    response = client.post(
        "/api/batch-inputs/fastq-single",
        json={"task_name": "release_fastq_batch", "short_left": str(release_data["input_dir"]), "short_right": ""},
        headers=headers,
    )

    assert response.status_code == 201
    payload = response.get_json()
    assert payload["rows"] == 1
    assert payload["path"].endswith(".tsv")


def test_release_fastq_directory_batch_generation_rejects_empty_directory(release_app_client, tmp_path) -> None:
    _app, client = release_app_client
    headers = login_as(client)
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()

    response = client.post(
        "/api/batch-inputs/fastq-single",
        json={"task_name": "release_empty_fastq", "short_left": str(empty_dir), "short_right": ""},
        headers=headers,
    )

    assert response.status_code == 400
    assert "未识别到可组成样本的 fastq 文件" in response.get_json()["error"]


def test_release_task_reference_upload_accepts_fasta_and_rejects_wrong_extension(release_app_client) -> None:
    _app, client = release_app_client
    headers = login_as(client)

    accepted = client.post(
        "/api/task-references/upload",
        data={"file": (io.BytesIO(b">ref\nACGT\n"), "ref.fa")},
        headers=headers,
        content_type="multipart/form-data",
    )
    rejected = client.post(
        "/api/task-references/upload",
        data={"file": (io.BytesIO(b"not fasta"), "ref.txt")},
        headers=headers,
        content_type="multipart/form-data",
    )

    assert accepted.status_code == 200
    assert accepted.get_json()["path"].endswith(".fa")
    assert rejected.status_code == 400
    assert "仅支持 .fa / .fasta / .fna" in rejected.get_json()["error"]


def test_release_sample_batch_precheck_accepts_valid_csv_and_rejects_missing_columns(
    release_app_client, release_data
) -> None:
    _app, client = release_app_client
    headers = login_as(client)
    valid_csv = (
        "sample_name,final_fasta_path,case_id,surveillance_source,suspected_syndrome,"
        "specimen_category,province,city,district,location_detail\n"
        f"S1,{release_data['fasta']},CASE-001,门急诊,肺炎/呼吸道感染,呼吸道,上海市,上海市,黄浦区,门诊采样点\n"
    ).encode("utf-8")
    invalid_csv = b"sample_name\nS1\n"

    valid = client.post(
        "/api/database/batch-import/precheck",
        data={"file": (io.BytesIO(valid_csv), "samples.csv")},
        headers=headers,
        content_type="multipart/form-data",
    )
    invalid = client.post(
        "/api/database/batch-import/precheck",
        data={"file": (io.BytesIO(invalid_csv), "samples.csv")},
        headers=headers,
        content_type="multipart/form-data",
    )

    assert valid.status_code == 200
    assert valid.get_json()["can_import"] is True
    assert invalid.status_code == 200
    assert invalid.get_json()["can_import"] is False
    assert "final_fasta_path" in invalid.get_json()["missing_columns"]


def test_release_reference_batch_precheck_accepts_valid_tsv_and_rejects_empty_file(
    release_app_client, release_data
) -> None:
    _app, client = release_app_client
    headers = login_as(client)
    valid_tsv = f"host_name\tgenome_name\tfasta_path\nEscherichia coli\tS1\t{release_data['fasta']}\n".encode("utf-8")

    valid = client.post(
        "/api/reference-database/batch-import/precheck?category=pathogen",
        data={"file": (io.BytesIO(valid_tsv), "references.tsv")},
        headers=headers,
        content_type="multipart/form-data",
    )
    empty = client.post(
        "/api/reference-database/batch-import/precheck?category=pathogen",
        data={"file": (io.BytesIO(b""), "references.tsv")},
        headers=headers,
        content_type="multipart/form-data",
    )

    assert valid.status_code == 200
    assert valid.get_json()["can_import"] is True
    assert empty.status_code == 400
    assert "缺少表头" in empty.get_json()["error"]
