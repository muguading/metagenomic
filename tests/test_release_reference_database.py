from __future__ import annotations

import io

import pytest

from conftest import login_as


pytestmark = pytest.mark.release


def test_release_pathogen_reference_upload_import_accepts_valid_fasta(release_app_client) -> None:
    _app, client = release_app_client
    headers = login_as(client)

    response = client.post(
        "/api/pathogen-database/upload-import",
        data={
            "host_name": "Escherichia coli",
            "genome_name": "S1",
            "taxid": "562",
            "file": (io.BytesIO(b">S1\nACGT\n"), "S1.fasta"),
        },
        headers=headers,
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    payload = response.get_json()
    assert payload["host_name"] == "Escherichia coli"
    assert payload["genome_name"] == "S1"
    assert payload["index_status"] == "pending"


def test_release_reference_upload_rejects_wrong_file_extension(release_app_client) -> None:
    _app, client = release_app_client
    headers = login_as(client)

    response = client.post(
        "/api/pathogen-database/upload-import",
        data={
            "host_name": "Escherichia coli",
            "genome_name": "S1",
            "file": (io.BytesIO(b"not fasta"), "S1.txt"),
        },
        headers=headers,
        content_type="multipart/form-data",
    )

    assert response.status_code == 400
    assert "仅支持 .fa / .fasta / .fna" in response.get_json()["error"]


def test_release_remote_reference_import_rejects_invalid_provider_and_taxid_mode(release_app_client) -> None:
    _app, client = release_app_client
    headers = login_as(client)

    bad_provider = client.post(
        "/api/pathogen-database/remote-import",
        json={"provider": "ensembl", "query": "GCF_000001405.40"},
        headers=headers,
    )
    bad_taxid = client.post(
        "/api/pathogen-database/remote-import",
        json={"provider": "ncbi", "ncbi_mode": "taxid_refs", "query": "not-a-taxid"},
        headers=headers,
    )

    assert bad_provider.status_code == 400
    assert "provider 只支持" in bad_provider.get_json()["error"]
    assert bad_taxid.status_code == 400
    assert "纯数字 TaxID" in bad_taxid.get_json()["error"]


def test_release_reference_upload_rejects_empty_fasta_content(release_app_client) -> None:
    _app, client = release_app_client
    headers = login_as(client)

    response = client.post(
        "/api/pathogen-database/upload-import",
        data={
            "host_name": "Escherichia coli",
            "genome_name": "Empty",
            "file": (io.BytesIO(b""), "empty.fa"),
        },
        headers=headers,
        content_type="multipart/form-data",
    )

    assert response.status_code == 400
