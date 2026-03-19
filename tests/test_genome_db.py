from __future__ import annotations

import time
from pathlib import Path

import pytest
from sqlalchemy import select

from genome_db.database import session_scope
from genome_db.genome_manager import DuplicateGenomeError, GenomeManager
from genome_db.models import AuditLog, Genome
from genome_db.validators import ValidationError


@pytest.fixture
def manager(tmp_path: Path) -> GenomeManager:
    db_path = tmp_path / "genomes.sqlite3"
    return GenomeManager(database_url=f"sqlite:///{db_path}")


def write_fasta(path: Path, header: str, sequence: str) -> Path:
    path.write_text(f">{header}\n{sequence}\n", encoding="utf-8")
    return path


def fetch_audit_rows(manager: GenomeManager) -> list[AuditLog]:
    with session_scope(manager.session_factory) as session:
        return list(session.scalars(select(AuditLog).order_by(AuditLog.id)).all())


def test_add_and_get_genome_write_success_audit_logs(manager: GenomeManager, tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "g1.fasta", "chr1", "ATGCATGC")

    added = manager.add_genome(
        genome_id="G001",
        sample_name="sample_001",
        species_name="Escherichia coli",
        taxid=562,
        genome_file_path=str(fasta),
        submitter="alice",
        description="reference genome",
    )

    fetched = manager.get_genome("G001", operator="alice")

    assert added["genome_id"] == "G001"
    assert added["genome_length"] == 8
    assert added["submit_time"] is not None
    assert added["last_modified_time"] is not None
    assert fetched["species_name"] == "Escherichia coli"
    assert fetched["genome_file_path"] == str(fasta.resolve())

    audit_rows = fetch_audit_rows(manager)
    assert [(row.operation, row.status) for row in audit_rows] == [
        ("add_genome", "SUCCESS"),
        ("get_genome", "SUCCESS"),
    ]


def test_duplicate_genome_is_rejected_and_failure_is_audited(manager: GenomeManager, tmp_path: Path) -> None:
    fasta1 = write_fasta(tmp_path / "g1.fasta", "chr1", "ATGCATGC")
    fasta2 = write_fasta(tmp_path / "g2.fasta", "chr1", "ATGC")

    manager.add_genome(
        genome_id="G001",
        sample_name="sample_001",
        species_name="Escherichia coli",
        taxid=562,
        genome_file_path=str(fasta1),
        submitter="alice",
    )

    with pytest.raises(DuplicateGenomeError):
        manager.add_genome(
            genome_id="G001",
            sample_name="sample_002",
            species_name="Shigella sonnei",
            taxid=624,
            genome_file_path=str(fasta2),
            submitter="bob",
        )

    audit_rows = fetch_audit_rows(manager)
    assert [(row.operation, row.status) for row in audit_rows] == [
        ("add_genome", "SUCCESS"),
        ("add_genome", "FAILED"),
    ]


def test_invalid_fasta_is_rejected_and_failure_is_audited(manager: GenomeManager, tmp_path: Path) -> None:
    invalid_fasta = tmp_path / "bad.fasta"
    invalid_fasta.write_text("chr1\nATGC\n", encoding="utf-8")

    with pytest.raises(ValidationError):
        manager.add_genome(
            genome_id="BAD001",
            sample_name="sample_bad",
            species_name="Bacillus subtilis",
            taxid=1423,
            genome_file_path=str(invalid_fasta),
            submitter="alice",
        )

    audit_rows = fetch_audit_rows(manager)
    assert len(audit_rows) == 1
    assert audit_rows[0].operation == "add_genome"
    assert audit_rows[0].status == "FAILED"


def test_update_genome_refreshes_data_and_audits(manager: GenomeManager, tmp_path: Path) -> None:
    fasta1 = write_fasta(tmp_path / "g1.fasta", "chr1", "ATGCATGC")
    fasta2 = write_fasta(tmp_path / "g2.fasta", "chr1", "ATGCATGCAT")

    manager.add_genome(
        genome_id="G001",
        sample_name="sample_001",
        species_name="Escherichia coli",
        taxid=562,
        genome_file_path=str(fasta1),
        submitter="alice",
    )
    original = manager.get_genome("G001", operator="alice")

    time.sleep(0.01)
    updated = manager.update_genome(
        "G001",
        operator="bob",
        species_name="Shigella flexneri",
        taxid=623,
        genome_file_path=str(fasta2),
        description="updated record",
    )

    assert updated["species_name"] == "Shigella flexneri"
    assert updated["taxid"] == 623
    assert updated["genome_length"] == 10
    assert updated["description"] == "updated record"
    assert updated["last_modified_time"] >= original["last_modified_time"]

    audit_rows = fetch_audit_rows(manager)
    assert [(row.operation, row.status) for row in audit_rows] == [
        ("add_genome", "SUCCESS"),
        ("get_genome", "SUCCESS"),
        ("update_genome", "SUCCESS"),
    ]


def test_delete_genome_removes_record_and_audits(manager: GenomeManager, tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "g1.fasta", "chr1", "ATGCATGC")
    manager.add_genome(
        genome_id="G001",
        sample_name="sample_001",
        species_name="Escherichia coli",
        taxid=562,
        genome_file_path=str(fasta),
        submitter="alice",
    )

    manager.delete_genome("G001", operator="alice")

    with pytest.raises(KeyError):
        manager.get_genome("G001", operator="alice")

    with session_scope(manager.session_factory) as session:
        assert session.scalar(select(Genome).where(Genome.genome_id == "G001")) is None

    audit_rows = fetch_audit_rows(manager)
    assert [(row.operation, row.status) for row in audit_rows] == [
        ("add_genome", "SUCCESS"),
        ("delete_genome", "SUCCESS"),
        ("get_genome", "FAILED"),
    ]


def test_search_supports_filters_pagination_and_audit(manager: GenomeManager, tmp_path: Path) -> None:
    records = [
        ("G001", "sample_001", "Klebsiella pneumoniae", 573, "alice", "ATGCATGC"),
        ("G002", "sample_002", "Klebsiella variicola", 244366, "alice", "ATGCAT"),
        ("G003", "sample_003", "Escherichia coli", 562, "bob", "ATGCATGCAT"),
    ]
    for genome_id, sample_name, species_name, taxid, submitter, sequence in records:
        fasta = write_fasta(tmp_path / f"{genome_id}.fasta", genome_id, sequence)
        manager.add_genome(
            genome_id=genome_id,
            sample_name=sample_name,
            species_name=species_name,
            taxid=taxid,
            genome_file_path=str(fasta),
            submitter=submitter,
        )

    result = manager.search_genomes(species_name="Klebsiella", page=1, page_size=1, operator="alice")
    page2 = manager.search_genomes(species_name="Klebsiella", page=2, page_size=1, operator="alice")
    by_taxid = manager.search_genomes(taxid=562, page=1, page_size=10, operator="bob")
    by_submitter = manager.search_genomes(submitter="alice", page=1, page_size=10, operator="alice")

    assert result.total == 2
    assert result.page == 1
    assert result.page_size == 1
    assert result.pages == 2
    assert len(result.items) == 1
    assert len(page2.items) == 1
    assert result.items[0]["genome_id"] != page2.items[0]["genome_id"]
    assert by_taxid.total == 1
    assert by_taxid.items[0]["species_name"] == "Escherichia coli"
    assert by_submitter.total == 2

    audit_rows = fetch_audit_rows(manager)
    search_rows = [row for row in audit_rows if row.operation == "search_genomes"]
    assert len(search_rows) == 4
    assert all(row.status == "SUCCESS" for row in search_rows)


def test_failed_search_validation_is_audited(manager: GenomeManager) -> None:
    with pytest.raises(ValidationError):
        manager.search_genomes(page=0, page_size=10, operator="alice")

    audit_rows = fetch_audit_rows(manager)
    assert len(audit_rows) == 1
    assert audit_rows[0].operation == "search_genomes"
    assert audit_rows[0].status == "FAILED"
