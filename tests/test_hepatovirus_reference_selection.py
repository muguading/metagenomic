from __future__ import annotations

import shlex
from pathlib import Path

from metagenomic_refactor import virus_analysis


def _write_fasta(path: Path, header: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{header}\n{sequence}\n", encoding="utf-8")


def test_resolve_hepatovirus_reference_refines_hav_subtype(tmp_path, monkeypatch) -> None:
    db_root = tmp_path / "database" / "virus" / "Hepatovirus"
    broad_dir = db_root / "reference_genomes"
    subtype_dir = db_root / "HAV_subtypes"

    _write_fasta(broad_dir / "M14707.fasta", "M14707.1 HM175", "ATGC" * 100)
    _write_fasta(broad_dir / "KT452742.fasta", "KT452742.1 HepV-C", "ATGC" * 100)
    (broad_dir / "gff3").mkdir(parents=True, exist_ok=True)
    (broad_dir / "gff3" / "M14707.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (broad_dir / "gff3" / "KT452742.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (broad_dir / "hepatovirus_typing_reference_genomes_manifest.tsv").write_text(
        "\t".join(
            [
                "genus",
                "species",
                "virus_name",
                "isolate",
                "accession",
                "available_sequence",
                "abbrev",
                "header",
                "sequence_length",
                "fasta_path",
                "gff3_path",
                "fasta_status",
                "gff3_status",
                "status",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "Hepatovirus",
                "Hepatovirus ahepa",
                "hepatovirus A1; hepatitis A virus",
                "HM175",
                "M14707",
                "Complete genome",
                "HAV",
                "M14707.1 HM175",
                "400",
                str((broad_dir / "M14707.fasta").resolve()),
                str((broad_dir / "gff3" / "M14707.gff3").resolve()),
                "downloaded",
                "downloaded",
                "ok",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "Hepatovirus",
                "Hepatovirus cemanavi",
                "hepatovirus C1; bat hepatovirus",
                "SMG18528Minmav2014",
                "KT452742",
                "Complete genome",
                "HepV-C",
                "KT452742.1 HepV-C",
                "400",
                str((broad_dir / "KT452742.fasta").resolve()),
                str((broad_dir / "gff3" / "KT452742.gff3").resolve()),
                "downloaded",
                "downloaded",
                "ok",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    _write_fasta(subtype_dir / "fasta" / "IA__K02990__LA.fasta", "K02990.1 LA", "ATGC" * 100)
    _write_fasta(subtype_dir / "fasta" / "IB__M14707__HM175.fasta", "M14707.1 HM175", "ATGC" * 100)
    (subtype_dir / "gff3").mkdir(parents=True, exist_ok=True)
    (subtype_dir / "gff3" / "IA__K02990__LA.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (subtype_dir / "gff3" / "IB__M14707__HM175.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (subtype_dir / "hav_subtype_complete_genomes_manifest.tsv").write_text(
        "\t".join(
            [
                "genotype",
                "accession",
                "isolate",
                "header",
                "sequence_length",
                "fasta_path",
                "gff3_path",
                "fasta_status",
                "gff3_status",
                "status",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "IA",
                "K02990",
                "LA",
                "K02990.1 LA",
                "400",
                str((subtype_dir / "fasta" / "IA__K02990__LA.fasta").resolve()),
                str((subtype_dir / "gff3" / "IA__K02990__LA.gff3").resolve()),
                "downloaded",
                "downloaded",
                "ok",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "IB",
                "M14707",
                "HM175",
                "M14707.1 HM175",
                "400",
                str((subtype_dir / "fasta" / "IB__M14707__HM175.fasta").resolve()),
                str((subtype_dir / "gff3" / "IB__M14707__HM175.gff3").resolve()),
                "downloaded",
                "downloaded",
                "ok",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    def fake_multi_reference_coverage(reference_fasta: Path, output_prefix: Path, **kwargs):
        del reference_fasta, kwargs
        if "broad_typing" in output_prefix.name:
            return [
                {"reference_name": "M14707.1", "coverage": 98.0, "mean_depth": 50.0, "covered_bases": 395.0, "num_reads": 120.0},
                {"reference_name": "KT452742.1", "coverage": 61.0, "mean_depth": 18.0, "covered_bases": 244.0, "num_reads": 70.0},
            ]
        return [
            {"reference_name": "K02990.1", "coverage": 97.0, "mean_depth": 42.0, "covered_bases": 388.0, "num_reads": 100.0},
            {"reference_name": "M14707.1", "coverage": 92.0, "mean_depth": 35.0, "covered_bases": 368.0, "num_reads": 88.0},
        ]

    monkeypatch.setenv("META_HEPATOVIRUS_DB_DIR", str(db_root))
    monkeypatch.setattr(virus_analysis, "_run_multi_reference_coverage", fake_multi_reference_coverage)
    monkeypatch.chdir(tmp_path)

    result = virus_analysis.resolve_hepatovirus_reference(
        "sample1",
        species="Hepatovirus",
        single_fastq="dummy.fastq",
    )

    assert result["status"] == "ready"
    assert result["broad_type"] == "HAV"
    assert result["hav_subtype"] == "IA"
    assert result["species_label"] == "Hepatitis A virus"
    assert result["reference_path"].endswith("IA.reference.fasta")
    assert result["gff_path"].endswith("IA__K02990__LA.gff3")


def test_run_hepatovirus_consensus_typing_keeps_selection_outputs_separate(tmp_path, monkeypatch) -> None:
    selection_dir = tmp_path / "sample1_hepatovirus_reference_selection"
    selection_dir.mkdir(parents=True, exist_ok=True)
    (selection_dir / "selection.tsv").write_text("sample\tstatus\nsample1\treads\n", encoding="utf-8")
    consensus_fasta = tmp_path / "sample1.consensus.fasta"
    consensus_fasta.write_text(">sample1\nATGC\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_resolve_hepatovirus_reference(pre: str, **kwargs):
        captured["pre"] = pre
        captured["output_dir"] = kwargs.get("output_dir")
        output_dir = Path(kwargs["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        selection_path = output_dir / "selection.tsv"
        selection_path.write_text("sample\tstatus\nsample1\tconsensus\n", encoding="utf-8")
        return {
            "status": "ready",
            "broad_type": "HAV",
            "hav_subtype": "IA",
            "summary_path": str(selection_path.resolve()),
        }

    monkeypatch.setattr(virus_analysis, "resolve_hepatovirus_reference", fake_resolve_hepatovirus_reference)
    monkeypatch.chdir(tmp_path)

    result = virus_analysis.run_hepatovirus_consensus_typing("sample1", consensus_fasta)

    assert result["status"] == "ready"
    assert result["broad_type"] == "HAV"
    assert result["hav_subtype"] == "IA"
    assert Path(captured["output_dir"]) == Path("sample1_hepatovirus_reference_selection") / "consensus_typing"
    assert (selection_dir / "selection.tsv").read_text(encoding="utf-8") == "sample\tstatus\nsample1\treads\n"
    assert (selection_dir / "consensus_typing" / "selection.tsv").is_file()
    assert (selection_dir / "consensus_typing" / "consensus_typing.tsv").is_file()


def test_resolve_hepatovirus_reference_refines_hev_subtype(tmp_path, monkeypatch) -> None:
    db_root = tmp_path / "database" / "virus" / "Hepatovirus"
    broad_dir = db_root / "broad_reference_genomes"
    subtype_dir = db_root / "typingE_reference_genomes"

    _write_fasta(broad_dir / "M73218.fasta", "M73218.1 HEV1a", "ATGC" * 100)
    _write_fasta(broad_dir / "AF082843.fasta", "AF082843.1 HEV3a", "ATGC" * 100)
    (broad_dir / "gff3").mkdir(parents=True, exist_ok=True)
    (broad_dir / "gff3" / "M73218.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (broad_dir / "gff3" / "AF082843.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (broad_dir / "hepatitis_broad_reference_genomes_manifest.tsv").write_text(
        "typing_source\tbroad_type\tspecies_label\tgenus\tspecies\tvirus_name\tisolate\taccession\tavailable_sequence\tabbrev\theader\tsequence_length\tfasta_path\tgff3_path\tfasta_status\tgff3_status\tstatus\n"
        f"E\tHEV\tHepatitis E virus\tPaslahepevirus\tPaslahepevirus balayani\thuman hepatitis E virus genotype 1a\tBurma\tM73218\tComplete genome\tHEV\tM73218.1 HEV1a\t400\t{(broad_dir / 'M73218.fasta').resolve()}\t{(broad_dir / 'gff3' / 'M73218.gff3').resolve()}\tdownloaded\tdownloaded\tok\n"
        f"E\tHEV\tHepatitis E virus\tPaslahepevirus\tPaslahepevirus balayani\tswine hepatitis E virus genotype 3a\tMeng\tAF082843\tComplete genome\tHEV\tAF082843.1 HEV3a\t400\t{(broad_dir / 'AF082843.fasta').resolve()}\t{(broad_dir / 'gff3' / 'AF082843.gff3').resolve()}\tdownloaded\tdownloaded\tok\n",
        encoding="utf-8",
    )

    _write_fasta(subtype_dir / "M73218.fasta", "M73218.1 HEV1a", "ATGC" * 100)
    _write_fasta(subtype_dir / "AF082843.fasta", "AF082843.1 HEV3a", "ATGC" * 100)
    (subtype_dir / "gff3").mkdir(parents=True, exist_ok=True)
    (subtype_dir / "gff3" / "M73218.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (subtype_dir / "gff3" / "AF082843.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (subtype_dir / "typingE_reference_genomes_manifest.tsv").write_text(
        "genus\tspecies\tvirus_name\tisolate\taccession\tavailable_sequence\tabbrev\theader\tsequence_length\tfasta_path\tgff3_path\tfasta_status\tgff3_status\tstatus\n"
        f"Paslahepevirus\tPaslahepevirus balayani\thuman hepatitis E virus genotype 1a\tBurma\tM73218\tComplete genome\tHEV\tM73218.1 HEV1a\t400\t{(subtype_dir / 'M73218.fasta').resolve()}\t{(subtype_dir / 'gff3' / 'M73218.gff3').resolve()}\tdownloaded\tdownloaded\tok\n"
        f"Paslahepevirus\tPaslahepevirus balayani\tswine hepatitis E virus genotype 3a\tMeng\tAF082843\tComplete genome\tHEV\tAF082843.1 HEV3a\t400\t{(subtype_dir / 'AF082843.fasta').resolve()}\t{(subtype_dir / 'gff3' / 'AF082843.gff3').resolve()}\tdownloaded\tdownloaded\tok\n",
        encoding="utf-8",
    )

    def fake_multi_reference_coverage(reference_fasta: Path, output_prefix: Path, **kwargs):
        del reference_fasta, kwargs
        if "broad_typing" in output_prefix.name:
            return [
                {"reference_name": "M73218.1", "coverage": 91.0, "mean_depth": 25.0, "covered_bases": 360.0, "num_reads": 90.0},
                {"reference_name": "AF082843.1", "coverage": 96.0, "mean_depth": 40.0, "covered_bases": 384.0, "num_reads": 120.0},
            ]
        return [
            {"reference_name": "M73218.1", "coverage": 84.0, "mean_depth": 20.0, "covered_bases": 336.0, "num_reads": 70.0},
            {"reference_name": "AF082843.1", "coverage": 99.0, "mean_depth": 52.0, "covered_bases": 396.0, "num_reads": 150.0},
        ]

    monkeypatch.setenv("META_HEPATOVIRUS_DB_DIR", str(db_root))
    monkeypatch.setattr(virus_analysis, "_run_multi_reference_coverage", fake_multi_reference_coverage)
    monkeypatch.chdir(tmp_path)

    result = virus_analysis.resolve_hepatovirus_reference(
        "sample_hev",
        species="Hepatitis E virus",
        single_fastq="dummy.fastq",
    )

    assert result["status"] == "ready"
    assert result["broad_type"] == "HEV"
    assert result["subtype"] == "3a"
    assert result["hav_subtype"] == ""
    assert result["species_label"] == "Hepatitis E virus"
    assert result["reference_path"].endswith("3a.reference.fasta")
    assert result["gff_path"].endswith("AF082843.gff3")
    summary_text = Path(result["summary_path"]).read_text(encoding="utf-8")
    assert "\tsubtype\t" in summary_text.splitlines()[0]
    assert "\tHEV\t3a\t" in summary_text


def test_virus_typing_writes_hev_subtype_result_not_influenza_placeholder(tmp_path, monkeypatch) -> None:
    sample = "sample_hev"
    selection_dir = tmp_path / f"{sample}_hepatovirus_reference_selection"
    selection_dir.mkdir(parents=True)
    selection_dir.joinpath("selection.tsv").write_text(
        "sample\tspecies_label\tbroad_type\tsubtype\thav_subtype\tisolate\tcoverage\tmean_depth\tcovered_bases\tnum_reads\treference_name\taccession\treference_path\tgff_path\tstatus\tnote\tbroad_summary_path\tsubtype_summary_path\n"
        f"{sample}\tHepatitis E virus\tHEV\t3a\t\tMeng\t99.0\t42.0\t7200\t120\tAF082843.1\tAF082843\t/tmp/3a.reference.fasta\t/tmp/AF082843.gff3\tready\tHEV subtype selected\t/tmp/broad.tsv\t/tmp/hev.tsv\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    virus_analysis.virus_typing(sample, "Hepatitis E virus")

    result_text = (tmp_path / f"{sample}_serotype_result.tsv").read_text(encoding="utf-8")
    assert "HA亚型" not in result_text
    assert "子亚型" in result_text
    assert "Hepatitis E virus\tHEV\t3A" in result_text


def test_multi_reference_blast_merges_split_hsps_for_reference_coverage(tmp_path, monkeypatch) -> None:
    query_fasta = tmp_path / "query.fasta"
    reference_fasta = tmp_path / "reference.fasta"
    query_fasta.write_text(">query\n" + ("A" * 180) + "\n", encoding="utf-8")
    reference_fasta.write_text(">refA\n" + ("A" * 200) + "\n>refB\n" + ("A" * 200) + "\n", encoding="utf-8")
    output_prefix = tmp_path / "hepatitis_subtype_typing"
    meta_by_id = {
        "refA": {"reference_length": 200},
        "refB": {"reference_length": 200},
    }

    def fake_run_command(cmd: str, logf=None) -> None:
        del logf
        parts = shlex.split(cmd)
        out_path = Path(parts[parts.index("-out") + 1])
        out_path.write_text(
            "\n".join(
                [
                    "query\trefA\t99.0\t80\t0\t0\t1\t80\t1\t80\t1e-80\t100",
                    "query\trefA\t98.0\t70\t0\t0\t91\t160\t91\t160\t1e-60\t90",
                    "query\trefA\t97.0\t40\t0\t0\t111\t150\t111\t150\t1e-30\t50",
                    "query\trefB\t100.0\t120\t0\t0\t1\t120\t1\t120\t1e-90\t150",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(virus_analysis, "run_command", fake_run_command)

    rows = virus_analysis._run_multi_reference_blast(query_fasta, reference_fasta, output_prefix, meta_by_id)

    assert rows[0]["reference_name"] == "refA"
    assert float(rows[0]["covered_bases"]) == 150.0
    assert float(rows[0]["coverage"]) == 75.0
    assert float(rows[1]["covered_bases"]) == 120.0
    coverage_text = output_prefix.with_suffix(".coverage.tsv").read_text(encoding="utf-8")
    assert "refA\t98.210526\t75.000000" in coverage_text


def test_resolve_hepatovirus_reference_preserves_versioned_fasta_id_for_snpeff(tmp_path, monkeypatch) -> None:
    db_root = tmp_path / "database" / "virus" / "Hepatovirus"
    broad_dir = db_root / "reference_genomes"
    subtype_dir = db_root / "HAV_subtypes"

    _write_fasta(broad_dir / "M14707.fasta", "M14707.1 HM175", "ATGC" * 50)
    (broad_dir / "gff3").mkdir(parents=True, exist_ok=True)
    (broad_dir / "gff3" / "M14707.gff3").write_text("##gff-version 3\nM14707.1\tGenbank\tregion\t1\t200\t.\t+\t.\tID=x\n", encoding="utf-8")
    (broad_dir / "hepatovirus_typing_reference_genomes_manifest.tsv").write_text(
        "genus\tspecies\tvirus_name\tisolate\taccession\tavailable_sequence\tabbrev\theader\tsequence_length\tfasta_path\tgff3_path\tfasta_status\tgff3_status\tstatus\n"
        f"Hepatovirus\tHepatovirus ahepa\thepatitis A virus\tHM175\tM14707\tComplete genome\tHAV\tM14707.1 HM175\t200\t{(broad_dir / 'M14707.fasta').resolve()}\t{(broad_dir / 'gff3' / 'M14707.gff3').resolve()}\tdownloaded\tdownloaded\tok\n",
        encoding="utf-8",
    )

    _write_fasta(subtype_dir / "fasta" / "IB__M14707__HM175.fasta", "M14707.1 HM175", "ATGC" * 50)
    (subtype_dir / "gff3").mkdir(parents=True, exist_ok=True)
    (subtype_dir / "gff3" / "IB__M14707__HM175.gff3").write_text("##gff-version 3\nM14707.1\tGenbank\tregion\t1\t200\t.\t+\t.\tID=x\n", encoding="utf-8")
    (subtype_dir / "hav_subtype_complete_genomes_manifest.tsv").write_text(
        "genotype\taccession\tisolate\theader\tsequence_length\tfasta_path\tgff3_path\tfasta_status\tgff3_status\tstatus\n"
        f"IB\tM14707\tHM175\tM14707.1 HM175\t200\t{(subtype_dir / 'fasta' / 'IB__M14707__HM175.fasta').resolve()}\t{(subtype_dir / 'gff3' / 'IB__M14707__HM175.gff3').resolve()}\tdownloaded\tdownloaded\tok\n",
        encoding="utf-8",
    )

    def fake_multi_reference_coverage(reference_fasta: Path, output_prefix: Path, **kwargs):
        del reference_fasta, kwargs
        return [{"reference_name": "M14707.1", "coverage": 99.0, "mean_depth": 12.0, "covered_bases": 198.0, "num_reads": 33.0}]

    monkeypatch.setenv("META_HEPATOVIRUS_DB_DIR", str(db_root))
    monkeypatch.setattr(virus_analysis, "_run_multi_reference_coverage", fake_multi_reference_coverage)
    monkeypatch.chdir(tmp_path)

    result = virus_analysis.resolve_hepatovirus_reference("sample2", species="Hepatovirus", single_fastq="dummy.fastq")

    reference_path = Path(result["reference_path"])
    assert reference_path.is_file()
    assert reference_path.read_text(encoding="utf-8").startswith(">M14707.1")
