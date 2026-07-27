from __future__ import annotations

from pathlib import Path

def _discover_influenza_igv_assets(report_dir: Path, sample_name: str) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    gff_path = report_dir / f"{sample_name}_virus_typing" / "vadr" / f"{sample_name}.vadr.gff3" if sample_name else Path("")
    gff_asset = ""
    if gff_path.is_file():
        try:
            if gff_path.stat().st_size > 0:
                gff_asset = str(gff_path.relative_to(report_dir))
        except OSError:
            gff_asset = ""
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": gff_asset,
        "viewer_label": "参考比对视图",
        "note": "点击上方变异位点后，IGV 会自动跳转到对应位置。",
    }

def _discover_hiv_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    payload = {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "viewer_label": "HIV 参考比对视图",
        "note": "展示当前 HIV 样本与自动选择代表株参考的比对结果，可结合耐药突变和分型证据一起查看。",
    }
    gff_candidates = [
        report_dir / "ref" / "genes.gff",
        report_dir / "vadr" / f"{report_dir.name}.vadr.gff3",
    ]
    gff_path = next((path for path in gff_candidates if path.is_file() and path.stat().st_size > 0), None)
    if gff_path is not None:
        payload["gff_asset"] = str(gff_path.relative_to(report_dir))
    return payload

def _discover_ncov_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = Path(__file__).resolve().parent.parent / "database" / "virus" / "ncov" / "ref.fna"
    reference_fai = Path(__file__).resolve().parent.parent / "database" / "virus" / "ncov" / "ref.fna.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    ncov_gff = Path(__file__).resolve().parent.parent / "database" / "virus" / "ncov" / "genomic.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and ncov_gff.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": "__ncov_ref.fna",
        "reference_index_asset": "__ncov_ref.fna.fai",
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": "__ncov_genomic.gff",
        "viewer_label": "新冠参考比对视图",
        "note": "点击上方变异位点后，IGV 会自动跳转到对应位置。",
    }

def _discover_monkeypox_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    local_gff = report_dir / "vadr" / f"{report_dir.name}.vadr.gff3"
    hmpxv_db_root = Path(__file__).resolve().parent.parent / "database" / "virus" / "nextclade" / "hMPXV"
    if not hmpxv_db_root.is_dir():
        hmpxv_db_root = Path(__file__).resolve().parent.parent / "database" / "nextclade_db" / "hMPXV"
    hmpxv_gff = hmpxv_db_root / "genome_annotation.gff3"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    if local_gff.is_file() and local_gff.stat().st_size > 0:
        gff_asset = str(local_gff.relative_to(report_dir))
        note = "下方比对视图优先使用本次猴痘 VADR 生成的 GFF3 注释，可联动查看当前样本的参考比对结果。"
    elif hmpxv_gff.is_file():
        gff_asset = "__hmpxv_genome_annotation.gff3"
        note = "下方比对视图使用猴痘参考基因组与注释文件，可联动查看当前样本的参考比对结果。"
    else:
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": gff_asset,
        "viewer_label": "猴痘参考比对视图",
        "note": note,
    }

def _discover_rsv_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "RSV 参考比对视图",
        "note": "展示自动选择的 RSV A/B 参考序列比对结果，可与上方突变位点表联动查看。",
    }

def _discover_hmpv_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "HMPV 参考比对视图",
        "note": "展示 HMPV 固定参考序列比对结果，可与上方突变位点表联动查看。",
    }

def _discover_denv_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "DENV 参考比对视图",
        "note": "展示自动选择的 DENV1-4 参考序列比对结果，可与上方突变位点表联动查看。",
    }

def _discover_zikav_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "ZIKV 参考比对视图",
        "alignment_visibility_window": 3000,
        "alignment_sampling_window_size": 25,
        "alignment_sampling_depth": 25,
        "alignment_height": 180,
        "note": "展示 Zika virus 参考序列比对结果，可与上方突变位点表联动查看。高深度区域已启用 IGV 下采样以避免卡顿。",
    }

def _discover_chikv_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "CHIKV 参考比对视图",
        "alignment_visibility_window": 3000,
        "alignment_sampling_window_size": 25,
        "alignment_sampling_depth": 25,
        "alignment_height": 180,
        "note": "展示 Chikungunya virus 参考序列比对结果，可与上方突变位点表联动查看。高深度区域已启用 IGV 下采样以避免卡顿。",
    }

def _discover_ebola_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "Ebola virus 参考比对视图",
        "alignment_visibility_window": 3000,
        "alignment_sampling_window_size": 25,
        "alignment_sampling_depth": 25,
        "alignment_height": 180,
        "note": "展示自动选择的 Orthoebolavirus 参考序列比对结果，可结合 Nextclade clade / lineage 判读一起查看。",
    }

def _discover_hpiv_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "HPIV 参考比对视图",
        "note": "展示自动选择的 HPIV 最优参考序列比对结果，可与上方突变位点表联动查看。",
    }

def _discover_norovirus_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "Norovirus 参考比对视图",
        "note": "展示自动选择的 Norovirus 最优参考序列比对结果，可与上方突变位点表联动查看。",
    }

def _discover_enterovirus_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_candidates = [
        report_dir / "ref" / "genes.gff",
        report_dir / "vadr" / f"{report_dir.name}.vadr.gff3",
    ]
    gff_path = next((path for path in gff_candidates if path.is_file()), None)
    if gff_path is None:
        vadr_root = report_dir / "vadr"
        if vadr_root.is_dir():
            gff_path = next(iter(sorted(vadr_root.glob("*.vadr.gff3"))), None)
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    payload = {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "viewer_label": "Enterovirus 参考比对视图",
        "note": "展示自动选择的 Enterovirus 最优参考序列比对结果；若存在注释文件会一并加载基因轨道。",
    }
    if gff_path and gff_path.is_file() and gff_path.stat().st_size > 0:
        payload["gff_asset"] = str(gff_path.relative_to(report_dir))
    return payload

def _discover_hepatovirus_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_candidates = [
        report_dir / "ref" / "genes.gff",
        report_dir / "vadr" / f"{report_dir.name}.vadr.gff3",
    ]
    gff_path = next((path for path in gff_candidates if path.is_file()), None)
    if gff_path is None:
        vadr_root = report_dir / "vadr"
        if vadr_root.is_dir():
            gff_path = next(iter(sorted(vadr_root.glob("*.vadr.gff3"))), None)
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    payload = {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "viewer_label": "Hepatovirus 参考比对视图",
        "note": "展示自动选择的肝炎病毒最优参考序列比对结果，可结合大亚型、子亚型/基因型筛选与突变位点一起查看。",
    }
    if gff_path and gff_path.is_file() and gff_path.stat().st_size > 0:
        payload["gff_asset"] = str(gff_path.relative_to(report_dir))
    return payload

def _discover_bandavirus_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_candidates = [
        report_dir / "ref" / "genes.gff",
        report_dir / "vadr" / f"{report_dir.name}.vadr.gff3",
    ]
    gff_path = next((path for path in gff_candidates if path.is_file()), None)
    if gff_path is None:
        vadr_root = report_dir / "vadr"
        if vadr_root.is_dir():
            gff_path = next(iter(sorted(vadr_root.glob("*.vadr.gff3"))), None)
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    payload = {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "viewer_label": "Bandavirus 参考比对视图",
        "note": "展示自动选择的 Bandavirus 最优参考序列比对结果，可结合 A_F/CJ 三片段分型结果一起查看。",
    }
    if gff_path and gff_path.is_file() and gff_path.stat().st_size > 0:
        payload["gff_asset"] = str(gff_path.relative_to(report_dir))
    return payload

def _discover_orthohantavirus_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_candidates = [
        report_dir / "ref" / "genes.gff",
        report_dir / "vadr" / f"{report_dir.name}.vadr.gff3",
    ]
    gff_path = next((path for path in gff_candidates if path.is_file()), None)
    if gff_path is None:
        vadr_root = report_dir / "vadr"
        if vadr_root.is_dir():
            gff_path = next(iter(sorted(vadr_root.glob("*.vadr.gff3"))), None)
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    payload = {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "viewer_label": "Orthohantavirus 参考比对视图",
        "note": "展示自动选择的 Orthohantavirus 最优参考序列比对结果，可结合 broad 筛选和 L/M/S 三片段证据一起查看。",
    }
    if gff_path and gff_path.is_file() and gff_path.stat().st_size > 0:
        payload["gff_asset"] = str(gff_path.relative_to(report_dir))
    return payload

def _discover_astroviridae_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_candidates = [
        report_dir / "ref" / "genes.gff",
        report_dir / "vadr" / f"{report_dir.name}.vadr.gff3",
    ]
    gff_path = next((path for path in gff_candidates if path.is_file()), None)
    if gff_path is None:
        vadr_root = report_dir / "vadr"
        if vadr_root.is_dir():
            gff_path = next(iter(sorted(vadr_root.glob("*.vadr.gff3"))), None)
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    payload = {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "viewer_label": "Astrovirus 参考比对视图",
        "note": "展示自动选择的 Astrovirus 最优参考序列比对结果；若存在注释文件会一并加载基因轨道。",
    }
    if gff_path and gff_path.is_file() and gff_path.stat().st_size > 0:
        payload["gff_asset"] = str(gff_path.relative_to(report_dir))
    return payload

def _discover_rhinovirus_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "Rhinovirus 参考比对视图",
        "note": "展示自动选择的 Rhinovirus 最优参考序列比对结果，可与上方突变位点表联动查看。",
    }

def _discover_seasonal_hcov_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file() and gff_path.is_file()):
        return {"status": "empty"}
    return {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "gff_asset": str(gff_path.relative_to(report_dir)),
        "viewer_label": "季节性冠状病毒参考比对视图",
        "note": "展示自动选择的季节性冠状病毒最优参考序列比对结果，可与上方突变位点和 S 基因系统树联动查看。",
    }

def _discover_hadv_igv_assets(report_dir: Path) -> dict:
    if not report_dir.is_dir():
        return {"status": "empty"}
    reference_fasta = report_dir / "genomes" / "ref.fa"
    reference_fai = report_dir / "genomes" / "ref.fa.fai"
    bam_path = report_dir / "ref.mapping.bam"
    bai_path = report_dir / "ref.mapping.bam.bai"
    gff_path = report_dir / "ref" / "genes.gff"
    if not (reference_fasta.is_file() and reference_fai.is_file() and bam_path.is_file() and bai_path.is_file()):
        return {"status": "empty"}
    payload = {
        "status": "ready",
        "reference_asset": str(reference_fasta.relative_to(report_dir)),
        "reference_index_asset": str(reference_fai.relative_to(report_dir)),
        "bam_asset": str(bam_path.relative_to(report_dir)),
        "bam_index_asset": str(bai_path.relative_to(report_dir)),
        "viewer_label": "HAdV 参考比对视图",
        "note": "展示自动选择的 HAdV 最优参考序列比对结果，可与上方突变位点表联动查看。",
    }
    if gff_path.is_file():
        payload["gff_asset"] = str(gff_path.relative_to(report_dir))
    return payload
