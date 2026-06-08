#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_mag_binning_pipeline.sh --manifest samples.tsv --outdir mag_binning_out [options]

Required:
  --manifest FILE            Tab-separated sample sheet with header:
                             sample<TAB>contigs<TAB>bam_files
                             bam_files must be a comma-separated list of sorted BAM files.
  --outdir DIR               Output directory.

Optional:
  --threads INT              Threads per sample run. Default: 16
  --min-contig-len INT       Minimum contig length for MetaBAT2/VAMB. Default: 1500
  --semibin-env NAME         SemiBin2 pretrained environment. Default: global
  --semibin-seq-type TYPE    SemiBin2 sequencing type: short_read or long_read. Default: short_read
  --vamb-minfasta INT        Minimum total bin size for VAMB FASTA export. Default: 200000
  --score-threshold FLOAT    DASTool score threshold. Default: 0.1
  --force                    Re-run finished steps if outputs already exist.
  -h, --help                 Show this help message.

Example manifest:
sample  contigs bam_files
S1      /path/S1.contigs.fa /path/S1.sorted.bam
S2      /path/S2.contigs.fa /path/S2_A.sorted.bam,/path/S2_B.sorted.bam
EOF
}

log() {
    local ts
    ts="$(date '+%F %T')"
    printf '[%s] %s\n' "$ts" "$*" >&2
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
}

realpath_compat() {
    local target="$1"
    if command -v realpath >/dev/null 2>&1; then
        realpath "$target"
    else
        perl -MCwd=abs_path -e 'print abs_path(shift), "\n"' "$target"
    fi
}

cat_any() {
    local file="$1"
    case "$file" in
        *.gz) gzip -dc "$file" ;;
        *) cat "$file" ;;
    esac
}

bins_dir_to_contigs2bin() {
    local bins_dir="$1"
    local out_tsv="$2"
    local found=0
    : >"$out_tsv"

    shopt -s nullglob
    local fasta
    for fasta in "$bins_dir"/*.fa "$bins_dir"/*.fna "$bins_dir"/*.fasta "$bins_dir"/*.fa.gz "$bins_dir"/*.fna.gz "$bins_dir"/*.fasta.gz; do
        found=1
        local bin_name
        bin_name="$(basename "$fasta")"
        bin_name="${bin_name%.gz}"
        bin_name="${bin_name%.fa}"
        bin_name="${bin_name%.fna}"
        bin_name="${bin_name%.fasta}"

        cat_any "$fasta" | awk -v bin="$bin_name" '
            /^>/ {
                header = substr($0, 2)
                split(header, parts, /[ \t]/)
                print parts[1] "\t" bin
            }
        ' >>"$out_tsv"
    done
    shopt -u nullglob

    [[ "$found" -eq 1 ]] || die "No FASTA bins found in: $bins_dir"
    [[ -s "$out_tsv" ]] || die "Generated empty contigs2bin table: $out_tsv"
}

vamb_clusters_to_contigs2bin() {
    local clusters_tsv="$1"
    local out_tsv="$2"
    awk -F'\t' 'NR == 1 { next } NF >= 2 { print $2 "\t" $1 }' "$clusters_tsv" >"$out_tsv"
    [[ -s "$out_tsv" ]] || die "Generated empty VAMB contigs2bin table: $out_tsv"
}

link_bams_into_dir() {
    local bam_csv="$1"
    local bam_dir="$2"
    mkdir -p "$bam_dir"

    local IFS=','
    local bam
    local idx=0
    for bam in $bam_csv; do
        bam="${bam#"${bam%%[![:space:]]*}"}"
        bam="${bam%"${bam##*[![:space:]]}"}"
        [[ -n "$bam" ]] || continue
        [[ -f "$bam" ]] || die "BAM file not found: $bam"
        idx=$((idx + 1))
        ln -sfn "$(realpath_compat "$bam")" "$bam_dir/$(printf '%03d' "$idx")_$(basename "$bam")"
    done
    [[ "$idx" -gt 0 ]] || die "No BAM files parsed from: $bam_csv"
}

run_metabat2() {
    local sample="$1"
    local contigs="$2"
    local bam_csv="$3"
    local sample_out="$4"

    local metabat_dir="$sample_out/metabat2"
    local depth_tsv="$metabat_dir/${sample}.depth.tsv"
    local bin_prefix="$metabat_dir/${sample}.bin"
    local done_flag="$metabat_dir/.done"
    mkdir -p "$metabat_dir"

    if [[ -f "$done_flag" && "$FORCE" -eq 0 ]]; then
        log "[$sample] Skip MetaBAT2: existing results found."
        return
    fi

    local IFS=','
    local bam_array=()
    local bam
    for bam in $bam_csv; do
        bam="${bam#"${bam%%[![:space:]]*}"}"
        bam="${bam%"${bam##*[![:space:]]}"}"
        [[ -n "$bam" ]] || continue
        bam_array+=("$bam")
    done
    [[ "${#bam_array[@]}" -gt 0 ]] || die "[$sample] No BAM files available for MetaBAT2"

    log "[$sample] Running jgi_summarize_bam_contig_depths"
    jgi_summarize_bam_contig_depths --outputDepth "$depth_tsv" "${bam_array[@]}" \
        >"$metabat_dir/jgi.stdout.log" 2>"$metabat_dir/jgi.stderr.log"

    log "[$sample] Running MetaBAT2"
    metabat2 \
        -i "$contigs" \
        -a "$depth_tsv" \
        -o "$bin_prefix" \
        -m "$MIN_CONTIG_LEN" \
        -t "$THREADS" \
        >"$metabat_dir/metabat2.stdout.log" 2>"$metabat_dir/metabat2.stderr.log"

    touch "$done_flag"
}

run_semibin2() {
    local sample="$1"
    local contigs="$2"
    local bam_csv="$3"
    local sample_out="$4"

    local semibin_dir="$sample_out/semibin2"
    local done_flag="$semibin_dir/.done"
    mkdir -p "$semibin_dir"

    if [[ -f "$done_flag" && "$FORCE" -eq 0 ]]; then
        log "[$sample] Skip SemiBin2: existing results found."
        return
    fi

    local IFS=','
    local bam_array=()
    local bam
    for bam in $bam_csv; do
        bam="${bam#"${bam%%[![:space:]]*}"}"
        bam="${bam%"${bam##*[![:space:]]}"}"
        [[ -n "$bam" ]] || continue
        bam_array+=("$bam")
    done
    [[ "${#bam_array[@]}" -gt 0 ]] || die "[$sample] No BAM files available for SemiBin2"

    log "[$sample] Running SemiBin2"
    SemiBin2 single_easy_bin \
        --input-fasta "$contigs" \
        --input-bam "${bam_array[@]}" \
        --output "$semibin_dir" \
        --environment "$SEMIBIN_ENV" \
        --threads "$THREADS" \
        --sequencing-type "$SEMIBIN_SEQ_TYPE" \
        >"$semibin_dir/semibin2.stdout.log" 2>"$semibin_dir/semibin2.stderr.log"

    touch "$done_flag"
}

run_vamb() {
    local sample="$1"
    local contigs="$2"
    local bam_csv="$3"
    local sample_out="$4"

    local vamb_parent="$sample_out"
    local vamb_dir="$vamb_parent/vamb"
    local bam_dir="$vamb_parent/vamb_bams"
    local done_flag="$vamb_parent/vamb.done"

    if [[ -f "$done_flag" && "$FORCE" -eq 0 ]]; then
        log "[$sample] Skip VAMB: existing results found."
        return
    fi

    rm -rf "$bam_dir"
    mkdir -p "$vamb_parent"
    link_bams_into_dir "$bam_csv" "$bam_dir"

    if [[ -d "$vamb_dir" ]]; then
        if [[ "$FORCE" -eq 1 ]]; then
            rm -rf "$vamb_dir"
        else
            die "[$sample] VAMB output directory already exists: $vamb_dir. Use --force to overwrite."
        fi
    fi

    log "[$sample] Running VAMB"
    vamb bin default \
        --outdir "$vamb_dir" \
        --fasta "$contigs" \
        --bamdir "$bam_dir" \
        --minfasta "$VAMB_MIN_FASTA" \
        -m "$MIN_CONTIG_LEN" \
        -p "$THREADS" \
        -o \
        >"$vamb_parent/vamb.stdout.log" 2>"$vamb_parent/vamb.stderr.log"

    touch "$done_flag"
}

run_dastool() {
    local sample="$1"
    local contigs="$2"
    local sample_out="$3"

    local dastool_dir="$sample_out/dastool"
    local done_flag="$dastool_dir/.done"
    mkdir -p "$dastool_dir"

    if [[ -f "$done_flag" && "$FORCE" -eq 0 ]]; then
        log "[$sample] Skip DASTool: existing results found."
        return
    fi

    local metabat_bins="$sample_out/metabat2"
    local semibin_bins="$sample_out/semibin2/output_bins"
    local vamb_clusters="$sample_out/vamb/vae_clusters_unsplit.tsv"

    [[ -d "$metabat_bins" ]] || die "[$sample] MetaBAT2 bin directory not found: $metabat_bins"
    [[ -d "$semibin_bins" ]] || die "[$sample] SemiBin2 output_bins not found: $semibin_bins"

    if [[ -f "$sample_out/vamb/vae_clusters_split.tsv" ]]; then
        vamb_clusters="$sample_out/vamb/vae_clusters_split.tsv"
    fi
    [[ -f "$vamb_clusters" ]] || die "[$sample] VAMB cluster file not found: $vamb_clusters"

    local metabat_tsv="$dastool_dir/metabat2_contigs2bin.tsv"
    local semibin_tsv="$dastool_dir/semibin2_contigs2bin.tsv"
    local vamb_tsv="$dastool_dir/vamb_contigs2bin.tsv"
    local out_prefix="$dastool_dir/${sample}"

    log "[$sample] Preparing contigs2bin tables for DASTool"
    bins_dir_to_contigs2bin "$metabat_bins" "$metabat_tsv"
    bins_dir_to_contigs2bin "$semibin_bins" "$semibin_tsv"
    vamb_clusters_to_contigs2bin "$vamb_clusters" "$vamb_tsv"

    log "[$sample] Running DASTool"
    DAS_Tool \
        -i "$metabat_tsv,$semibin_tsv,$vamb_tsv" \
        -l "metabat2,semibin2,vamb" \
        -c "$contigs" \
        -o "$out_prefix" \
        --write_bins \
        --score_threshold "$SCORE_THRESHOLD" \
        --threads "$THREADS" \
        >"$dastool_dir/dastool.stdout.log" 2>"$dastool_dir/dastool.stderr.log"

    touch "$done_flag"
}

MANIFEST=""
OUTDIR=""
THREADS=16
MIN_CONTIG_LEN=1500
SEMIBIN_ENV="global"
SEMIBIN_SEQ_TYPE="short_read"
VAMB_MIN_FASTA=200000
SCORE_THRESHOLD=0.1
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)
            MANIFEST="${2:-}"
            shift 2
            ;;
        --outdir)
            OUTDIR="${2:-}"
            shift 2
            ;;
        --threads)
            THREADS="${2:-}"
            shift 2
            ;;
        --min-contig-len)
            MIN_CONTIG_LEN="${2:-}"
            shift 2
            ;;
        --semibin-env)
            SEMIBIN_ENV="${2:-}"
            shift 2
            ;;
        --semibin-seq-type)
            SEMIBIN_SEQ_TYPE="${2:-}"
            shift 2
            ;;
        --vamb-minfasta)
            VAMB_MIN_FASTA="${2:-}"
            shift 2
            ;;
        --score-threshold)
            SCORE_THRESHOLD="${2:-}"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown argument: $1"
            ;;
    esac
done

[[ -n "$MANIFEST" ]] || {
    usage
    die "--manifest is required"
}
[[ -n "$OUTDIR" ]] || {
    usage
    die "--outdir is required"
}
[[ -f "$MANIFEST" ]] || die "Manifest file not found: $MANIFEST"

require_cmd jgi_summarize_bam_contig_depths
require_cmd metabat2
require_cmd SemiBin2
require_cmd vamb
require_cmd DAS_Tool
require_cmd awk
require_cmd date

mkdir -p "$OUTDIR"
MANIFEST="$(realpath_compat "$MANIFEST")"
OUTDIR="$(realpath_compat "$OUTDIR")"

log "MAG binning pipeline started"
log "Manifest: $MANIFEST"
log "Output directory: $OUTDIR"

{
    printf 'sample\tstatus\toutput_dir\n'

    local_line_no=0
    while IFS=$'\t' read -r sample contigs bam_files _rest; do
        local_line_no=$((local_line_no + 1))

        [[ -n "${sample:-}" ]] || continue
        if [[ "$local_line_no" -eq 1 ]]; then
            case "$sample" in
                sample|Sample|SAMPLE) continue ;;
            esac
        fi

        [[ -n "${contigs:-}" ]] || die "Line $local_line_no: missing contigs path"
        [[ -n "${bam_files:-}" ]] || die "Line $local_line_no: missing bam_files column"
        [[ -f "$contigs" ]] || die "[$sample] Contig FASTA not found: $contigs"

        sample="$(printf '%s' "$sample" | tr ' /' '__')"
        contigs="$(realpath_compat "$contigs")"

        sample_out="$OUTDIR/$sample"
        mkdir -p "$sample_out"

        log "[$sample] Processing sample"
        run_metabat2 "$sample" "$contigs" "$bam_files" "$sample_out"
        run_semibin2 "$sample" "$contigs" "$bam_files" "$sample_out"
        run_vamb "$sample" "$contigs" "$bam_files" "$sample_out"
        run_dastool "$sample" "$contigs" "$sample_out"
        printf '%s\t%s\t%s\n' "$sample" "done" "$sample_out"
    done <"$MANIFEST"
} >"$OUTDIR/run_summary.tsv"

log "MAG binning pipeline finished"
