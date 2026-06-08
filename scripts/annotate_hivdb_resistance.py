#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path


DEFAULT_HIVDB_XML = Path(
    "/Users/wuhhh/Desktop/徐老师/代码/metagenomic/database/virus/HIV/HIVDB_10.2.xml"
)
DEFAULT_HXB2_FASTA = Path(
    "/Users/wuhhh/Desktop/徐老师/代码/metagenomic/database/virus/HIV/HXB2_K03455.fasta"
)

GENE_REGIONS = {
    "CA": (1186, 1878),
    "PR": (2253, 2549),
    "RT": (2550, 4229),
    "IN": (4230, 5093),
}

AA_REFERENCE_SEQUENCES = {
    "PR": "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
    "RT": (
        "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGL"
        "KKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRT"
        "KIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYAGIKVKQLCKLLRGTKALTEVIPLTEEAE"
        "LELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKIATESIVIWGKTPKFKLPIQKETWEA"
        "WWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTDRGRQKVVSLTDTTNQKTELQAIHLALQDSGLEVNIVTDSQ"
        "YALGIIQAQPDKSESELVSQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVL"
    ),
    "IN": (
        "FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYF"
        "LLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERI"
        "VDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED"
    ),
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

GENE_ALIASES = {
    "CA": {"CA", "CAPSID"},
    "IN": {"IN", "INT", "INTEGRASE"},
    "PR": {"PR", "PRO", "PROTEASE"},
    "RT": {"RT", "POL", "REVERSETRANSCRIPTASE", "REVERSE_TRANSCRIPTASE"},
}

MUTATION_RE = re.compile(r"^(?:([A-Za-z*]))?(\d+)([A-Za-z*]+)$")


@dataclass
class LevelDefinition:
    order: int
    name: str
    sir: str


@dataclass
class CommentDefinition:
    id: str
    text: str
    sort_tag: str


@dataclass
class DrugResult:
    drug_class: str
    drug: str
    fullname: str
    score: int
    level: int
    level_name: str
    sir: str
    triggered_rules: list[str]
    result_comments: list[str]


@dataclass
class MutationComment:
    gene: str
    condition: str
    comment_id: str
    text: str


@dataclass
class ResultCommentRule:
    drug_name: str
    eq: int
    comment_id: str


@dataclass
class DrugDefinition:
    name: str
    fullname: str
    expression: str


@dataclass
class MutationRule:
    gene: str
    condition: str
    comment_id: str


@dataclass
class Algorithm:
    name: str
    version: str
    date: str
    default_level: int
    levels: dict[int, LevelDefinition]
    global_range: list[tuple[float, float, int]]
    drug_classes: dict[str, list[str]]
    gene_to_classes: dict[str, list[str]]
    drugs: list[DrugDefinition]
    comment_definitions: dict[str, CommentDefinition]
    mutation_rules: list[MutationRule]
    result_comment_rules: list[ResultCommentRule]


@dataclass
class FastaSampleResult:
    sample: str
    mutations: dict[str, dict[int, set[str]]]
    drug_results: list[DrugResult]
    mutation_comments: list[MutationComment]
    sequence_alerts: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Interpret HIV amino-acid mutations with the local Stanford HIVDB ASI XML "
            "and report per-drug resistance calls."
        )
    )
    parser.add_argument(
        "xml_or_mutation",
        nargs="?",
        help=(
            "Optional HIVDB ASI XML path. If omitted, the built-in default XML is used. "
            "If this token is not an XML path, it is treated as the first mutation group."
        ),
    )
    parser.add_argument(
        "mutations",
        nargs="*",
        help=(
            "Mutation groups such as 'RT:M184V,K65R' 'IN:G118R' 'PR:V82A' 'CA:Q67H,N74D'. "
            "Semicolon-delimited strings are also accepted."
        ),
    )
    parser.add_argument(
        "--mutation-file",
        type=Path,
        default=None,
        help="Optional text file containing mutation groups, one per line.",
    )
    parser.add_argument(
        "--xml",
        type=Path,
        default=None,
        help="Path to HIVDB ASI XML file.",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=None,
        help="Optional HIV consensus FASTA. When provided, amino-acid mutations are inferred automatically.",
    )
    parser.add_argument(
        "--hxb2-fasta",
        type=Path,
        default=DEFAULT_HXB2_FASTA,
        help="Local HXB2 reference FASTA used for nucleotide alignment.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit JSON instead of plain text.",
    )
    return parser.parse_args()


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:
        raise SystemExit(f"未找到 HIVDB XML 文件: {path}") from exc


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:].strip() or "unnamed"
                seq_parts = []
                continue
            seq_parts.append(re.sub(r"[^ACGTNacgtn-]", "", line))
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    if not records:
        raise SystemExit(f"FASTA 文件为空或格式无效: {path}")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def translate_codon(codon: str) -> str:
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3 or any(base not in {"A", "C", "G", "T"} for base in codon):
        return "X"
    return CODON_TABLE.get(codon, "X")


def translate_sequence(sequence: str, frame: int = 0) -> str:
    usable = sequence[frame:]
    limit = len(usable) - (len(usable) % 3)
    return "".join(
        translate_codon(usable[index:index + 3])
        for index in range(0, limit, 3)
    )


def resolve_mafft_path() -> str:
    candidates = [
        shutil.which("mafft"),
        "/opt/homebrew/Caskroom/mambaforge/base/bin/mafft",
        "/opt/homebrew/bin/mafft",
        "/usr/local/bin/mafft",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate))
    raise SystemExit("未找到 mafft，无法从 FASTA 自动推断 HIV 突变。")


def strip_xml_text(value: str) -> str:
    text = re.sub(r"<!\[CDATA\[(.*?)\]\]>", r"\1", value, flags=re.S)
    return re.sub(r"\s+", " ", text).strip()


def extract_tag_text(block: str, tag: str) -> str:
    match = re.search(rf"<{tag}>(.*?)</{tag}>", block, flags=re.S)
    if not match:
        raise ValueError(f"Missing <{tag}> block")
    return strip_xml_text(match.group(1))


def extract_section(text: str, tag: str) -> str:
    match = re.search(rf"<{tag}>(.*?)</{tag}>", text, flags=re.S)
    if not match:
        raise ValueError(f"Missing <{tag}> section")
    return match.group(1)


def parse_global_range(text: str) -> list[tuple[float, float, int]]:
    ranges: list[tuple[float, float, int]] = []
    for piece in text.strip().strip("()").split(","):
        piece = piece.strip()
        if not piece:
            continue
        part_match = re.match(r"(.+?)\s+TO\s+(.+?)\s*=>\s*(\d+)$", piece)
        if not part_match:
            raise ValueError(f"Unsupported global range entry: {piece}")
        lower_raw, upper_raw, level_raw = part_match.groups()
        lower = float("-inf") if lower_raw == "-INF" else float(lower_raw)
        upper = float("inf") if upper_raw == "INF" else float(upper_raw)
        ranges.append((lower, upper, int(level_raw)))
    return ranges


def parse_algorithm(xml_path: Path) -> Algorithm:
    text = read_text(xml_path)
    definitions = extract_section(text, "DEFINITIONS")

    gene_to_classes: dict[str, list[str]] = {}
    for block in re.findall(r"<GENE_DEFINITION>(.*?)</GENE_DEFINITION>", definitions, flags=re.S):
        gene = extract_tag_text(block, "NAME")
        classes = [item.strip() for item in extract_tag_text(block, "DRUGCLASSLIST").split(",") if item.strip()]
        gene_to_classes[gene] = classes

    levels: dict[int, LevelDefinition] = {}
    for block in re.findall(r"<LEVEL_DEFINITION>(.*?)</LEVEL_DEFINITION>", definitions, flags=re.S):
        level = LevelDefinition(
            order=int(extract_tag_text(block, "ORDER")),
            name=extract_tag_text(block, "ORIGINAL"),
            sir=extract_tag_text(block, "SIR"),
        )
        levels[level.order] = level

    drug_classes: dict[str, list[str]] = {}
    for match in re.finditer(
        r"<DRUGCLASS>\s*<NAME>(.*?)</NAME>\s*<DRUGLIST>(.*?)</DRUGLIST>\s*</DRUGCLASS>",
        definitions,
        flags=re.S,
    ):
        class_name = strip_xml_text(match.group(1))
        drugs = [item.strip() for item in strip_xml_text(match.group(2)).split(",") if item.strip()]
        drug_classes[class_name] = drugs

    comments_section = extract_section(definitions, "COMMENT_DEFINITIONS")
    comment_definitions: dict[str, CommentDefinition] = {}
    for match in re.finditer(
        r'<COMMENT_STRING id="([^"]+)">(.*?)</COMMENT_STRING>',
        comments_section,
        flags=re.S,
    ):
        comment_id, block = match.groups()
        comment_definitions[comment_id] = CommentDefinition(
            id=comment_id,
            text=extract_tag_text(block, "TEXT"),
            sort_tag=extract_tag_text(block, "SORT_TAG"),
        )

    drugs: list[DrugDefinition] = []
    drug_region = text.split("<MUTATION_COMMENTS>", 1)[0]
    for block in re.findall(r"<DRUG>(.*?)</DRUG>", drug_region, flags=re.S):
        expression = extract_tag_text(block, "CONDITION")
        if not expression.startswith("SCORE FROM"):
            continue
        drugs.append(
            DrugDefinition(
                name=extract_tag_text(block, "NAME"),
                fullname=extract_tag_text(block, "FULLNAME"),
                expression=expression,
            )
        )

    mutation_rules: list[MutationRule] = []
    mutation_section = extract_section(text, "MUTATION_COMMENTS")
    for gene_block in re.findall(r"<GENE>(.*?)</GENE>", mutation_section, flags=re.S):
        gene = extract_tag_text(gene_block, "NAME")
        for rule_block in re.findall(r"<RULE>(.*?)</RULE>", gene_block, flags=re.S):
            condition = extract_tag_text(rule_block, "CONDITION")
            comment_id = re.search(r'<COMMENT ref="([^"]+)"', rule_block)
            if not comment_id:
                continue
            mutation_rules.append(
                MutationRule(gene=gene, condition=condition, comment_id=comment_id.group(1))
            )

    result_comment_rules: list[ResultCommentRule] = []
    result_section = extract_section(text, "RESULT_COMMENTS")
    for rule_block in re.findall(r"<RESULT_COMMENT_RULE>(.*?)</RESULT_COMMENT_RULE>", result_section, flags=re.S):
        drug_name = extract_tag_text(rule_block, "DRUG_NAME")
        eq = int(extract_tag_text(rule_block, "EQ"))
        comment_match = re.search(r'<COMMENT ref="([^"]+)"', rule_block)
        if not comment_match:
            continue
        result_comment_rules.append(
            ResultCommentRule(drug_name=drug_name, eq=eq, comment_id=comment_match.group(1))
        )

    return Algorithm(
        name=extract_tag_text(text, "ALGNAME"),
        version=extract_tag_text(text, "ALGVERSION"),
        date=extract_tag_text(text, "ALGDATE"),
        default_level=int(extract_tag_text(definitions, "DEFAULT_LEVEL")),
        levels=levels,
        global_range=parse_global_range(extract_tag_text(definitions, "GLOBALRANGE")),
        drug_classes=drug_classes,
        gene_to_classes=gene_to_classes,
        drugs=drugs,
        comment_definitions=comment_definitions,
        mutation_rules=mutation_rules,
        result_comment_rules=result_comment_rules,
    )


def normalize_gene_name(raw: str) -> str:
    token = re.sub(r"[^A-Za-z]", "", raw).upper()
    for canonical, aliases in GENE_ALIASES.items():
        if token in aliases:
            return canonical
    raise ValueError(f"无法识别基因名称: {raw}")


def normalize_alt_token(raw_alt: str, *, has_reference_prefix: bool) -> set[str]:
    token_upper = raw_alt.upper()
    token_lower = raw_alt.lower()
    if raw_alt == "i":
        return {"i"}
    if raw_alt == "d":
        return {"d"}
    if token_lower == "ins":
        return {"i"}
    if token_lower == "del":
        return {"d"}
    if not has_reference_prefix and token_lower == "d":
        return {"d"}
    if not has_reference_prefix and token_lower == "i":
        return {"i"}
    return set(token_upper)


def parse_mutation_token(token: str) -> tuple[int, set[str]]:
    cleaned = token.strip().replace("/", "").replace(".", "")
    match = MUTATION_RE.match(cleaned)
    if not match:
        raise ValueError(f"无法解析突变: {token}")
    has_reference_prefix = bool(match.group(1))
    position = int(match.group(2))
    alts = normalize_alt_token(match.group(3), has_reference_prefix=has_reference_prefix)
    return position, alts


def collect_input_groups(args: argparse.Namespace) -> list[str]:
    groups: list[str] = []
    if args.xml_or_mutation and not str(args.xml_or_mutation).lower().endswith(".xml") and not is_fasta_like(str(args.xml_or_mutation)):
        groups.extend(part.strip() for part in args.xml_or_mutation.split(";") if part.strip())
    mutation_items = args.mutations
    if mutation_items and is_fasta_like(str(mutation_items[0])) and not args.fasta:
        mutation_items = mutation_items[1:]
    for item in mutation_items:
        groups.extend(part.strip() for part in item.split(";") if part.strip())
    if args.mutation_file:
        file_text = args.mutation_file.read_text(encoding="utf-8")
        for line in file_text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            groups.extend(part.strip() for part in line.split(";") if part.strip())
    if not groups:
        raise SystemExit("请至少提供一组突变，例如: RT:M184V,K65R IN:G118R")
    return groups


def resolve_xml_path(args: argparse.Namespace) -> Path:
    if args.xml:
        return args.xml
    if args.xml_or_mutation and str(args.xml_or_mutation).lower().endswith(".xml"):
        return Path(args.xml_or_mutation)
    return DEFAULT_HIVDB_XML


def is_fasta_like(path_text: str) -> bool:
    return path_text.lower().endswith((".fa", ".fasta", ".fna"))


def resolve_fasta_path(args: argparse.Namespace) -> Path | None:
    if args.fasta:
        return args.fasta
    if args.xml_or_mutation and is_fasta_like(str(args.xml_or_mutation)):
        return Path(args.xml_or_mutation)
    if args.mutations and is_fasta_like(str(args.mutations[0])) and not args.mutation_file:
        return Path(args.mutations[0])
    return None


def run_mafft_alignment(reference_fasta: Path, query_fasta: Path) -> dict[str, str]:
    mafft = resolve_mafft_path()
    result = subprocess.run(
        [mafft, "--quiet", "--addfragments", str(query_fasta), str(reference_fasta)],
        check=True,
        capture_output=True,
        text=True,
    )
    records: dict[str, str] = {}
    header: str | None = None
    seq_parts: list[str] = []
    for raw_line in result.stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(seq_parts).upper()
            header = line[1:].strip()
            seq_parts = []
            continue
        seq_parts.append(line)
    if header is not None:
        records[header] = "".join(seq_parts).upper()
    return records


def build_reference_position_map(aligned_reference: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    ref_position = 0
    for column_index, base in enumerate(aligned_reference):
        if base != "-":
            ref_position += 1
            mapping[ref_position] = column_index
    return mapping


def extract_gene_alignment_slice(
    aligned_reference: str,
    aligned_query: str,
    start_nt: int,
    end_nt: int,
) -> tuple[str, str]:
    position_map = build_reference_position_map(aligned_reference)
    gene_start_col = position_map[start_nt]
    gene_end_col = position_map[end_nt]
    return (
        aligned_reference[gene_start_col:gene_end_col + 1],
        aligned_query[gene_start_col:gene_end_col + 1],
    )


def infer_gene_mutations_and_alerts(
    gene: str,
    ref_slice: str,
    query_slice: str,
) -> tuple[list[str], list[str]]:
    ref_nt = "".join(base for base in ref_slice if base != "-")
    ref_aa = AA_REFERENCE_SEQUENCES.get(gene, translate_sequence(ref_nt, 0))
    mutations: list[str] = []
    alerts: list[str] = []
    seen_mutations: set[str] = set()

    codon_query_bases: list[str] = []
    codon_ref_bases: list[str] = []
    pending_insertion_nts = ""
    aa_position = 0

    def flush_insertion(before_position: int) -> None:
        nonlocal pending_insertion_nts
        if not pending_insertion_nts or before_position <= 0:
            pending_insertion_nts = ""
            return
        inserted_length = len(pending_insertion_nts)
        if inserted_length % 3 == 0:
            token = f"{before_position}i"
            if token not in seen_mutations:
                mutations.append(token)
                seen_mutations.add(token)
        else:
            alerts.append(f"{gene}:{before_position}ins{inserted_length}bp")
        pending_insertion_nts = ""

    for ref_base, query_base in zip(ref_slice, query_slice):
        if ref_base == "-":
            if query_base != "-":
                pending_insertion_nts += query_base
            continue

        codon_ref_bases.append(ref_base)
        codon_query_bases.append(query_base)
        if len(codon_ref_bases) < 3:
            continue

        aa_position += 1
        flush_insertion(aa_position - 1)

        ref_aa_base = ref_aa[aa_position - 1] if aa_position - 1 < len(ref_aa) else translate_codon("".join(codon_ref_bases))
        query_codon = "".join(codon_query_bases)
        if query_codon == "---":
            token = f"{aa_position}d"
            if token not in seen_mutations:
                mutations.append(token)
                seen_mutations.add(token)
        elif "-" in query_codon:
            alerts.append(f"{gene}:{aa_position}del{query_codon.count('-')}bp")
        else:
            query_aa_base = translate_codon(query_codon)
            if query_aa_base not in {"X", ref_aa_base}:
                token = f"{ref_aa_base}{aa_position}{query_aa_base}"
                if token not in seen_mutations:
                    mutations.append(token)
                    seen_mutations.add(token)

        codon_ref_bases = []
        codon_query_bases = []

    flush_insertion(aa_position)

    return mutations, alerts


def infer_mutations_from_fasta(
    fasta_path: Path,
    hxb2_fasta: Path,
) -> list[tuple[str, dict[str, dict[int, set[str]]], list[str]]]:
    reference_records = parse_fasta(hxb2_fasta)
    if len(reference_records) != 1:
        raise SystemExit(f"HXB2 参考 FASTA 应只包含一条序列: {hxb2_fasta}")

    query_records = parse_fasta(fasta_path)
    _ref_name, ref_sequence = reference_records[0]
    ref_name = "__HXB2_REF__"
    renamed_query_records = [
        (f"query_{index}", sequence)
        for index, (_query_name, sequence) in enumerate(query_records, start=1)
    ]
    query_name_lookup = {
        f"query_{index}": original_name
        for index, (original_name, _sequence) in enumerate(query_records, start=1)
    }

    with tempfile.TemporaryDirectory(prefix="hivdb_align_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        ref_path = tmpdir_path / "reference.fasta"
        query_path = tmpdir_path / "query.fasta"
        write_fasta(ref_path, [(ref_name, ref_sequence)])
        write_fasta(query_path, renamed_query_records)
        aligned_records = run_mafft_alignment(ref_path, query_path)

    if ref_name not in aligned_records:
        raise SystemExit("MAFFT 输出中缺少 HXB2 参考序列。")

    aligned_reference = aligned_records[ref_name]
    results: list[tuple[str, dict[str, dict[int, set[str]]], list[str]]] = []
    for aligned_query_name, _ in renamed_query_records:
        query_name = query_name_lookup[aligned_query_name]
        aligned_query = aligned_records.get(aligned_query_name)
        if not aligned_query:
            continue
        rendered_groups: list[str] = []
        sequence_alerts: list[str] = []
        for gene, (start_nt, end_nt) in GENE_REGIONS.items():
            ref_slice, query_slice = extract_gene_alignment_slice(
                aligned_reference,
                aligned_query,
                start_nt,
                end_nt,
            )
            gene_mutations, gene_alerts = infer_gene_mutations_and_alerts(
                gene,
                ref_slice,
                query_slice,
            )
            if gene_mutations:
                rendered_groups.append(f"{gene}:{','.join(gene_mutations)}")
            sequence_alerts.extend(gene_alerts)
        results.append(
            (
                query_name,
                parse_mutations(rendered_groups) if rendered_groups else {},
                sorted(set(sequence_alerts)),
            )
        )
    return results


def parse_mutations(groups: list[str]) -> dict[str, dict[int, set[str]]]:
    parsed: dict[str, dict[int, set[str]]] = {}
    for group in groups:
        if ":" not in group:
            raise SystemExit(f"突变分组缺少基因前缀: {group}")
        gene_raw, mutation_text = group.split(":", 1)
        gene = normalize_gene_name(gene_raw)
        gene_map = parsed.setdefault(gene, {})
        for item in re.split(r"[\s,]+", mutation_text.strip()):
            if not item:
                continue
            position, alts = parse_mutation_token(item)
            gene_map.setdefault(position, set()).update(alts)
    return parsed


def format_mutation_map(mutations: dict[str, dict[int, set[str]]]) -> dict[str, list[str]]:
    rendered: dict[str, list[str]] = {}
    for gene in sorted(mutations):
        values: list[str] = []
        ref_aa = AA_REFERENCE_SEQUENCES.get(gene, "")
        for position in sorted(mutations[gene]):
            alts = mutations[gene][position]
            letters = "".join(
                sorted(
                    alts,
                    key=lambda item: (item not in {"d", "i"}, item),
                )
            )
            if letters in {"d", "i"}:
                values.append(f"{position}{letters}")
                continue
            ref_prefix = ref_aa[position - 1] if 0 < position <= len(ref_aa) else ""
            values.append(f"{ref_prefix}{position}{letters}")
        rendered[gene] = values
    return rendered


def split_top_level(text: str, delimiter: str = ",") -> list[str]:
    parts: list[str] = []
    depth = 0
    start = 0
    for index, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == delimiter and depth == 0:
            parts.append(text[start:index].strip())
            start = index + 1
    tail = text[start:].strip()
    if tail:
        parts.append(tail)
    return parts


def split_boolean(expression: str, operator: str) -> list[str]:
    parts: list[str] = []
    needle = f" {operator} "
    start = 0
    while True:
        index = expression.find(needle, start)
        if index < 0:
            parts.append(expression[start:].strip())
            return parts
        parts.append(expression[start:index].strip())
        start = index + len(needle)


def match_condition_token(token: str, observed: dict[int, set[str]]) -> bool:
    token = token.strip()
    match = re.match(r"^(\d+)([A-Za-z]+)$", token)
    if not match:
        raise ValueError(f"Unsupported condition token: {token}")
    position = int(match.group(1))
    choices = match.group(2)
    observed_alts = observed.get(position, set())
    if not observed_alts:
        return False
    allowed = set("d" if char == "d" else "i" if char == "i" else char for char in choices)
    return bool(observed_alts & allowed)


def evaluate_condition(expression: str, observed: dict[int, set[str]]) -> bool:
    for or_part in split_boolean(expression, "OR"):
        and_parts = split_boolean(or_part, "AND")
        if all(match_condition_token(part, observed) for part in and_parts):
            return True
    return False


def parse_score_expression(expression: str) -> list[str]:
    prefix = "SCORE FROM ("
    if not expression.startswith(prefix) or not expression.endswith(")"):
        raise ValueError(f"Unsupported score expression: {expression}")
    inner = expression[len(prefix):-1].strip()
    return split_top_level(inner, ",")


def parse_scored_term(term: str) -> tuple[str, int]:
    condition, score_text = term.rsplit("=>", 1)
    return condition.strip(), int(score_text.strip())


def score_entry(entry: str, observed: dict[int, set[str]]) -> tuple[int, list[str]]:
    entry = entry.strip()
    if entry.startswith("MAX(") and entry.endswith(")"):
        options = split_top_level(entry[4:-1].strip(), ",")
        best_score = 0
        best_conditions: list[str] = []
        for option in options:
            condition, score = parse_scored_term(option)
            if evaluate_condition(condition, observed) and score >= best_score:
                if score > best_score:
                    best_conditions = [condition]
                    best_score = score
                else:
                    best_conditions.append(condition)
        return best_score, best_conditions
    condition, score = parse_scored_term(entry)
    if evaluate_condition(condition, observed):
        return score, [condition]
    return 0, []


def score_to_level(score: int, algorithm: Algorithm) -> int:
    for lower, upper, level in algorithm.global_range:
        if lower <= score <= upper:
            return level
    return algorithm.default_level


def infer_drug_class(drug_name: str, algorithm: Algorithm) -> str:
    for class_name, drug_names in algorithm.drug_classes.items():
        if drug_name in drug_names:
            return class_name
    return ""


def evaluate_drugs(
    algorithm: Algorithm,
    mutations: dict[str, dict[int, set[str]]],
) -> list[DrugResult]:
    results: list[DrugResult] = []
    for drug in algorithm.drugs:
        drug_class = infer_drug_class(drug.name, algorithm)
        relevant_genes = [
            gene for gene, classes in algorithm.gene_to_classes.items() if drug_class in classes
        ]
        observed: dict[int, set[str]] = {}
        for gene in relevant_genes:
            for position, alts in mutations.get(gene, {}).items():
                observed.setdefault(position, set()).update(alts)

        total_score = 0
        triggered: list[str] = []
        for entry in parse_score_expression(drug.expression):
            entry_score, matched_conditions = score_entry(entry, observed)
            total_score += entry_score
            triggered.extend(f"{condition} => {entry_score}" for condition in matched_conditions)

        level = score_to_level(total_score, algorithm)
        level_definition = algorithm.levels[level]
        result_comments: list[str] = []
        for rule in algorithm.result_comment_rules:
            if rule.drug_name == drug.name and rule.eq == level:
                result_comments.append(algorithm.comment_definitions[rule.comment_id].text)

        results.append(
            DrugResult(
                drug_class=drug_class,
                drug=drug.name,
                fullname=drug.fullname,
                score=total_score,
                level=level,
                level_name=level_definition.name,
                sir=level_definition.sir,
                triggered_rules=triggered,
                result_comments=result_comments,
            )
        )
    return results


def evaluate_mutation_comments(
    algorithm: Algorithm,
    mutations: dict[str, dict[int, set[str]]],
) -> list[MutationComment]:
    comments: list[MutationComment] = []
    seen: set[tuple[str, str]] = set()
    for rule in algorithm.mutation_rules:
        observed = mutations.get(rule.gene, {})
        if not observed:
            continue
        if not evaluate_condition(rule.condition, observed):
            continue
        key = (rule.gene, rule.comment_id)
        if key in seen:
            continue
        seen.add(key)
        comments.append(
            MutationComment(
                gene=rule.gene,
                condition=rule.condition,
                comment_id=rule.comment_id,
                text=algorithm.comment_definitions[rule.comment_id].text,
            )
        )
    comments.sort(key=lambda item: (item.gene, item.comment_id))
    return comments


def render_text(
    algorithm: Algorithm,
    mutations: dict[str, dict[int, set[str]]],
    drug_results: list[DrugResult],
    mutation_comments: list[MutationComment],
) -> str:
    lines: list[str] = []
    lines.append(f"{algorithm.name} {algorithm.version} ({algorithm.date})")
    lines.append("")
    lines.append("Input mutations:")
    for gene, items in format_mutation_map(mutations).items():
        lines.append(f"  {gene}: {', '.join(items)}")

    lines.append("")
    lines.append("Drug results:")
    for result in sorted(drug_results, key=lambda item: (item.drug_class, item.drug)):
        lines.append(
            f"  {result.drug:5} [{result.drug_class}] score={result.score:>3} "
            f"level={result.level} {result.level_name} ({result.sir})"
        )
        for comment in result.result_comments:
            lines.append(f"    note: {comment}")
        if result.triggered_rules:
            lines.append(f"    matched: {'; '.join(result.triggered_rules)}")

    if mutation_comments:
        lines.append("")
        lines.append("Mutation comments:")
        for item in mutation_comments:
            lines.append(f"  {item.gene} {item.condition}: {item.text}")

    return "\n".join(lines)


def render_json(
    algorithm: Algorithm,
    mutations: dict[str, dict[int, set[str]]],
    drug_results: list[DrugResult],
    mutation_comments: list[MutationComment],
) -> str:
    payload = {
        "algorithm": {
            "name": algorithm.name,
            "version": algorithm.version,
            "date": algorithm.date,
        },
        "input_mutations": format_mutation_map(mutations),
        "drug_results": [asdict(item) for item in drug_results],
        "mutation_comments": [asdict(item) for item in mutation_comments],
    }
    return json.dumps(payload, ensure_ascii=False, indent=2)


def render_fasta_text(
    algorithm: Algorithm,
    sample_results: list[FastaSampleResult],
) -> str:
    sections: list[str] = [f"{algorithm.name} {algorithm.version} ({algorithm.date})"]
    for sample_result in sample_results:
        sections.append("")
        sections.append(f"Sample: {sample_result.sample}")
        if not sample_result.mutations:
            sections.append("No HIVDB-covered mutations inferred from FASTA.")
            continue
        detailed_lines = render_text(
            algorithm,
            sample_result.mutations,
            sample_result.drug_results,
            sample_result.mutation_comments,
        )
        sections.extend(detailed_lines.splitlines()[2:])
        if sample_result.sequence_alerts:
            sections.append("")
            sections.append("Sequence alerts:")
            for alert in sample_result.sequence_alerts:
                sections.append(f"  {alert}")
    return "\n".join(sections)


def render_fasta_json(
    algorithm: Algorithm,
    sample_results: list[FastaSampleResult],
) -> str:
    payload = {
        "algorithm": {
            "name": algorithm.name,
            "version": algorithm.version,
            "date": algorithm.date,
        },
        "samples": [
            {
                "sample": sample_result.sample,
                "input_mutations": format_mutation_map(sample_result.mutations),
                "drug_results": [asdict(item) for item in sample_result.drug_results],
                "mutation_comments": [asdict(item) for item in sample_result.mutation_comments],
                "sequence_alerts": sample_result.sequence_alerts,
            }
            for sample_result in sample_results
        ],
    }
    return json.dumps(payload, ensure_ascii=False, indent=2)


def main() -> int:
    args = parse_args()
    algorithm = parse_algorithm(resolve_xml_path(args))
    fasta_path = resolve_fasta_path(args)

    if fasta_path is not None:
        sample_mutations = infer_mutations_from_fasta(fasta_path, args.hxb2_fasta)
        sample_results = [
            FastaSampleResult(
                sample=sample_name,
                mutations=mutations,
                drug_results=evaluate_drugs(algorithm, mutations),
                mutation_comments=evaluate_mutation_comments(algorithm, mutations),
                sequence_alerts=sequence_alerts,
            )
            for sample_name, mutations, sequence_alerts in sample_mutations
        ]
        output = (
            render_fasta_json(algorithm, sample_results)
            if args.json
            else render_fasta_text(algorithm, sample_results)
        )
    else:
        mutations = parse_mutations(collect_input_groups(args))
        drug_results = evaluate_drugs(algorithm, mutations)
        mutation_comments = evaluate_mutation_comments(algorithm, mutations)
        output = (
            render_json(algorithm, mutations, drug_results, mutation_comments)
            if args.json
            else render_text(algorithm, mutations, drug_results, mutation_comments)
        )
    print(output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
