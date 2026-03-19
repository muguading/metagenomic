from __future__ import annotations

from pathlib import Path

from genome_db import GenomeManager


def main() -> None:
    demo_fasta = Path("example_genome.fasta")
    if not demo_fasta.exists():
        demo_fasta.write_text(">chr1\nATGCGTACGTAGCTAGCTAG\n", encoding="utf-8")

    manager = GenomeManager("sqlite:///example_genomes.sqlite3")

    added = manager.add_genome(
        genome_id="G001",
        sample_name="patient_001",
        species_name="Klebsiella pneumoniae",
        taxid=573,
        genome_file_path=str(demo_fasta),
        submitter="bioinfo_team",
        description="Reference assembly for validation",
    )
    print("Added:", added)

    updated = manager.update_genome(added["id"], operator="bioinfo_team", description="Updated description")
    print("Updated:", updated)

    fetched = manager.get_genome(added["id"], operator="bioinfo_team")
    print("Fetched:", fetched)

    search_result = manager.search_genomes(species_name="Klebsiella", page=1, page_size=10, operator="bioinfo_team")
    print("Search total:", search_result.total)
    print("Search items:", search_result.items)

    manager.delete_genome(added["id"], operator="bioinfo_team")
    print("Deleted G001")


if __name__ == "__main__":
    main()
