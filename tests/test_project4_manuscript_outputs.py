"""Focused integrity tests for the Step 91 revised manuscript and figures."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import re

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript" / "project4_insilico_manuscript_draft.md"
FIGDIR = ROOT / "results" / "figures" / "project4_manuscript"
PROVENANCE = FIGDIR / "figure_provenance.csv"
PDF = ROOT / "manuscript" / "Tharuniya_Anpalakan_Polyherbal_Computational_Manuscript.pdf"
PDF_BUILDER = ROOT / "scripts" / "build_project4_manuscript_pdf.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_manuscript_structure_and_word_count() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    body = text.split("## References", 1)[0]
    words = re.findall(r"\b[\w′’\-]+\b", body, flags=re.UNICODE)
    assert 4000 <= len(words) <= 6500
    for heading in ("## Abstract", "## Introduction", "## Methods", "## Results", "## Discussion", "## Limitations", "## Conclusion", "## Data and code availability", "## References"):
        assert heading in text
    correspondence = re.findall(r"^Correspondence:\s*(\S+)$", text, flags=re.MULTILINE)
    assert correspondence == ["tharuniya08@gmail.com"]


def test_exact_archived_counts_and_boundaries_are_reported() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for required in (
        "130 analyzable searches",
        "2,559 poses",
        "126 searches",
        "40 registered",
        "10 separate enlarged-box",
        "50 archived native-control calculations",
        "10 PASS and 16 INCONCLUSIVE",
        "four of five",
    ):
        assert required in text
    for method in ("OpenMM 8.6.1", "PROPKA 3.5.1", "ETKDGv3", "MMFF94s", "PCG64", "Meeko 0.8.0"):
        assert method in text
    lower = text.lower()
    for forbidden in ("docking proves", "docking confirms inhibition", "proves antidiabetic", "confirms clinical efficacy", "demonstrates synergy"):
        assert forbidden not in lower


def test_six_figures_have_png_and_svg_exports() -> None:
    pngs = sorted(FIGDIR.glob("figure_*.png"))
    svgs = sorted(FIGDIR.glob("figure_*.svg"))
    assert len(pngs) == 6
    assert len(svgs) == 6
    for png in pngs:
        with Image.open(png) as image:
            assert image.width >= 2000
            assert image.height >= 1000
    for svg in svgs:
        assert svg.stat().st_size > 10_000


def test_figure_provenance_sources_match_current_hashes() -> None:
    with PROVENANCE.open(encoding="utf-8", newline="") as stream:
        records = list(csv.DictReader(stream))
    assert {row["figure"] for row in records} == {f"Figure {n}" for n in range(1, 7)}
    for row in records:
        source = ROOT / row["source_path"]
        assert source.is_file(), source
        assert sha256(source) == row["source_sha256"], source
        assert row["selection_rule"]
    figure_2_sources = {row["source_path"] for row in records if row["figure"] == "Figure 2"}
    assert "data/metadata/vina_3L4W_large_box_mig_results.csv" in figure_2_sources
    assert "docking/validation/vina/3L4W_diagnostic_redocking/seed_1101_docked.pdbqt" in figure_2_sources
    figure_4_sources = {row["source_path"] for row in records if row["figure"] == "Figure 4"}
    assert any(path.startswith("structures/proteins/vina_prepared/") for path in figure_4_sources)


def test_structural_identity_and_no_hbond_overstatement() -> None:
    captions = (ROOT / "manuscript" / "figure_captions.md").read_text(encoding="utf-8")
    for required in ("PDB 3L4W", "seed 1101", "rank 1", "neutral LIG08", "Asp542−", "waters", "N-glycans"):
        assert required in captions
    assert "do not assign hydrogen bonds" in captions


def test_native_control_figure_uses_observed_statistics_only() -> None:
    script = (ROOT / "scripts" / "build_project4_manuscript_figures.py").read_text(encoding="utf-8")
    assert "(mins[i] + maxs[i]) / 2" not in script
    captions = (ROOT / "manuscript" / "figure_captions.md").read_text(encoding="utf-8")
    assert "No midpoint or mean is plotted" in captions


def test_reference_register_is_complete_and_verified() -> None:
    path = ROOT / "manuscript" / "reference_verification.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        records = list(csv.DictReader(stream))
    assert len(records) == 27
    assert all(row["verification_status"] == "VERIFIED" for row in records)
    assert len({row["doi_or_url"].lower() for row in records}) == 27


def test_publication_pdf_and_reproducible_builder_exist() -> None:
    assert PDF_BUILDER.is_file()
    content = PDF.read_bytes()
    assert content.startswith(b"%PDF-")
    assert content.rstrip().endswith(b"%%EOF")
    assert len(content) > 500_000


def test_readme_links_final_manuscript_and_synthesis() -> None:
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    for required in (
        "Tharuniya_Anpalakan_Polyherbal_Computational_Manuscript.pdf",
        "10 PASS and 16 INCONCLUSIVE",
        "2,559 returned poses",
        "neutral luteolin passed the geometric screen **5/5 in both receptor states**",
        "not peer-reviewed, accepted or published",
    ):
        assert required in text


def test_publication_availability_wording_is_current() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    assert "This release contains the revised manuscript" in text
    assert "local untracked release-candidate artifacts" not in text
