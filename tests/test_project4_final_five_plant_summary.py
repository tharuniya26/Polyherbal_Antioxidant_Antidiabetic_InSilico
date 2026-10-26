"""Integrity tests for the deterministic Project 4 five-plant synthesis."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/build_project4_final_five_plant_summary.py"
SPEC = importlib.util.spec_from_file_location("project4_summary", SCRIPT)
assert SPEC and SPEC.loader
SUMMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMMARY)

FINAL_CSV = ROOT / "data/metadata/project4_final_five_plant_docking_summary.csv"
FINAL_REPORT = ROOT / "docs/protocols/project4_final_five_plant_scientific_synthesis.md"
EXPECTED_SHA256 = "5e3168dd3e5ed289e53eb8b48922a4c809a942f6b80a8dd7f294f184b20da3a1"
EXPECTED_SCHEMA = (
    "plant_associations", "botanical_evidence_tier", "structure_id", "compound_name",
    "ligand_microstate", "formal_charge", "receptor_state", "docking_box_A",
    "experiment_context", "analyzable_completed_searches", "exit_zero_searches",
    "exit_unconfirmed_complete_outputs", "returned_poses", "rank1_pocket_localized",
    "rank1_full_geometric_passes", "rank1_severe_clash_count",
    "rank1_vina_score_mean_kcal_mol", "rank1_vina_score_median_kcal_mol",
    "rank1_core_RMSD_to_medoid_min_A", "rank1_core_RMSD_to_medoid_max_A",
    "geometric_classification", "water_overlap_diagnostic", "docking_box_sensitivity",
    "main_limitations", "botanical_detection_boundary", "shared_experiment_note",
    "source_results_csv", "source_all_pose_csv", "source_manifest",
    "source_results_sha256", "source_all_pose_sha256", "source_manifest_sha256",
    "source_columns", "water_diagnostic_provenance", "botanical_evidence_provenance",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        return list(reader.fieldnames or ()), list(reader)


class FinalFivePlantSummaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.schema, cls.rows = read_csv(FINAL_CSV)

    def test_01_final_schema_row_count_and_identities(self) -> None:
        self.assertEqual(tuple(self.schema), EXPECTED_SCHEMA)
        self.assertEqual(len(self.rows), 26)
        identities = {
            (r["structure_id"], r["ligand_microstate"], r["receptor_state"], r["docking_box_A"])
            for r in self.rows
        }
        self.assertEqual(len(identities), 26)
        expected_counts = {"LIG01": 4, "LIG02": 2, "LIG04": 8, "LIG06": 8, "LIG08": 2, "LIG09": 2}
        self.assertEqual(
            {key: sum(r["structure_id"] == key for r in self.rows) for key in expected_counts},
            expected_counts,
        )

    def test_02_deterministic_rebuild_and_expected_hash(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            first = Path(temp_dir) / "first.csv"
            second = Path(temp_dir) / "second.csv"
            SUMMARY.write_summary(first)
            SUMMARY.write_summary(second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(first.read_bytes(), FINAL_CSV.read_bytes())
            self.assertEqual(sha256(first), EXPECTED_SHA256)

    def test_03_summary_source_provenance(self) -> None:
        for row in self.rows:
            for path_column, hash_column in (
                ("source_results_csv", "source_results_sha256"),
                ("source_all_pose_csv", "source_all_pose_sha256"),
                ("source_manifest", "source_manifest_sha256"),
            ):
                source = ROOT / row[path_column]
                self.assertTrue(source.is_file(), source)
                self.assertEqual(sha256(source), row[hash_column], source)

    def test_04_archived_candidate_outputs_and_attempt_records(self) -> None:
        result_rows: list[dict[str, str]] = []
        pose_count = 0
        for source in SUMMARY.SOURCES.values():
            _, rows = read_csv(ROOT / source["results"])
            result_rows.extend(rows)
            _, poses = read_csv(ROOT / source["poses"])
            pose_count += len(poses)
        self.assertEqual(len(result_rows), 130)
        self.assertEqual(pose_count, 2559)
        self.assertEqual(sum(row["exit_code"] == "0" for row in result_rows), 126)
        self.assertEqual(sum(row["exit_code"] == "unconfirmed" for row in result_rows), 4)

        for row in result_rows:
            for path_column, hash_column in (
                ("receptor_path", "receptor_sha256"),
                ("ligand_path", "ligand_sha256"),
                ("source_sdf_path", "source_sdf_sha256"),
                ("crystal_reference_path", "crystal_reference_sha256"),
                ("start_path", "start_sha256"),
                ("output_path", "output_sha256"),
            ):
                path = ROOT / row[path_column]
                self.assertTrue(path.is_file(), path)
                self.assertEqual(sha256(path), row[hash_column], path)
            output_text = (ROOT / row["output_path"]).read_text(encoding="utf-8", errors="replace")
            self.assertEqual(sum(line.startswith("MODEL ") for line in output_text.splitlines()), int(row["n_modes_returned"]))
            self.assertEqual(sum(line.startswith("ENDMDL") for line in output_text.splitlines()), int(row["n_modes_returned"]))
            log_paths = row.get("attempt_log_paths", row["log_path"]).split(";")
            log_hashes = row.get("attempt_log_hashes", row["log_sha256"]).split(";")
            self.assertEqual(len(log_paths), len(log_hashes))
            for relative, expected in zip(log_paths, log_hashes, strict=True):
                self.assertEqual(sha256(ROOT / relative), expected, relative)

    def test_05_native_control_outcomes_and_run_count(self) -> None:
        SUMMARY.validate_target_controls()
        _, controls = read_csv(ROOT / "data/metadata/vina_target_validation_summary.csv")
        self.assertEqual(len(controls), 8)
        self.assertEqual(sum(int(row["completed_runs"]) for row in controls), 40)
        self.assertEqual(sum(row["PDB_ID"] == "3L4W" and row["rank1_passes_at_2A"] == "5" for row in controls), 5)
        self.assertEqual(sum(row["PDB_ID"] == "4GQR" and row["rank1_passes_at_2A"] == "0" for row in controls), 3)
        expected_by_source: dict[str, int] = {}
        for control in controls:
            expected_by_source[control["source_results_csv"]] = (
                expected_by_source.get(control["source_results_csv"], 0) + int(control["completed_runs"])
            )
        for source, expected in expected_by_source.items():
            _, run_rows = read_csv(ROOT / source)
            self.assertEqual(len(run_rows), expected, source)

    def test_06_scientific_classifications_are_preserved(self) -> None:
        by_ligand = {}
        for row in self.rows:
            by_ligand.setdefault(row["structure_id"], []).append(row)
        self.assertTrue(all(r["geometric_classification"] == "INCONCLUSIVE" for r in by_ligand["LIG01"] + by_ligand["LIG02"] + by_ligand["LIG09"]))
        self.assertTrue(all(r["geometric_classification"] == "PASS" and r["rank1_full_geometric_passes"] == "5" for r in by_ligand["LIG08"]))
        self.assertGreater(len({r["geometric_classification"] for r in by_ligand["LIG04"]}), 1)
        self.assertGreater(len({r["geometric_classification"] for r in by_ligand["LIG06"]}), 1)

    def test_07_shared_quercetin_is_not_duplicated(self) -> None:
        quercetin = [row for row in self.rows if row["structure_id"] == "LIG06"]
        self.assertEqual(len(quercetin), 8)
        self.assertEqual({row["plant_associations"] for row in quercetin}, {"Boerhavia diffusa; Senna auriculata"})
        identities = {(r["ligand_microstate"], r["receptor_state"], r["docking_box_A"]) for r in quercetin}
        self.assertEqual(len(identities), 8)
        self.assertTrue(all("reused for two botanical associations" in r["shared_experiment_note"] for r in quercetin))

    def test_08_exit_unconfirmed_evidence_is_preserved(self) -> None:
        affected = [row for row in self.rows if int(row["exit_unconfirmed_complete_outputs"])]
        self.assertEqual(sum(int(row["exit_unconfirmed_complete_outputs"]) for row in affected), 4)
        self.assertEqual(
            {(r["structure_id"], r["ligand_microstate"], r["receptor_state"], r["docking_box_A"], r["exit_zero_searches"], r["exit_unconfirmed_complete_outputs"]) for r in affected},
            {
                ("LIG04", "neutral_keto_reference", "Asp542-", "24 x 18 x 18", "3", "2"),
                ("LIG04", "neutral_keto_reference", "Asp542H(OD2)", "30 x 24 x 24", "3", "2"),
            },
        )

    def test_09_no_unsupported_biological_efficacy_claims(self) -> None:
        report = FINAL_REPORT.read_text(encoding="utf-8").lower()
        self.assertIn("does not mean experimentally demonstrated binding, inhibition", report)
        self.assertIn("does not support claims", report)
        combined = report + "\n" + FINAL_CSV.read_text(encoding="utf-8").lower()
        for forbidden in (
            "docking proves inhibition", "docking confirms inhibition", "docking demonstrates inhibition",
            "proves antidiabetic activity", "confirms clinical efficacy", "demonstrates clinical efficacy",
            "demonstrates synergy", "proves synergy",
        ):
            self.assertNotIn(forbidden, combined)


if __name__ == "__main__":
    unittest.main()
