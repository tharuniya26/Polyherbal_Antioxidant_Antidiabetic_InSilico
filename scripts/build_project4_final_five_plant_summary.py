"""Build the final five-plant docking summary from verified run tables.

This script performs no docking or rescoring. It recomputes only the frozen
rank-1 geometric decision from archived run-level fields, checks the stored
classification, and writes a deterministic synthesis table.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from statistics import fmean, median


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/metadata/project4_final_five_plant_docking_summary.csv"

SOURCES = {
    "gallic_enlarged": {
        "results": "data/metadata/vina_3L4W_gallic_acid_pilot_results.csv",
        "poses": "docking/pilot/vina/3L4W_gallic_acid_pose_only/all_pose_metrics.csv",
        "manifest": "docking/pilot/vina/3L4W_gallic_acid_pose_only/experiment_manifest.json",
        "box": "30 x 24 x 24",
        "experiment": "enlarged-box gallic-acid pilot",
    },
    "alternative_enlarged": {
        "results": "data/metadata/vina_3L4W_alternative_ligand_pilot_results.csv",
        "poses": "docking/pilot/vina/3L4W_alternative_ligands_pose_only/all_pose_metrics.csv",
        "manifest": "docking/pilot/vina/3L4W_alternative_ligands_pose_only/experiment_manifest.json",
        "box": "30 x 24 x 24",
        "experiment": "enlarged-box rutin/luteolin pilot",
    },
    "phytochemical_original": {
        "results": "data/metadata/vina_3L4W_phytochemical_pilot_results.csv",
        "poses": "docking/pilot/vina/3L4W_phytochemical_pose_only/all_pose_metrics.csv",
        "manifest": "docking/pilot/vina/3L4W_phytochemical_pose_only/experiment_manifest.json",
        "box": "24 x 18 x 18",
        "experiment": "original-box mangiferin/quercetin pilot",
    },
    "phytochemical_enlarged": {
        "results": "data/metadata/vina_3L4W_box_adequacy_results.csv",
        "poses": "docking/pilot/vina/3L4W_box_adequacy_sensitivity/all_pose_metrics.csv",
        "manifest": "docking/pilot/vina/3L4W_box_adequacy_sensitivity/experiment_manifest.json",
        "box": "30 x 24 x 24",
        "experiment": "preregistered enlarged-box comparison",
    },
    "agnuside_enlarged": {
        "results": "data/metadata/vina_3L4W_agnuside_pilot_results.csv",
        "poses": "docking/pilot/vina/3L4W_agnuside_pose_only/all_pose_metrics.csv",
        "manifest": "docking/pilot/vina/3L4W_agnuside_pose_only/experiment_manifest.json",
        "box": "30 x 24 x 24",
        "experiment": "enlarged-box agnuside pilot",
    },
}

COMPOUND_METADATA = {
    "LIG01": {
        "plants": "Aegle marmelos",
        "tier": "B",
        "evidence": "docs/protocols/project4_alternative_compound_selection.md",
        "boundary": "Ethanol-leaf HPLC-PDA standard match; not detected in the original aqueous extract.",
        "shared": "Unique gallic-acid experiment for the Aegle association.",
    },
    "LIG02": {
        "plants": "Aegle marmelos",
        "tier": "B",
        "evidence": "docs/protocols/project4_alternative_ligand_input_qualification.md",
        "boundary": "Organic/aqueous-organic leaf literature; not detected in the original aqueous extract.",
        "shared": "Rutin is also literature-associated with Senna, but this synthesis assigns the completed pilot to the prespecified Aegle question and does not count another experiment.",
    },
    "LIG04": {
        "plants": "Salacia reticulata",
        "tier": "A",
        "evidence": "data/metadata/project4_five_plant_evidence_summary.csv",
        "boundary": "Published compound-specific water-leaf HPLC evidence; not detected in the original project extract.",
        "shared": "Unique mangiferin molecular experiment for the Salacia association.",
    },
    "LIG06": {
        "plants": "Boerhavia diffusa; Senna auriculata",
        "tier": "A with analytical qualification (Boerhavia); B (Senna)",
        "evidence": "docs/protocols/project4_alternative_compound_selection.md",
        "boundary": "Boerhavia boiled-water leaf HPLC-DAD assignment and Senna aqueous-acetone leaf isolation; neither establishes detection in the original project extracts.",
        "shared": "One archived LIG06 molecular experiment reused for two botanical associations; searches must not be summed across plants or described as a new Boerhavia experiment.",
    },
    "LIG08": {
        "plants": "Vitex negundo",
        "tier": "B",
        "evidence": "docs/protocols/project4_alternative_ligand_input_qualification.md",
        "boundary": "Organic-solvent leaf HPLC standard match; not detected in the original aqueous extract.",
        "shared": "Luteolin is also literature-associated with Senna, but no duplicate molecular experiment is counted.",
    },
    "LIG09": {
        "plants": "Vitex negundo",
        "tier": "B / abstract-only aqueous evidence",
        "evidence": "data/metadata/vina_remaining_plants_input_manifest.csv",
        "boundary": "Methanol-first leaf evidence plus an accessible water-study abstract; not detected in the original aqueous extract.",
        "shared": "Unique agnuside experiment for the Vitex association.",
    },
}

FORMAL_CHARGE = {
    "neutral_acid_reference": "0",
    "carboxylate_pH6p5_hypothesis": "-1",
    "neutral_reference": "0",
    "neutral_keto_reference": "0",
    "xanthone_6_OH_monoanion_hypothesis": "-1",
    "flavonol_4prime_OH_monoanion_hypothesis": "-1",
    "neutral_agnuside_reference": "0",
}

ORDER = {"LIG01": 1, "LIG02": 2, "LIG06": 3, "LIG04": 4, "LIG09": 5, "LIG08": 6}
MICROSTATE_ORDER = {
    "neutral_acid_reference": 1,
    "carboxylate_pH6p5_hypothesis": 2,
    "neutral_reference": 1,
    "neutral_keto_reference": 1,
    "xanthone_6_OH_monoanion_hypothesis": 2,
    "neutral_agnuside_reference": 1,
    "flavonol_4prime_OH_monoanion_hypothesis": 2,
}
SOURCE_COLUMNS = (
    "run_status;exit_code;n_modes_returned;rank1_score;rank1_pocket_member;"
    "rank1_severe_clash;rank1_core_RMSD_to_medoid_A;passing_rank1_of_5;"
    "combination_status;output_sha256"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(relative: str) -> list[dict[str, str]]:
    with (ROOT / relative).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def bool_field(value: str) -> bool:
    if value not in {"True", "False"}:
        raise ValueError(f"Expected boolean text, found {value!r}")
    return value == "True"


def water_diagnostic(structure: str, microstate: str, receptor: str) -> tuple[str, str]:
    if structure == "LIG01":
        rows = [
            row for row in read_csv("data/metadata/vina_3L4W_gallic_acid_off_pocket_investigation.csv")
            if row["microstate"] == microstate and row["receptor_state"] == receptor
        ]
        overlaps = sum(bool(row["rank1_original_waters_within1p8A"].strip()) for row in rows)
        alternatives = sum(bool(row["first_MIG_pocket_original_waters_within1p8A"].strip()) for row in rows)
        return (
            f"{overlaps}/5 rank-1 poses overlap deposited waters <1.8 A; {alternatives}/5 first pocket alternatives do.",
            "data/metadata/vina_3L4W_gallic_acid_off_pocket_investigation.csv:rank1_original_waters_within1p8A;first_MIG_pocket_original_waters_within1p8A",
        )
    if structure in {"LIG02", "LIG08"}:
        poses = [
            row for row in read_csv(SOURCES["alternative_enlarged"]["poses"])
            if row["structure_id"] == structure and row["receptor_state"] == receptor
        ]
        rank1 = [row for row in poses if row["pose_rank"] == "1"]
        return (
            f"{sum(bool_field(row['water_under_1p8']) for row in rank1)}/5 rank-1 and {sum(bool_field(row['water_under_1p8']) for row in poses)}/{len(poses)} all-mode poses overlap deposited waters <1.8 A.",
            f"{SOURCES['alternative_enlarged']['poses']}:water_under_1p8",
        )
    if structure == "LIG09":
        poses = [
            row for row in read_csv(SOURCES["agnuside_enlarged"]["poses"])
            if row["receptor_state"] == receptor
        ]
        rank1 = [row for row in poses if row["pose_rank"] == "1"]
        return (
            f"{sum(bool_field(row['water_under_1p8']) for row in rank1)}/5 rank-1 and {sum(bool_field(row['water_under_1p8']) for row in poses)}/{len(poses)} all-mode poses overlap deposited waters <1.8 A.",
            f"{SOURCES['agnuside_enlarged']['poses']}:water_under_1p8",
        )
    return (
        "NOT EVALUATED in the Step 58/59 pose tables; crystallographic waters were omitted from the dry receptors.",
        "not_available",
    )


def box_sensitivity(structure: str, microstate: str, receptor: str, box: str) -> str:
    if structure == "LIG04":
        mapping = {
            ("neutral_keto_reference", "Asp542-"): "PASS 5/5 in both boxes; same-seed cores shift 1.828-1.868 A.",
            ("neutral_keto_reference", "Asp542H(OD2)"): "PASS 5/5 original; INCONCLUSIVE 3/5 enlarged (3+2 core split).",
            ("xanthone_6_OH_monoanion_hypothesis", "Asp542-"): "INCONCLUSIVE 0/5 original; PASS 4/5 enlarged; strongly box-sensitive.",
            ("xanthone_6_OH_monoanion_hypothesis", "Asp542H(OD2)"): "PASS 5/5 in both boxes, but cores shift 2.433-2.449 A.",
        }
        return mapping[(microstate, receptor)]
    if structure == "LIG06":
        mapping = {
            ("neutral_keto_reference", "Asp542-"): "INCONCLUSIVE 0/5 in both boxes; stable off-rule placement.",
            ("neutral_keto_reference", "Asp542H(OD2)"): "INCONCLUSIVE 0/5 in both boxes; stable off-rule placement.",
            ("flavonol_4prime_OH_monoanion_hypothesis", "Asp542-"): "INCONCLUSIVE 0/5 in both boxes.",
            ("flavonol_4prime_OH_monoanion_hypothesis", "Asp542H(OD2)"): "PASS 5/5 in both boxes; same-seed shifts 0.005-0.074 A.",
        }
        return mapping[(microstate, receptor)]
    if structure == "LIG01":
        return "Only enlarged box tested; all rank-1 poses are at the Z- face, so box independence is not established."
    if structure == "LIG02":
        return "Only enlarged box tested; 3/5 rank-1 face warnings in each receptor state."
    if structure == "LIG08":
        return "Only enlarged box tested; no rank-1 face warning and minimum margin exceeds 4.3 A."
    if structure == "LIG09":
        return "Only enlarged box tested; one rank-1 face warning; no independent box comparison."
    raise KeyError(structure)


def main_limitations(structure: str) -> str:
    specific = {
        "LIG01": "Top-ranked off-pocket box-face site conflicts with six deposited water positions; lower-ranked pocket poses do not satisfy the rank-1 rule.",
        "LIG02": "Flexible 16-torsion glycoside; five conformer/seed starts are limited; rank-1 pocket localization fails despite lower-ranked pocket modes.",
        "LIG04": "Classifications and detailed cores depend on ligand state, Asp542 state, and box; two complete outputs per box experiment have unconfirmed process exits.",
        "LIG06": "Neutral state is consistently off-rule; monoanion PASS is restricted to Asp542H; two complete outputs per box experiment have unconfirmed process exits.",
        "LIG08": "Only the neutral reference was tested; phenolate hypotheses and physical bound microstate remain unresolved.",
        "LIG09": "Five conformers are confounded one-to-one with seeds; 13-torsion glycoside sampling is sparse and deposited-water overlap is extensive.",
    }
    return specific[structure] + " Dry/no-NAG rigid receptors omit crystallographic waters and use provisional protonation; Vina scores are not affinities."


def build_rows() -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []
    expected_group_counts = {"LIG01": 4, "LIG02": 2, "LIG04": 8, "LIG06": 8, "LIG08": 2, "LIG09": 2}
    observed_group_counts = {key: 0 for key in expected_group_counts}

    for source_name, spec in SOURCES.items():
        rows = read_csv(spec["results"])
        groups: dict[tuple[str, str, str], list[dict[str, str]]] = {}
        for row in rows:
            structure = row["structure_id"]
            if structure not in COMPOUND_METADATA:
                continue
            key = (structure, row["microstate"], row["receptor_state"])
            groups.setdefault(key, []).append(row)

        for (structure, microstate, receptor), group in groups.items():
            if len(group) != 5:
                raise RuntimeError(f"{source_name} {structure}/{microstate}/{receptor}: expected 5 rows, found {len(group)}")
            observed_group_counts[structure] += 1
            full_passes = sum(
                bool_field(row["rank1_pocket_member"])
                and not bool_field(row["rank1_severe_clash"])
                and float(row["rank1_core_RMSD_to_medoid_A"]) <= 2.0
                for row in group
            )
            decision = "PASS" if full_passes >= 4 else "INCONCLUSIVE"
            stored_decisions = {row["combination_status"] for row in group}
            stored_passes = {int(row["passing_rank1_of_5"]) for row in group}
            if stored_decisions != {decision} or stored_passes != {full_passes}:
                raise RuntimeError(
                    f"Frozen-rule mismatch for {source_name} {structure}/{microstate}/{receptor}: "
                    f"recomputed {decision} {full_passes}/5, stored {stored_decisions} {stored_passes}"
                )
            scores = [float(row["rank1_score"]) for row in group]
            rmsds = [float(row["rank1_core_RMSD_to_medoid_A"]) for row in group]
            exit_zero = sum(row["exit_code"] == "0" for row in group)
            exit_unconfirmed = sum(row["exit_code"] == "unconfirmed" for row in group)
            if exit_zero + exit_unconfirmed != 5:
                raise RuntimeError(f"Unexpected exit evidence for {source_name} {structure}/{microstate}/{receptor}")
            water_text, water_source = water_diagnostic(structure, microstate, receptor)
            metadata = COMPOUND_METADATA[structure]
            rows_out.append({
                "plant_associations": metadata["plants"],
                "botanical_evidence_tier": metadata["tier"],
                "structure_id": structure,
                "compound_name": group[0]["compound_name"],
                "ligand_microstate": microstate,
                "formal_charge": group[0].get("formal_charge") or group[0].get("state_formal_charge") or FORMAL_CHARGE[microstate],
                "receptor_state": receptor,
                "docking_box_A": spec["box"],
                "experiment_context": spec["experiment"],
                "analyzable_completed_searches": "5",
                "exit_zero_searches": str(exit_zero),
                "exit_unconfirmed_complete_outputs": str(exit_unconfirmed),
                "returned_poses": str(sum(int(row["n_modes_returned"]) for row in group)),
                "rank1_pocket_localized": str(sum(bool_field(row["rank1_pocket_member"]) for row in group)),
                "rank1_full_geometric_passes": str(full_passes),
                "rank1_severe_clash_count": str(sum(bool_field(row["rank1_severe_clash"]) for row in group)),
                "rank1_vina_score_mean_kcal_mol": f"{fmean(scores):.3f}",
                "rank1_vina_score_median_kcal_mol": f"{median(scores):.3f}",
                "rank1_core_RMSD_to_medoid_min_A": f"{min(rmsds):.4f}",
                "rank1_core_RMSD_to_medoid_max_A": f"{max(rmsds):.4f}",
                "geometric_classification": decision,
                "water_overlap_diagnostic": water_text,
                "docking_box_sensitivity": box_sensitivity(structure, microstate, receptor, spec["box"]),
                "main_limitations": main_limitations(structure),
                "botanical_detection_boundary": metadata["boundary"],
                "shared_experiment_note": metadata["shared"],
                "source_results_csv": spec["results"],
                "source_all_pose_csv": spec["poses"],
                "source_manifest": spec["manifest"],
                "source_results_sha256": digest(ROOT / spec["results"]),
                "source_all_pose_sha256": digest(ROOT / spec["poses"]),
                "source_manifest_sha256": digest(ROOT / spec["manifest"]),
                "source_columns": SOURCE_COLUMNS,
                "water_diagnostic_provenance": water_source,
                "botanical_evidence_provenance": metadata["evidence"],
            })

    if observed_group_counts != expected_group_counts:
        raise RuntimeError(f"Unexpected condition coverage: {observed_group_counts}")
    if len(rows_out) != 26:
        raise RuntimeError(f"Expected 26 unique molecular-condition rows, found {len(rows_out)}")
    rows_out.sort(key=lambda row: (
        ORDER[row["structure_id"]],
        MICROSTATE_ORDER[row["ligand_microstate"]],
        0 if row["receptor_state"] == "Asp542-" else 1,
        0 if row["docking_box_A"] == "24 x 18 x 18" else 1,
    ))
    return rows_out


def validate_target_controls() -> None:
    rows = read_csv("data/metadata/vina_target_validation_summary.csv")
    for row in rows:
        passes = int(row["rank1_passes_at_2A"])
        expected = "passed_rank1_pose_screen" if passes == 5 else "failed_rank1_pose_screen"
        if row["diagnostic_outcome"] != expected:
            raise RuntimeError(f"Target-control mismatch: {row['PDB_ID']} {row['experiment']}")
    if sum(row["PDB_ID"] == "3L4W" and int(row["rank1_passes_at_2A"]) == 5 for row in rows) != 5:
        raise RuntimeError("Expected five successful 3L4W native-control conditions")
    if sum(row["PDB_ID"] == "4GQR" and int(row["rank1_passes_at_2A"]) == 0 for row in rows) != 3:
        raise RuntimeError("Expected three failed 4GQR native-control conditions")


def write_summary(output: Path = OUTPUT) -> list[dict[str, str]]:
    """Validate archived inputs and write the deterministic 26-row summary."""
    validate_target_controls()
    rows = build_rows()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    rows = write_summary()
    print(f"Wrote {len(rows)} verified rows to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
