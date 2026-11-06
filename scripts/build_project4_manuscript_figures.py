"""Build the six Step 91 revised manuscript figures from archived Project 4 records.

This is an analysis-only renderer. It never invokes AutoDock Vina and never
changes a receptor, ligand, pose, run record, or frozen classification.
"""
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figures" / "project4_manuscript"
SUMMARY = ROOT / "data" / "metadata" / "project4_final_five_plant_docking_summary.csv"
ALT_RESULTS = ROOT / "data" / "metadata" / "vina_3L4W_alternative_ligand_pilot_results.csv"
ALT_POSES = ROOT / "docking" / "pilot" / "vina" / "3L4W_alternative_ligands_pose_only" / "all_pose_metrics.csv"
TARGETS = ROOT / "data" / "metadata" / "vina_target_validation_summary.csv"
REDOCK = ROOT / "data" / "metadata" / "vina_3L4W_redocking_results.csv"
LARGE_BOX_MIG = ROOT / "data" / "metadata" / "vina_3L4W_large_box_mig_results.csv"
BOX = ROOT / "data" / "metadata" / "vina_3L4W_box_adequacy_results.csv"
CRYSTAL_MIG = ROOT / "docking" / "validation" / "vina" / "3L4W_MIG_native_pose.sdf"
REDOCK_SEED_1101 = ROOT / "docking" / "validation" / "vina" / "3L4W_diagnostic_redocking" / "seed_1101_docked.pdbqt"
LUTEOLIN_RECEPTOR = ROOT / "structures" / "proteins" / "vina_prepared" / "3L4W_noNAG_D443m_D542m_H600HID_pH6p5_dry_sensitivity_Meeko0.8.0.pdbqt"
LUTEOLIN_SEED_1101 = ROOT / "docking" / "pilot" / "vina" / "3L4W_alternative_ligands_pose_only" / "LIG08_neutral_D542minus" / "seed_1101_docked.pdbqt"

PASS = "#147D73"
INC = "#D9882B"
NAVY = "#174A67"
BLUE = "#3C78A8"
PALE = "#E9F1F5"
GREY = "#66717A"
RED = "#B6423C"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "svg.hashsalt": "project4-step89",
    }
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(
        OUT / f"{stem}.svg",
        bbox_inches="tight",
        facecolor="white",
        metadata={"Date": None, "Creator": "Project 4 Step 91 reproducible figure script"},
    )
    plt.close(fig)


def figure_1_workflow() -> None:
    fig, ax = plt.subplots(figsize=(11.5, 5.4))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")
    stages = [
        (0.25, 3.55, 2.05, 1.25, "Botanical evidence", "5 plants; 6 unique\nreference structures"),
        (2.65, 3.55, 2.05, 1.25, "Input qualification", "Explicit ligand and\nreceptor microstates"),
        (5.05, 3.55, 2.05, 1.25, "Frozen searches", "130 candidate searches\n5 prespecified seeds"),
        (7.45, 3.55, 2.05, 1.25, "Pose audit", "2,559 poses; pocket,\nclash and convergence"),
        (9.85, 3.55, 1.9, 1.25, "Synthesis", "26 conditions\n10 PASS\n16 INCONCLUSIVE"),
    ]
    for i, (x, y, w, h, title, body) in enumerate(stages):
        ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03,rounding_size=0.08",
                                            facecolor=PALE if i < 4 else "#E8F4F1", edgecolor=NAVY, lw=1.4))
        ax.text(x + w / 2, y + 0.88, title, ha="center", va="center", weight="bold", color=NAVY)
        ax.text(x + w / 2, y + 0.39, body, ha="center", va="center", color="#26343D", linespacing=1.15,
                fontsize=8.4 if i == 4 else 9.0)
        if i < len(stages) - 1:
            ax.annotate("", xy=(stages[i + 1][0] - 0.08, y + h / 2), xytext=(x + w + 0.08, y + h / 2),
                        arrowprops=dict(arrowstyle="-|>", lw=1.5, color=GREY))
    controls = [
        (0.7, 1.25, 2.3, "Native controls", "40 registered searches\n+ 10 separate large-box MIG"),
        (3.3, 1.25, 2.3, "Frozen geometric rule", "At least 4/5 qualifying\nrank-1 poses"),
        (5.9, 1.25, 2.3, "Sensitivity analyses", "Receptor state, ligand\nstate and box"),
        (8.5, 1.25, 2.8, "Interpretation boundary", "Pose reproducibility—not affinity,\ninhibition or clinical efficacy"),
    ]
    for x, y, w, title, body in controls:
        ax.add_patch(patches.FancyBboxPatch((x, y), w, 1.05, boxstyle="round,pad=0.03,rounding_size=0.08",
                                            facecolor="white", edgecolor="#8BA4B2", lw=1.1))
        ax.text(x + w / 2, y + 0.72, title, ha="center", va="center", weight="bold", color=NAVY)
        ax.text(x + w / 2, y + 0.31, body, ha="center", va="center", fontsize=7.7, color="#34444E")
    ax.annotate("validation and prespecified decision rules", xy=(6, 3.45), xytext=(6, 2.52), ha="center",
                arrowprops=dict(arrowstyle="-|>", lw=1.2, color=GREY), fontsize=8.3, color=GREY)
    ax.set_title("Project 4 analysis workflow and evidential boundaries", pad=10, color="#20323D")
    save(fig, "figure_01_workflow")


def _rdkit_from_pdbqt(path: Path) -> Chem.Mol:
    made = [m for m in RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(path.read_text(encoding="utf-8"))) if m]
    if len(made) != 1:
        raise ValueError(f"Expected one molecule in {path}; found {len(made)}")
    return made[0]


def _heavy_coords(mol: Chem.Mol, conformer: int = 0) -> tuple[np.ndarray, list[int]]:
    conf = mol.GetConformer(conformer)
    ids = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.asarray([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in ids])
    return xyz, ids


def _project(*arrays: np.ndarray) -> list[np.ndarray]:
    all_xyz = np.vstack(arrays)
    center = all_xyz.mean(axis=0)
    _, _, vh = np.linalg.svd(all_xyz - center, full_matrices=False)
    basis = vh[:2].T
    return [(arr - center) @ basis for arr in arrays]


def _draw_molecule(ax: plt.Axes, mol: Chem.Mol, xy: np.ndarray, atom_ids: list[int], color: str,
                   label: str, alpha: float = 1.0, zorder: int = 3) -> None:
    index = {atom_id: n for n, atom_id in enumerate(atom_ids)}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in index and b in index:
            p, q = xy[index[a]], xy[index[b]]
            ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=1.4, alpha=alpha, zorder=zorder)
    ax.scatter(xy[:, 0], xy[:, 1], s=24, c=color, edgecolors="white", linewidths=0.45,
               alpha=alpha, label=label, zorder=zorder + 1)


def figure_2_native_validation() -> None:
    controls = rows(TARGETS)
    large_box = rows(LARGE_BOX_MIG)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    labels, mins, maxs, bests, colors = [], [], [], [], []
    for r in controls:
        pdb = r["PDB_ID"]
        if pdb == "3L4W":
            state = "D542H" if "D542H" in r["receptor_variant"] else "D542−"
            lig = "MIG+" if r["ligand_microstate"] == "amine_plus1" else "MIG"
            labels.append(f"3L4W {state}\n{lig}; {r['experiment'].replace('_', ' ')}")
            colors.append(PASS)
        else:
            prot = "Glu233H" if "Glu233H" in r["receptor_variant"] else "Glu233−"
            labels.append(f"4GQR {prot}\n{r['experiment'].replace('_', ' ')}")
            colors.append(INC)
        mins.append(float(r["rank1_RMSD_min_A"]))
        maxs.append(float(r["rank1_RMSD_max_A"]))
        bests.append(float(r["best_returned_RMSD_A"]))
    for state in ("Asp542-", "Asp542H(OD2)"):
        subset = [r for r in large_box if r["receptor_state"] == state]
        if len(subset) != 5:
            raise ValueError(f"Expected five enlarged-box MIG controls for {state}")
        labels.append(f"3L4W {state}\nseparate enlarged-box control")
        colors.append(PASS)
        mins.append(min(float(r["rank1_RMSD_A"]) for r in subset))
        maxs.append(max(float(r["rank1_RMSD_A"]) for r in subset))
        bests.append(min(float(r["best_pose_RMSD_A"]) for r in subset))
    y = np.arange(len(labels))
    for i in range(len(labels)):
        ax.plot([mins[i], maxs[i]], [i, i], color=colors[i], lw=5, solid_capstyle="round")
        ax.scatter([bests[i]], [i], facecolors="white", edgecolors=colors[i], marker="D", s=34, zorder=4)
    ax.axvline(2.0, ls="--", lw=1.2, color=RED, label="2.0 Å rank-1 threshold")
    ax.set_yticks(y, labels, fontsize=7.4)
    ax.invert_yaxis()
    ax.set_xlabel("RMSD to deposited native ligand (Å)")
    ax.set_title("A  Registered and separate enlarged-box controls")
    ax.grid(axis="x", color="#D8E0E4", lw=0.7)
    ax.legend(handles=[Line2D([0], [0], color=PASS, lw=5, label="Rank-1 range"),
                       Line2D([0], [0], marker="D", markerfacecolor="white", markeredgecolor=NAVY,
                              color="none", label="Best returned mode"),
                       Line2D([0], [0], color=RED, ls="--", label="2.0 Å threshold")],
              loc="lower right", fontsize=7.5, frameon=False)

    redock = rows(REDOCK)
    chosen = sorted(redock, key=lambda r: float(r["top_RMSD_A"]))[len(redock) // 2]
    docked = _rdkit_from_pdbqt(ROOT / chosen["output_path"])
    crystal = Chem.SDMolSupplier(str(CRYSTAL_MIG), removeHs=False)[0]
    if crystal is None:
        raise ValueError("Could not read native MIG reference")
    cxyz, cids = _heavy_coords(crystal)
    dxyz, dids = _heavy_coords(docked, 0)
    cxy, dxy = _project(cxyz, dxyz)
    ax = axes[1]
    _draw_molecule(ax, crystal, cxy, cids, "#30343B", "Deposited MIG")
    _draw_molecule(ax, docked, dxy, dids, BLUE, f"Recovered rank 1 (seed {chosen['seed']})", alpha=0.78)
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    ax.set_title("B  Representative 3L4W rank-1 recovery")
    ax.text(0.02, 0.02, f"PDB 3L4W · Asp542− receptor · neutral MIG\nseed {chosen['seed']}, rank 1 · RMSD {float(chosen['top_RMSD_A']):.4f} Å\northographic PCA view; coordinates not fitted for display",
            transform=ax.transAxes, fontsize=8.0, va="bottom", color="#33434C",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#BCC8CE", alpha=0.94))
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    fig.suptitle("Native-control performance and crystallographic-pose recovery", y=1.01, weight="bold")
    fig.tight_layout()
    save(fig, "figure_02_native_validation")


def _short_state(row: dict[str, str]) -> str:
    state = "D542H" if "H" in row["receptor_state"] else "D542−"
    micro = row["ligand_microstate"].replace("neutral_keto_reference", "neutral").replace("neutral_reference", "neutral")
    micro = micro.replace("neutral_acid_reference", "neutral acid").replace("neutral_agnuside_reference", "neutral")
    micro = micro.replace("carboxylate_pH6p5_hypothesis", "carboxylate").replace("xanthone_6_OH_monoanion_hypothesis", "6-O−")
    micro = micro.replace("flavonol_4prime_OH_monoanion_hypothesis", "4′-O−")
    box = row["docking_box_A"].replace(" x ", "×")
    return f"{micro} · {state} · {box} Å"


def figure_3_final_outcomes() -> None:
    data = rows(SUMMARY)
    fig, ax = plt.subplots(figsize=(10.5, 11.2))
    y = np.arange(len(data))
    passes = [int(r["rank1_full_geometric_passes"]) for r in data]
    colors = [PASS if r["geometric_classification"] == "PASS" else INC for r in data]
    ax.barh(y, passes, color=colors, height=0.66, edgecolor="white")
    ax.axvline(4, color=RED, ls="--", lw=1.1)
    labels = [f"{r['structure_id']} {r['compound_name']}  |  {_short_state(r)}" for r in data]
    ax.set_yticks(y, labels, fontsize=7.8)
    ax.invert_yaxis()
    ax.set_xlim(0, 5.75)
    ax.set_xticks(range(6))
    ax.set_xlabel("Rank-1 poses meeting all frozen geometric criteria (of 5 seeds)")
    ax.set_title("Final condition-level geometric classifications (26 prespecified rows)", pad=12)
    ax.grid(axis="x", color="#DEE5E8", lw=0.65)
    for i, r in enumerate(data):
        n = int(r["rank1_full_geometric_passes"])
        ax.text(n + 0.08, i, f"{r['geometric_classification']}  {n}/5", va="center", fontsize=7.6,
                color="#24333B", weight="bold" if r["geometric_classification"] == "PASS" else "normal")
    ax.legend(handles=[patches.Patch(color=PASS, label="PASS (≥4/5)"),
                       patches.Patch(color=INC, label="INCONCLUSIVE"),
                       Line2D([0], [0], color=RED, ls="--", label="Prespecified 4/5 threshold")],
              loc="upper right", frameon=False)
    fig.text(0.99, 0.008, "PASS is a reproducible pose-geometry classification, not evidence of affinity or inhibition.",
             ha="right", fontsize=8, color=GREY)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    save(fig, "figure_03_final_outcomes")


def _read_receptor_atoms(path: Path) -> list[dict[str, object]]:
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atom_name = line[12:16].strip()
            if atom_name.upper().startswith("H"):
                continue
            atoms.append({"atom": atom_name, "res": line[17:20].strip(), "chain": line[21].strip(),
                          "resid": int(line[22:26]), "xyz": np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])])})
        except ValueError:
            continue
    return atoms


def figure_4_luteolin_pose() -> None:
    data = [r for r in rows(ALT_RESULTS) if r["structure_id"] == "LIG08" and r["receptor_state"] == "Asp542-" and r["seed"] == "1101"]
    if len(data) != 1 or data[0]["medoid_seed"] != "1101":
        raise ValueError("Prespecified luteolin medoid row was not uniquely identified")
    row = data[0]
    ligand = _rdkit_from_pdbqt(ROOT / row["output_path"])
    lxyz, lids = _heavy_coords(ligand, 0)
    mig = Chem.SDMolSupplier(str(CRYSTAL_MIG), removeHs=False)[0]
    mxyz, mids = _heavy_coords(mig)
    receptor = _read_receptor_atoms(ROOT / row["receptor_path"])
    pocket = [a for a in receptor if min(np.linalg.norm(a["xyz"] - x) for x in lxyz) <= 4.5]
    pxyz = np.asarray([a["xyz"] for a in pocket])
    lxy, mxy, pxy = _project(lxyz, mxyz, pxyz)
    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    residue_colors = {443: "#B6423C", 542: "#8E5AAE", 600: "#2D6F9D"}
    other = [i for i, a in enumerate(pocket) if a["resid"] not in residue_colors]
    if other:
        ax.scatter(pxy[other, 0], pxy[other, 1], s=13, c="#C7D0D5", alpha=0.55, label="Other pocket atoms", zorder=1)
    for resid, color in residue_colors.items():
        indices = [i for i, a in enumerate(pocket) if a["resid"] == resid]
        if not indices:
            continue
        pts = pxy[indices]
        ax.scatter(pts[:, 0], pts[:, 1], s=34, c=color, alpha=0.75, edgecolors="white", linewidths=0.4, zorder=2)
        original = np.asarray([pocket[i]["xyz"] for i in indices])
        distance = min(np.linalg.norm(p - q) for p in original for q in lxyz)
        centroid = pts.mean(axis=0)
        ax.annotate(f"{pocket[indices[0]]['res'].title()}{resid}\nnearest heavy atoms {distance:.2f} Å",
                    xy=centroid, xytext=(centroid[0] + 1.0, centroid[1] + 0.8), fontsize=8, color=color,
                    arrowprops=dict(arrowstyle="-", color=color, lw=0.8))
    ax.scatter(mxy[:, 0], mxy[:, 1], s=17, facecolors="none", edgecolors="#59636A", linewidths=0.8,
               label="Deposited MIG pocket reference", zorder=2)
    _draw_molecule(ax, ligand, lxy, lids, "#E49A22", "Luteolin rank 1", zorder=4)
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    ax.set_title("Representative convergent luteolin pose in the 3L4W MIG pocket", pad=14)
    ax.text(0.01, 0.01, "PDB 3L4W · dry/no-NAG Asp443−/Asp542−/His600-HID receptor\nLIG08 neutral luteolin · seed 1101 · rank 1 (within-state medoid)\nNo hydrogen bonds are assigned; labels report nearest heavy-atom distances.\nCrystallographic waters were omitted from docking; rank-1 water-overlap diagnostic: no.",
            transform=ax.transAxes, fontsize=8.2, va="bottom", color="#31414A",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#BCC8CE", alpha=0.96))
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    fig.tight_layout()
    save(fig, "figure_04_luteolin_pose")


def figure_5_sensitivity() -> None:
    data = [r for r in rows(SUMMARY) if r["structure_id"] in {"LIG04", "LIG06"}]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.8), sharex=True)
    for ax, sid, title in zip(axes, ["LIG04", "LIG06"], ["Mangiferin", "Quercetin"], strict=True):
        subset = [r for r in data if r["structure_id"] == sid]
        subset.sort(key=lambda r: (r["ligand_microstate"], r["receptor_state"], r["docking_box_A"]))
        labels = [_short_state(r) for r in subset]
        vals = [int(r["rank1_full_geometric_passes"]) for r in subset]
        cols = [PASS if r["geometric_classification"] == "PASS" else INC for r in subset]
        y = np.arange(len(subset))
        ax.barh(y, vals, color=cols, height=0.62)
        ax.set_yticks(y, labels, fontsize=8)
        ax.invert_yaxis()
        ax.axvline(4, color=RED, ls="--", lw=1)
        ax.set_xlim(0, 5.8)
        ax.set_xticks(range(6))
        ax.set_title(title)
        ax.grid(axis="x", color="#DFE6E9", lw=0.65)
        for i, r in enumerate(subset):
            marker = "†" if int(r["exit_unconfirmed_complete_outputs"]) else ""
            ax.text(vals[i] + 0.08, i, f"{r['geometric_classification']} {vals[i]}/5{marker}", va="center", fontsize=7.7)
    axes[0].set_xlabel("Full geometric passes among five seeds")
    axes[1].set_xlabel("Full geometric passes among five seeds")
    fig.suptitle("Receptor-, ligand-state- and docking-box sensitivity", y=1.01, weight="bold")
    fig.text(0.5, 0.01, "† Two complete, chemistry-checked pose outputs in that condition have unconfirmed process exits.\nThe original 24×18×18 Å and enlarged 30×24×24 Å boxes are separate prespecified analyses.",
             ha="center", fontsize=8.2, color=GREY)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, "figure_05_mangiferin_quercetin_sensitivity")


def figure_6_rutin() -> None:
    data = [r for r in rows(ALT_RESULTS) if r["structure_id"] == "LIG02"]
    states = ["Asp542-", "Asp542H(OD2)"]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.4), sharex="col")
    for col, state in enumerate(states):
        subset = sorted([r for r in data if r["receptor_state"] == state], key=lambda r: int(r["seed"]))
        seeds = [r["seed"] for r in subset]
        distances = [float(r["rank1_core_centroid_MIG_distance_A"]) for r in subset]
        rmsd = [float(r["rank1_core_RMSD_to_medoid_A"]) for r in subset]
        pocket = [r["rank1_pocket_member"] == "True" for r in subset]
        water = [r["rank1_water_under_1p8"] == "True" for r in subset]
        ax = axes[0, col]
        for i, (d, ok, wet) in enumerate(zip(distances, pocket, water, strict=True)):
            ax.scatter(i, d, s=75, c=PASS if ok else INC, marker="o" if not wet else "s", edgecolor="#24333B", linewidth=0.6)
        ax.axhline(6.0, color=RED, ls="--", lw=1)
        ax.set_title(f"3L4W {state}")
        ax.set_ylabel("Core-centroid distance to MIG (Å)" if col == 0 else "")
        ax.grid(axis="y", color="#DFE6E9", lw=0.65)
        ax = axes[1, col]
        ax.scatter(range(5), rmsd, s=68, c=[PASS if x <= 2 else INC for x in rmsd], edgecolor="#24333B", linewidth=0.6)
        ax.axhline(2.0, color=RED, ls="--", lw=1)
        ax.set_xticks(range(5), seeds)
        ax.set_xlabel("Prespecified Vina seed")
        ax.set_ylabel("Unfitted core RMSD to medoid (Å)" if col == 0 else "")
        ax.grid(axis="y", color="#DFE6E9", lw=0.65)
    axes[0, 0].legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=PASS, markeredgecolor="#24333B", label="MIG-pocket member"),
                                  Line2D([0], [0], marker="o", color="none", markerfacecolor=INC, markeredgecolor="#24333B", label="Off-pocket"),
                                  Line2D([0], [0], marker="s", color="none", markerfacecolor="#BBBBBB", markeredgecolor="#24333B", label="Deposited-water overlap <1.8 Å"),
                                  Line2D([0], [0], color=RED, ls="--", label="Frozen threshold")], fontsize=7.5, frameon=False, loc="upper right")
    fig.suptitle("Rutin remains inconclusive across all prespecified seeds", y=1.01, weight="bold")
    fig.text(0.5, 0.012, "LIG02 neutral reference · enlarged 30×24×24 Å box · rank 1 shown for every seed.\nNeither receptor state reached the required four qualifying rank-1 poses; no favorable subset was selected.",
             ha="center", fontsize=8.4, color=GREY)
    fig.tight_layout(rect=(0, 0.065, 1, 1))
    save(fig, "figure_06_rutin_inconclusive")


def write_provenance() -> None:
    specs = [
        ("Figure 1", "figure_01_workflow.png;figure_01_workflow.svg", [SUMMARY, TARGETS, LARGE_BOX_MIG], "Complete documented workflow; counts distinguish the 40-search registered control set from 10 separate enlarged-box 3L4W MIG controls.", "No structural exemplar."),
        ("Figure 2", "figure_02_native_validation.png;figure_02_native_validation.svg", [TARGETS, LARGE_BOX_MIG, REDOCK, CRYSTAL_MIG, REDOCK_SEED_1101], "All eight registered conditions and two separate enlarged-box 3L4W MIG conditions shown; range bars are observed rank-1 minima to maxima, white diamonds are best returned modes, and no midpoint statistic is plotted. Structural panel uses the median rank-1 RMSD among the five initial 3L4W neutral-MIG runs.", "PDB 3L4W; neutral MIG; Asp542-; seed selected by median rank-1 RMSD; rank 1."),
        ("Figure 3", "figure_03_final_outcomes.png;figure_03_final_outcomes.svg", [SUMMARY], "All 26 final synthesis rows shown in authoritative CSV order.", "No structural exemplar."),
        ("Figure 4", "figure_04_luteolin_pose.png;figure_04_luteolin_pose.svg", [ALT_RESULTS, CRYSTAL_MIG, LUTEOLIN_RECEPTOR, LUTEOLIN_SEED_1101], "Within-condition medoid for neutral luteolin/Asp542-; medoid seed was frozen before figure production. Nearest-distance annotations use receptor and ligand heavy atoms only.", "PDB 3L4W; LIG08 neutral luteolin; Asp542-; seed 1101; rank 1."),
        ("Figure 5", "figure_05_mangiferin_quercetin_sensitivity.png;figure_05_mangiferin_quercetin_sensitivity.svg", [SUMMARY, BOX], "All 16 mangiferin/quercetin ligand-state × receptor-state × box rows shown.", "No single pose selected; four exit-unconfirmed complete outputs remain marked by condition."),
        ("Figure 6", "figure_06_rutin_inconclusive.png;figure_06_rutin_inconclusive.svg", [ALT_RESULTS, ALT_POSES], "Every rank-1 rutin result across five seeds and both receptor states shown.", "PDB 3L4W; LIG02 neutral rutin; Asp542- and Asp542H(OD2); seeds 1101/2202/3303/4404/5505; rank 1."),
    ]
    with (OUT / "figure_provenance.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["figure", "exported_files", "source_path", "source_sha256", "selection_rule", "structural_identifiers"])
        writer.writeheader()
        for fig, exports, sources, rule, identifiers in specs:
            for source in sources:
                writer.writerow({"figure": fig, "exported_files": exports, "source_path": source.relative_to(ROOT).as_posix(),
                                 "source_sha256": sha256(source), "selection_rule": rule, "structural_identifiers": identifiers})


def main() -> None:
    figure_1_workflow()
    figure_2_native_validation()
    figure_3_final_outcomes()
    figure_4_luteolin_pose()
    figure_5_sensitivity()
    figure_6_rutin()
    write_provenance()
    print(f"Wrote six PNG and six SVG figures plus provenance to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
