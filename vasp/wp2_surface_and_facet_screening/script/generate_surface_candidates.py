"""Generate only WP2 Stage 0 structures and registries from frozen WP1 parents.

Run with the existing .venv-wp1-plots interpreter. No import-time writes occur.
Existing different output bytes cause an error; regeneration can be checked in
a temporary directory with --output-root. Nothing under calculation/ is written.
"""

import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from importlib.metadata import version
from pathlib import Path

import numpy as np
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp import Poscar

from surface_model_utils import (
    FACE_DEPTH, PLANE_TOL, VACUUM, bond_graph, composition_record, connectivity,
    coordination_reference, cut_diagnostics, diagnose_slab, present_slab, slabs_match,
)

WP2 = Path(__file__).resolve().parents[1]
REPO = WP2.parents[1]
WP1 = REPO / "vasp/wp1_polymorph_polytype_benchmark"
AUDIT = WP1 / "results/01_geometry_optimisation/FINAL_GEOMETRY_OPT_AUDIT.md"
FACETS = {"beta": ("001", "100", "110"), "IIb": ("001", "100"), "spinel": ("001", "110")}
ACCEPTED = {
    "beta": dict(sha256="e52ccce041a7efac9eabdbdd11b26f9dfe29b5be8d1047d7dddceef691c8480d",
                 units=1, a=3.878284, c=12.148433, gamma=120, group="P3m1 (156)"),
    "IIb": dict(sha256="2ec9c4fefc2d97f83ba94f70e1b048061fcbed0bd2efd8b8db5e1e59d23a1023",
                units=2, a=3.883347, c=24.351774, gamma=120, group="P6_3mc (186)"),
    "spinel": dict(sha256="89fed03927ad1d221915383fbfd61927e9d16798c1f1597737508cc9fa19ab05",
                   units=8, a=10.617497, c=10.617497, gamma=90, group="Fd-3m (227)"),
}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def relative(path):
    return path.relative_to(REPO).as_posix()


def json_bytes(data):
    return (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def csv_bytes(rows):
    """Serialize nested diagnostics as JSON cells, with empty values for NA."""
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v
                         for k, v in row.items()})
    return output.getvalue().encode()


def parent_record(phase):
    """Check pinned Stage 01 bytes, audited metrics, and exact Stage 03 geometry handoff."""
    expected = ACCEPTED[phase]
    source = WP1 / "calculation/01_geometry_optimisation" / phase / "CONTCAR"
    downstream = WP1 / "calculation/03_static_scf" / phase / "POSCAR"
    data = source.read_bytes()
    if sha256(data) != expected["sha256"]:
        raise ValueError(f"{phase}: accepted Stage 01 bytes changed")
    parent = Poscar.from_file(source, check_for_potcar=False).structure
    reference = Poscar.from_file(downstream, check_for_potcar=False).structure
    parameters = (expected["a"], expected["a"], expected["c"], 90, 90, expected["gamma"])
    if (parent.species != reference.species
            or not np.allclose(parent.lattice.matrix, reference.lattice.matrix, atol=1e-8, rtol=0)
            or not np.allclose(parent.frac_coords, reference.frac_coords, atol=1e-8, rtol=0)
            or not np.allclose(parent.lattice.parameters, parameters, atol=2e-6, rtol=0)
            or composition_record(parent)["n_formula_units_if_integer"] != expected["units"]):
        raise ValueError(f"{phase}: audited identity/downstream geometry mismatch")
    criteria = coordination_reference(parent, phase)
    topology = connectivity(parent, bond_graph(parent, criteria))
    expected_rank = 3 if phase == "spinel" else 2
    expected_components = 1 if phase == "spinel" else expected["units"]
    if len(topology) != expected_components or any(c["periodic_rank"] != expected_rank for c in topology):
        raise ValueError(f"{phase}: parent connectivity disagrees with WP1")
    record = dict(phase=phase, source_wp1_path=relative(source), source_sha256=sha256(data),
                  destination_sha256=sha256(data), source_role="accepted relaxed WP1 parent",
                  **composition_record(parent), lattice_vectors_A=parent.lattice.matrix.tolist(),
                  lattice_parameters_A_deg=list(parent.lattice.parameters),
                  accepted_space_group=expected["group"], wp1_acceptance_audit=relative(AUDIT),
                  downstream_geometry_check=relative(downstream), downstream_match=True,
                  downstream_audit=relative(WP1 / "results/03_static_scf/STATIC_SCF_ANALYSIS.md"),
                  coordination_criterion=criteria, connectivity=topology)
    parent.add_site_property("parent_index", list(range(len(parent))))
    return parent, data, record


def construction_settings(basal):
    """One OUC for basal stacking; >=12 A repeat span for other seed families."""
    return dict(min_slab_size=0.1 if basal else 12.0, min_vacuum_size=VACUUM,
                primitive=False, lll_reduce=False, center_slab=False,
                reorient_lattice=False, max_normal_search=2)


def surface_row(slab_id, phase, facet, number, shift, diagnostics, transform, cut, parent_metadata):
    """Keep one complete regeneration/diagnostic record for every unfiltered cut."""
    d = diagnostics
    row = dict(slab_id=slab_id, phase=phase, h=int(facet[0]), k=int(facet[1]), l=int(facet[2]),
               facet_label=f"({facet})", termination_id=f"t{number:02d}", raw_cut_shift=float(shift),
               parent_structure_path=relative(WP2 / "structure/bulk_parents" / phase / "POSCAR"),
               parent_sha256=parent_metadata["source_sha256"], transformation=transform)
    for key in ("composition", "atom_count", "n_formula_units_if_integer", "stoichiometric", "nominal_charge",
                "cell_normal_length_A", "slab_thickness_A", "vacuum_A", "surface_area_A2"):
        row[key] = d[key]
    row.update(upper_face_id=slab_id + "_upper", lower_face_id=slab_id + "_lower",
               faces_equivalent=d["faces_equivalent"])
    for side in ("upper", "lower"):
        for element in ("Zn", "In", "S"):
            row[f"undercoord_{element}_{side}"] = d["faces"][side][f"undercoord_{element}"]
        row[f"{side}_face_composition"] = d["faces"][side]["face_composition"]
        row[f"{side}_termination_fingerprint"] = d["faces"][side]["surface_species_fingerprint"]
    row.update(broken_ZnS_estimate=d["coordination"]["broken_bonds"]["Zn"],
               broken_InS_estimate=d["coordination"]["broken_bonds"]["In"],
               minimum_distance_A=d["minimum_distance_A"], connectivity_status=d["connectivity_status"],
               connectivity_components=d["components"], polarity_screen=d["polarity_screen"],
               formal_dipole_eA=d["formal_dipole_eA"],
               formal_dipole_per_area_e_per_A=d["formal_dipole_per_area_e_per_A"],
               charged_plane_sequence=d["charged_plane_sequence"], cut_diagnostics=cut,
               complete_layer_count=d["complete_layer_count"], duplicate_of="", materialized_path="",
               status="", reason="", notes="Unrelaxed seed; vacuum/thickness not converged; no energy ranking.")
    return row


def build_candidates():
    """Construct and validate all raw seeds in memory before any persistent write."""
    files, rows, faces, parents = {}, [], [], {}
    software = {name: version(name) for name in ("pymatgen", "numpy", "scipy", "spglib")}
    scripts = {p.name: sha256(p.read_bytes()) for p in
               (Path(__file__), Path(__file__).with_name("surface_model_utils.py"))}
    for phase, facets in FACETS.items():
        parent, parent_bytes, parent_metadata = parent_record(phase)
        parents[phase] = parent_metadata
        parent_dir = Path("structure/bulk_parents") / phase
        files[parent_dir / "POSCAR"] = parent_bytes
        files[parent_dir / "metadata.json"] = json_bytes(parent_metadata)
        criteria = parent_metadata["coordination_criterion"]
        for facet in facets:
            basal = facet == "001" and phase != "spinel"
            settings = construction_settings(basal)
            generator = SlabGenerator(parent, tuple(map(int, facet)), **settings)
            # Direct shifts preserve even invalid/duplicate cuts in the raw registry.
            shifts = sorted(float(x % 1) for x in generator.gen_possible_terminations(ftol=PLANE_TOL))
            canonical = []
            for number, shift in enumerate(shifts, start=1):
                slab_id = f"{phase}_{facet}_t{number:02d}"
                raw = generator.get_slab(shift=shift)
                slab, transform = present_slab(raw, generator)
                normal = np.cross(*generator.oriented_unit_cell.lattice.matrix[:2])
                normal /= np.linalg.norm(normal)
                reciprocal = np.array(tuple(map(int, facet))) @ parent.lattice.reciprocal_lattice.matrix
                alignment = float(np.dot(normal, reciprocal / np.linalg.norm(reciprocal)))
                if abs(abs(alignment) - 1) > 1e-8:
                    raise ValueError("Oriented normal does not match requested Miller family")
                transform["normal_in_parent_cartesian"] = normal.tolist()
                transform["normal_sign_relative_to_hkl"] = int(round(alignment))
                cut = cut_diagnostics(generator, shift, criteria)
                diagnostics = diagnose_slab(slab, parent, criteria, basal)
                if basal and not cut["no_first_shell_bond_crosses"]:
                    diagnostics["validation_errors"].append("basal cut crosses bulk first-shell bonds")
                row = surface_row(slab_id, phase, facet, number, shift, diagnostics, transform, cut, parent_metadata)
                errors = diagnostics["validation_errors"]
                if errors:
                    row.update(status="REJECT", reason="; ".join(errors))
                elif phase == "spinel" and diagnostics["polarity_screen"] == "HIGH_CONCERN":
                    row.update(status="DEFERRED", reason="Uncompensated formal dipole exceeds bounded spinel gate; reconstruction/electronic compensation unresolved")
                elif diagnostics["polarity_screen"] == "HIGH_CONCERN":
                    row.update(status="HOLD", reason="Structurally valid seed; high formal-electrostatics concern requires human inspection")
                else:
                    row.update(status="PASS", reason="Structural seed checks pass; electronic stability and convergence untested")
                for previous_slab, previous_diagnostics, previous_id in canonical:
                    if slabs_match(slab, previous_slab, diagnostics, previous_diagnostics):
                        row["duplicate_of"] = previous_id
                        # Preserve independent rejection/defer reasons; suppress valid duplicates with HOLD.
                        if row["status"] == "PASS":
                            row["status"] = "HOLD"
                        row["reason"] = f"Equivalent slab-pair duplicate of {previous_id}; " + row["reason"]
                        break
                if not row["duplicate_of"]:
                    canonical.append((slab, diagnostics, slab_id))
                if row["status"] in ("PASS", "HOLD") and not row["duplicate_of"]:
                    directory = Path("structure/slab_candidates") / phase / facet / slab_id
                    row["materialized_path"] = relative(WP2 / directory / "POSCAR")
                    serialized = Poscar(slab, comment=slab_id, sort_structure=False).get_str(significant_figures=16)
                    parsed = Poscar.from_str(serialized).structure
                    if (parsed.species != slab.species
                            or not np.allclose(parsed.cart_coords, slab.cart_coords, atol=1e-12, rtol=0)
                            or Poscar(parsed, comment=slab_id).get_str(significant_figures=16) != serialized):
                        raise ValueError(f"{slab_id}: unstable POSCAR serialization")
                    files[directory / "POSCAR"] = serialized.encode()
                    metadata = dict(row, hkl=list(map(int, facet)), source_parent=parent_metadata["source_wp1_path"],
                                    construction_settings=settings, plane_tolerance_A=PLANE_TOL,
                                    actual_vacuum_target_A=VACUUM, face_depth_A=FACE_DEPTH,
                                    structure_sha256=sha256(serialized.encode()), diagnostics=diagnostics,
                                    software=software, generation_version="wp2-stage0-v1",
                                    generation_script_hashes=scripts,
                                    generation_script=relative(Path(__file__)),
                                    stacking_extension=(dict(layers_per_oriented_repeat=ACCEPTED[phase]["units"],
                                                             rule="Repeat the same OUC and shift; preserve full stacking repeat",
                                                             later_complete_layer_targets=[2, 4, 6], generated_now=False)
                                                        if basal else None))
                    files[directory / "metadata.json"] = json_bytes(metadata)
                    for side in ("upper", "lower"):
                        face = diagnostics["faces"][side]
                        faces.append(dict(face_id=slab_id + "_" + side, slab_id=slab_id, phase=phase,
                                          facet_label=f"({facet})", termination_id=f"t{number:02d}", side=side,
                                          surface_area_A2=diagnostics["surface_area_A2"], **face,
                                          polarity_screen=diagnostics["polarity_screen"],
                                          equivalent_to_face_id=(slab_id + "_upper" if side == "lower" and
                                                                 diagnostics["faces_equivalent"] else ""),
                                          reaction_face_candidate="PROVISIONAL" if row["status"] == "PASS" else "HOLD",
                                          notes="Distinct side identity retained; reaction suitability untested. Face charge uses outer 3.6 A."))
                rows.append(row)
    results = Path("results/00_surface_registry")
    files[results / "surface_registry.csv"] = csv_bytes(rows)
    files[results / "face_registry.csv"] = csv_bytes(faces)
    return files, rows, faces, parents


def write_outputs(files, output_root):
    """Preflight every output, refusing changed existing bytes before writing any."""
    for path, data in files.items():
        target = output_root / path
        if target.exists() and target.read_bytes() != data:
            raise FileExistsError(f"Existing output differs: {target}; compare explicitly before replacing")
    for phase, facets in FACETS.items():
        for facet in facets:
            (output_root / "structure/slab_candidates" / phase / facet).mkdir(parents=True, exist_ok=True)
    for path, data in files.items():
        target = output_root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            target.write_bytes(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=WP2,
                        help="WP2 output root; use a temporary directory for regeneration checks")
    args = parser.parse_args()
    files, rows, faces, _ = build_candidates()
    write_outputs(files, args.output_root)
    for phase, facets in FACETS.items():
        for facet in facets:
            group = [r for r in rows if r["phase"] == phase and r["facet_label"] == f"({facet})"]
            print(phase, facet, "raw", len(group), "unique", sum(not r["duplicate_of"] for r in group),
                  "materialized", sum(bool(r["materialized_path"]) for r in group),
                  dict(Counter(r["status"] for r in group)))
    print("TOTAL", len(rows), "raw;", len(faces), "face rows;", len(files), "files")


if __name__ == "__main__":
    main()
