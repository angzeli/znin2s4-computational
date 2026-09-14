"""Archive the completed beta HSE06 optional path and plot its portable CSV.

Default: validate native XML, export all 180 x 40 eigenvalues, then plot the CSV.
With --csv: read only the portable CSV; no native VASP files or JSON are needed.
Use the existing .venv-wp1-plots environment and its read-only CMW analysis
dependency (or put that checkout's src directory on PYTHONPATH).

Indices are one-based; fractional k coordinates and raw energies retain their
printed XML precision. One global band-31 path VBM defines the energy zero.
Physical distances use the accepted reciprocal basis, including 2*pi, and carry
forward across branches without adding a disconnected jump. Display-only gaps
come from the frozen PBE plotter. No smoothing or extra eigenvalues are added.

Bands 38-40 are retained unchanged but are not uniformly quantitatively converged:
repeated-point discrepancies reach 0.1421 eV. They lie above +5.4992 eV relative
to the VBM, outside the canonical -4 to +4 eV figure. The band-31/32 sampled gap
is accepted; neither continuous-zone extrema nor an optical gap is established.
Import and --help do not analyse data or write outputs. No calculations are run.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import analyse_wp1_electronic_structure as accepted
from plot_wp1_bands import BRANCH_GAPS, display_distances, plot_bands
from wp1_electronic_plot_data import WP1, prepare_path
from wp1_electronic_plot_style import RC

ROOT = Path(__file__).resolve().parents[3]
SOURCE = WP1 / "calculation/05_hybrid_validation/02_hse06_band/beta"
OUTPUT = WP1 / "results/05_hybrid_validation"
CSV_NAME = "beta_hse06_kpoints_opt_eigenvalues.csv"
JSON_NAME = "beta_hse06_kpoints_opt_metadata.json"
PDF_NAME = "beta_hse06_band_structure.pdf"
LABELS = accepted.PHASE_META["beta"]["segments"]
POINTS_PER_SEGMENT = 20
NBANDS = 40
WINDOW = (-4., 4.)
COLUMNS = (
    "path_point_index", "branch_index", "segment_index", "segment_start_label",
    "segment_end_label", "point_index_in_segment", "is_segment_start",
    "is_segment_end", "is_repeated_boundary_point", "kx_frac", "ky_frac", "kz_frac",
    "k_distance_Ainv", "display_x", "band_index", "eigenvalue_eV", "energy_minus_vbm_eV",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    """Hash only the small, explicitly archived source/code identities."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def band_edges(raw):
    """Return sampled edges in eV for nonmagnetic NELECT=62, bands 31/32."""
    require(raw.shape == (180, NBANDS) and np.isfinite(raw).all(),
            "Expected exactly 180 x 40 finite optional-path eigenvalues")
    require(np.all(np.diff(raw, axis=1) >= 0), "Eigenvalues are not band-ordered")
    vbm, cbm = float(raw[:, 30].max()), float(raw[:, 31].min())
    gap, direct = cbm - vbm, float((raw[:, 31] - raw[:, 30]).min())
    require(np.allclose([gap, direct], [1.0907, 1.0919], atol=1e-8, rtol=0),
            "STOP: optional-path gaps differ from the accepted four-decimal XML result")
    require(float(raw[:, 37:].min()) - vbm > WINDOW[1],
            "STOP: a problematic band 38-40 enters the +4 eV figure window")
    return {"vbm_band_index": 31, "cbm_band_index": 32,
            "sampled_vbm_raw_eV": vbm, "sampled_cbm_raw_eV": cbm,
            "fundamental_gap_eV": round(gap, 4), "minimum_direct_gap_eV": round(direct, 4),
            "direct_minus_fundamental_meV": round((direct - gap) * 1000, 1)}


def repeated_points(kpoints, raw):
    """Compare every periodic-equivalent point pair; never merge its energies."""
    delta = kpoints[:, None, :] - kpoints[None, :, :]
    same = np.max(np.abs(delta - np.rint(delta)), axis=2) <= 1e-8
    pairs = np.argwhere(np.triu(same, k=1))
    require(len(pairs) == 18, "Expected the accepted 18 repeated-point pairs")
    differences = np.abs(raw[pairs[:, 0]] - raw[pairs[:, 1]])
    maxima = differences.max(axis=0)
    require(np.all(maxima[:37] == 0), "Repeated points disagree within bands 1-37")
    return pairs, {
        "pair_count": len(pairs), "path_point_pairs_one_based": (pairs + 1).tolist(),
        "maximum_difference_bands_1_to_37_eV": float(maxima[:37].max()),
        "maximum_difference_by_band_eV": {
            str(b): round(float(maxima[b - 1]), 4) for b in (31, 32, 38, 39, 40)},
        "minimum_energy_minus_vbm_bands_38_to_40_eV":
            round(float(raw[:, 37:].min() - raw[:, 30].max()), 4),
    }


def dataset(kpoints, raw, path, labels):
    """Adapt validated arrays to the unchanged Stage 04 plot_bands interface."""
    edges = band_edges(raw)
    require(kpoints.shape == (180, 3) and np.isfinite(kpoints).all(), "Invalid path coordinates")
    require(tuple(labels) == LABELS and path.disconnected_before == (7, 8),
            "Unexpected hP2 labels or disconnected branches")
    repeated_points(kpoints, raw)
    return SimpleNamespace(phase="beta", reference_mode="vbm",
        reference_ev=edges["sampled_vbm_raw_eV"], raw_energy=raw, kpoints=kpoints,
        band_energy=raw - edges["sampled_vbm_raw_eV"], path=path,
        endpoint_labels=tuple(labels), summary={"primitive_cells": 1})


def read_source(root=ROOT):
    """Validate the specific completed XML and its accepted geometry/path lineage.

    Only eigenvalues_kpoints_opt supplies path energies. The regular XML array
    is checked separately for the occupied/empty convention, not used as a path.
    Eight-decimal XML coordinates allow 5.1e-9 rounding error per component;
    four-decimal eigenvalues must reproduce the already accepted rounded gaps.
    """
    root = Path(root).resolve()
    source = root / SOURCE
    tree = ET.parse(source / "vasprun.xml").getroot()
    blocks = tree.findall("calculation/eigenvalues_kpoints_opt")
    require(tree.tag == "modeling" and len(blocks) == 1, "Missing/ambiguous optional eigenvalue block")
    params = {e.attrib.get("name", "").strip(): e.text.strip()
              for e in tree.findall("parameters//i")}
    expected = {"NELECT": 62, "ISPIN": 1, "NBANDS": 40, "ENMAX": 500,
                "AEXX": .25, "HFSCREEN": .2}
    require(all(float(params[k]) == v for k, v in expected.items()), "Unexpected resolved HSE settings")
    require(params["LHFCALC"] == "T" and params["LSORBIT"] == "F", "Unexpected hybrid/spin mode")
    require(tree.findtext("generator/i[@name='version']").strip() == "6.6.1", "Unexpected VASP version")
    spins = blocks[0].findall("eigenvalues/array/set/set")
    require(len(spins) == 1 and spins[0].get("comment") == "spin 1", "Expected one spin channel")
    point_sets = spins[0].findall("set")
    require([s.get("comment") for s in point_sets] == [f"kpoint {i}" for i in range(1, 181)],
            "Optional path count/order is incomplete")
    tokens = [[r.text.split() for r in s.findall("r")] for s in point_sets]
    require(all(len(rows) == NBANDS and all(len(r) == 1 for r in rows) for rows in tokens),
            "Expected 40 scalar eigenvalues per optional point")
    require(all(len(r[0].split(".")[-1]) == 4 for rows in tokens for r in rows),
            "Optional eigenvalue printed precision changed")
    raw = np.array([[float(r[0]) for r in rows] for rows in tokens])
    kpoints = np.array([list(map(float, v.text.split())) for v in
                       blocks[0].findall("kpoints/varray[@name='kpointlist']/v")])
    regular_spins = tree.findall("calculation/eigenvalues/array/set/set")
    require(len(regular_spins) == 1, "Expected one regular-mesh spin channel")
    regular = np.array([[list(map(float, r.text.split())) for r in s.findall("r")]
                        for s in regular_spins[0].findall("set")])
    require(regular.shape == (69, 40, 2) and np.isfinite(regular).all(), "Invalid regular XML array")
    require(np.all(regular[:, :31, 1] == 1) and np.all(regular[:, 31:, 1] == 0),
            "Regular occupations do not confirm bands 31/32 as the edge pair")
    mesh, shift = accepted.parse_uniform_kpoints(source / "KPOINTS")
    require(mesh == (12, 12, 4) and shift == (0., 0., 0.), "Unexpected regular Gamma mesh")

    structure = accepted.parse_poscar(source / "POSCAR")
    require(dict(zip(structure.species, structure.counts)) == {"S": 4, "In": 2, "Zn": 1},
            "Unexpected beta composition/cell")
    geometry_references = [WP1 / "calculation/01_geometry_optimisation/beta/CONTCAR",
                           WP1 / "calculation/04_electronic_structure/02_band/beta/POSCAR",
                           SOURCE / "CONTCAR"]
    for reference in geometry_references:
        require(accepted.structures_equivalent(structure, accepted.parse_poscar(root / reference)),
                f"Accepted geometry mismatch: {reference}")
    elements = [r.findtext("c").strip() for r in tree.findall("atominfo/array[@name='atoms']/set/rc")]
    require(elements == structure.site_elements, "XML atom order differs from POSCAR")
    for name in ("initialpos", "finalpos"):
        node = tree.find(f"structure[@name='{name}']")
        require(node is not None, f"Missing XML {name}")
        lattice = np.array([list(map(float, v.text.split())) for v in node.findall("crystal/varray[@name='basis']/v")])
        coords = np.array([list(map(float, v.text.split())) for v in node.findall("varray[@name='positions']/v")])
        delta = coords - np.asarray(structure.frac_coords)
        require(np.allclose(lattice, structure.lattice, atol=5.1e-9, rtol=0)
                and np.max(np.abs(delta - np.rint(delta))) <= 5.1e-9, "XML geometry differs from POSCAR")

    definition = accepted.parse_path_kpoints(source / "KPOINTS_OPT")
    require(definition.points_per_segment == POINTS_PER_SEGMENT, "Expected 20 points per HSE segment")
    labels = tuple((a, b) for _, a, _, b in definition.segments)
    endpoints = [(a, b) for a, _, b, _ in definition.segments]
    path = prepare_path(structure.lattice, endpoints, kpoints, definition.points_per_segment)
    require(np.max(np.abs(kpoints - path.intended_fractional)) <= 5.1e-9,
            "Optional coordinates differ from requested path at XML precision")
    data = dataset(kpoints, raw, path, labels)
    _, repeats = repeated_points(kpoints, raw)
    identities = {name: {"path": str(SOURCE / name), "sha256": sha256(source / name)}
                  for name in ("vasprun.xml", "KPOINTS_OPT", "POSCAR")}
    metadata = {
        "phase": "beta", "system": "ZnIn2S4", "method": "HSE06", "vasp_version": "6.6.1",
        "source_files": identities, "xml_dataset": "modeling/calculation/eigenvalues_kpoints_opt",
        "geometry": {"atom_count": 7, "composition": "ZnIn2S4",
            "lattice_vectors_angstrom": structure.lattice,
            "checked_against": [str(p) for p in geometry_references],
            "poscar_tolerance": 1e-10, "xml_rounding_tolerance": 5.1e-9, "transformed": False},
        "settings": {"ENCUT_eV": 500, "AEXX": .25, "HFSCREEN_Ainv": .2,
                     "ISPIN": 1, "NELECT": 62, "NBANDS": 40},
        "regular_sampling": {"mesh": mesh, "centering": "Gamma", "shift": shift, "NKPTS": 69},
        "path": {"definition": "GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K",
            "segments": [{"start_label": a, "end_label": b, "start_fractional": u, "end_fractional": v}
                         for u, a, v, b in definition.segments],
            "points_per_segment": POINTS_PER_SEGMENT, "path_entries": 180,
            "exported_eigenvalues": 7200, "disconnected_before_segment_one_based": [8, 9]},
        "edges": band_edges(raw),
        "energy_zero": "One global optional-path band-31 maximum; E_minus_VBM = raw_eigenvalue - VBM; not E_F",
        "printed_precision": {"eigenvalue_eV": .0001, "fractional_k_component": 1e-8},
        "csv_conventions": {"index_base": 1, "row_order": "path point, then band",
            "boolean_encoding": "0/1",
            "repeated_boundary": "All endpoint occurrences in an equivalent pair, including the first occurrence",
            "k_distance_Ainv": "Cumulative segment arc length from POSCAR 2*pi reciprocal basis and KPOINTS_OPT endpoints; no disconnected jump; branch starts inherit the prior end distance",
            "display_x": "Physical distance plus a display-only offset of 0.055 times the total physical span per break"},
        "repeated_point_validation": repeats,
        "high_unoccupied_band_caveat": "All 40 bands are exported unchanged. Bands 38-40 are not uniformly quantitatively converged; repeated-point differences reach 0.1421 eV. Their minimum is +5.4992 eV relative to VBM, outside the -4 to +4 eV figure. The band-31/32 gap is unaffected; no high-energy optical validation is claimed.",
        "figure": {"path": str(OUTPUT / "figures" / PDF_NAME), "figsize_inches": [7, 5],
                   "energy_window_eV": WINDOW, "hP2_display_gap_fraction": BRANCH_GAPS["hP2"]},
        "implementation": {"script": str(WP1 / "script" / Path(__file__).name),
            "sha256": sha256(__file__), "numpy": np.__version__, "matplotlib": mpl.__version__,
            "shared_plotter": str(WP1 / "script/plot_wp1_bands.py"),
            "shared_style": str(WP1 / "script/wp1_electronic_plot_style.py")},
    }
    return data, metadata


def point_records(data):
    """One metadata record per original point; repeated endpoints stay distinct."""
    pairs, _ = repeated_points(data.kpoints, data.raw_energy)
    repeated = set(pairs.ravel().tolist())
    shown = display_distances(data)
    records, branch = [], 1
    for segment_index, (segment, labels) in enumerate(zip(data.path.segments, data.endpoint_labels)):
        if segment_index in data.path.disconnected_before:
            branch += 1
        for local, index in enumerate(range(segment.start, segment.stop), start=1):
            start, end = local == 1, index == segment.stop - 1
            records.append((index + 1, branch, segment_index + 1, *labels, local,
                int(start), int(end), int((start or end) and index in repeated),
                *(f"{v:.8f}" for v in data.kpoints[index]),
                repr(float(data.path.distances[index])), repr(float(shown[index]))))
    return records


def write_csv(path, data):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(COLUMNS)
        for index, record in enumerate(point_records(data)):
            for band in range(NBANDS):
                writer.writerow((*record, band + 1, f"{data.raw_energy[index, band]:.4f}",
                                 f"{data.band_energy[index, band]:.4f}"))


def read_csv(path):
    """Reconstruct the plotting dataset without XML, JSON, POSCAR or EIGENVAL."""
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        require(tuple(next(reader)) == COLUMNS, "Unexpected portable CSV schema")
        rows = list(reader)
    require(len(rows) == 7200 and all(len(r) == len(COLUMNS) for r in rows), "Expected 7200 complete CSV rows")
    require(all(int(r[0]) == i // NBANDS + 1 and int(r[14]) == i % NBANDS + 1
                for i, r in enumerate(rows)), "Missing, duplicate or reordered point/band indices")
    points = rows[::NBANDS]
    require(all(r[:14] == points[i // NBANDS][:14] for i, r in enumerate(rows)),
            "Point metadata differs between bands")
    raw = np.array([float(r[15]) for r in rows]).reshape(180, NBANDS)
    shifted = np.array([float(r[16]) for r in rows]).reshape(180, NBANDS)
    coords = np.array([[float(v) for v in r[9:12]] for r in points])
    distances = np.array([float(r[12]) for r in points])
    require(np.isfinite(distances).all(), "Non-finite physical distance")
    segments = tuple(slice(i, i + POINTS_PER_SEGMENT) for i in range(0, 180, POINTS_PER_SEGMENT))
    labels = tuple(tuple(points[s.start][3:5]) for s in segments)
    breaks = tuple(i for i, s in enumerate(segments) if i and points[s.start][1] != points[s.start - 1][1])
    path_data = SimpleNamespace(distances=distances, segments=segments, disconnected_before=breaks)
    data = dataset(coords, raw, path_data, labels)
    require(np.allclose(shifted, data.band_energy, atol=1e-12, rtol=0), "CSV does not use one common VBM shift")
    # Regenerating only metadata checks branch/segment/endpoint flags and display
    # offsets; energies are never interpolated or repaired during this check.
    expected_points = point_records(data)
    require(all(r[:14] == [str(v) for v in expected] for r, expected in zip(points, expected_points)),
            "CSV path metadata or display spacing is inconsistent")
    require(distances[0] == 0, "Physical path must begin at zero")
    for segment in segments:
        steps = np.diff(distances[segment])
        require(np.all(steps > 0) and np.allclose(steps, steps[0], atol=1e-12, rtol=0),
                "Invalid within-segment reciprocal distances")
        if segment.start:
            require(distances[segment.start] == distances[segment.start - 1],
                    "Physical distance bridges a segment boundary")
    return data


def save_pdf(data, destination):
    """Reuse the frozen plotter unchanged; only the HSE basename/export differs."""
    fig, _, _ = plot_bands(data, window=WINDOW)
    try:
        with mpl.rc_context(RC):
            fig.canvas.draw()
            # Unlike the PBE saver's extra tight-bbox padding, preserve an exact
            # 7 x 5 inch PDF page as requested; constrained layout keeps labels in.
            fig.savefig(destination, format="pdf", facecolor="white", transparent=False,
                        bbox_inches=None, metadata={"CreationDate": None, "ModDate": None})
    finally:
        plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=ROOT, help="Scientific repository root")
    parser.add_argument("--csv", type=Path, help="Plot only this portable CSV; do not read native source files")
    parser.add_argument("--output-dir", type=Path, help="Default: WP1 results/05_hybrid_validation")
    parser.add_argument("--overwrite", action="store_true", help="Replace only the three explicitly named archival outputs")
    args = parser.parse_args(argv)
    output = args.output_dir or args.root / OUTPUT
    csv_path, metadata_path = output / CSV_NAME, output / JSON_NAME
    pdf_path = output / "figures" / PDF_NAME
    targets = [pdf_path] if args.csv else [csv_path, metadata_path, pdf_path]
    if not args.overwrite and any(p.exists() for p in targets):
        raise FileExistsError("Archival output exists; use --overwrite for these named artifacts")
    if args.csv:
        data = read_csv(args.csv)
    else:
        data, metadata = read_source(args.root)
        output.mkdir(parents=True, exist_ok=True)
        write_csv(csv_path, data)
        exported = read_csv(csv_path)
        require(np.array_equal(exported.raw_energy, data.raw_energy), "CSV round trip changed raw energies")
        metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        data = exported
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    save_pdf(data, pdf_path)
    print(json.dumps({"rows": 7200, "path_points": 180, "bands": 40,
                      **band_edges(data.raw_energy), "figure": str(pdf_path)}, indent=2))


if __name__ == "__main__":
    main()
