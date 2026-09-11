"""Read-only preparation of accepted WP1 PBE spectra and folded band paths.

The accepted numerical analysis supplies formula counts, eigenvalue-derived
edges and path definitions. This module never runs its file-writing entrypoint.
Both plot families share the selected dense-mesh reference: the audited VBM
for publication, or DOSCAR E_F for diagnostics; no line-path Fermi zero is used.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Sequence

import numpy as np
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from cmw.analysis.electronic import aggregate_projections, finite_array, prepare_dos, prepare_path, shift_energy
import analyse_wp1_electronic_structure as accepted

WP1 = Path("vasp/wp1_polymorph_polytype_benchmark")
PHASES = accepted.PHASE_ORDER
LABELS = {"spinel": "Spinel", "alpha1": "α₁", "beta": "β", "IIa_prime": "IIa′", "IIb": "IIb"}
EXPECTED_FU = dict(zip(PHASES, (8, 3, 1, 2, 2)))
SELECTED = ("S-p", "In-s", "S-s", "Zn-d", "Zn-s")
GROUPS = {group: tuple(orbital.name for orbital in Orbital if orbital.orbital_type.name == group)
          for group in ("s", "p", "d")}
METHOD = "PBE electronic structures on PBE+D3(BJ)-relaxed geometries"


def project_root(override: str | Path | None = None) -> Path:
    """Find the explicit WP1 provenance marker, from cwd or this module's location."""
    candidates = [Path(override).expanduser().resolve()] if override else [*Path.cwd().resolve().parents, Path.cwd().resolve(), *Path(__file__).resolve().parents]
    for root in candidates:
        if (root / WP1 / "structure/RAW_STRUCTURE_PROVENANCE.md").is_file():
            return root
    raise ValueError("Cannot find WP1 repository; supply --root or an explicit project_root override")


def source_mapping(root: Path, phases: Sequence[str] = PHASES) -> dict:
    """Resolve one deterministic, calculation-first source per phase and branch."""
    if not phases or len(set(phases)) != len(phases) or any(p not in PHASES for p in phases):
        raise ValueError("Select distinct known WP1 phase identifiers")
    calculation = root / WP1 / "calculation"
    result = {}
    for phase in phases:
        sources = {"dos": calculation / "04_electronic_structure/01_dos_pdos" / phase,
                   "band": calculation / "04_electronic_structure/02_band" / phase,
                   "parent": calculation / "03_static_scf" / phase}
        for branch in ("dos_pdos", "band", "beta", "alpha1"):
            if (calculation / "04_electronic_structure" / phase / branch).exists():
                raise ValueError(f"Ambiguous legacy source exists for {phase}/{branch}")
        for path in sources.values():
            if not path.is_dir():
                raise FileNotFoundError(path)
        result[phase] = sources
    return result


def table_rows(path: Path) -> list[dict]:
    """Read an accepted table without changing it."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def close_value(value: float, expected: float, name: str, atol: float = 1.1e-6) -> None:
    """Fail on a raw/accepted mismatch beyond printed numerical precision."""
    if not np.isfinite(value) or abs(value - expected) > atol:
        raise ValueError(f"{name}: raw {value} differs from accepted {expected} (atol={atol})")


def validate_branch(directory: Path, structure, branch: str):
    """Validate completed fixed-density PBE output, then read existing eigenvalues."""
    incar = accepted.read_incar(directory / "INCAR")
    expected = {"ISPIN": 1, "ICHARG": 11, "NSW": 0, "IBRION": -1, "ENCUT": 500, "LORBIT": 11}
    for key, value in expected.items():
        close_value(float(incar[key]), value, f"{directory}/{key}")
    if incar.get("GGA") != "PE":
        raise ValueError(f"Not PBE: {directory}")
    if not accepted.structures_equivalent(structure, accepted.parse_poscar(directory / "CONTCAR"), tolerance=1e-8):
        raise ValueError(f"Output geometry changed: {directory}")
    evidence = accepted.scan_outcar(directory / "OUTCAR")
    if not evidence["normal_termination"] or not evidence["ediff_reached"] or evidence["fatal_marker"]:
        raise ValueError(f"Not a completed accepted static branch: {directory}")
    text = (directory / "OUTCAR").read_text(errors="strict")
    for key, value in {**expected, "GGA": "PE", "LHFCALC": "F", "LSORBIT": "F", "LNONCOLLINEAR": "F"}.items():
        matches = re.findall(rf"\b{key}\s*=\s*(\S+)", text)
        if not matches or (matches[-1] != value if isinstance(value, str) else float(matches[-1]) != value):
            raise ValueError(f"Unexpected resolved {key}: {directory}")
    eig = accepted.parse_eigenval(directory / "EIGENVAL")
    finite_array(eig.energies, ndim=2, name="eigenvalues")
    finite_array(eig.occupancies, ndim=2, name="occupancies")
    finite_array(eig.kpoints, ndim=2, name="k points")
    for key in ("nkpoints", "nbands", "nelect"):
        close_value(getattr(eig, key), evidence[key], f"{branch}/{key}")
    if accepted.final_electronic_iteration(directory / "OSZICAR") >= int(incar["NELM"]):
        raise ValueError(f"NELM exhausted: {directory}")
    return eig, evidence, incar


@dataclass
class PhaseData:
    """Validated raw arrays and a single explicit reference shared by plot families."""
    phase: str
    sources: dict
    formula_units: int
    reference_ev: float
    reference_mode: str
    dos_energy: np.ndarray
    curves: dict
    band_energy: np.ndarray
    path: object
    endpoint_labels: tuple
    summary: dict
    dos_coverage: tuple[float, float]
    band_coverage: tuple[float, float]


def load_phase(root: Path, phase: str, reference: str = "fermi") -> PhaseData:
    """Validate mapped raw sources against accepted results and prepare one phase.

    Only the DOS XML's spectra/metadata are requested; projected band eigenvalues
    and POTCAR files are not loaded. Band energies come from the existing reader.
    Charge-density identity relies on the accepted audit, not new CHGCAR reads.
    """
    if reference not in ("fermi", "vbm"):
        raise ValueError("Reference must be fermi or vbm")
    sources = source_mapping(root, (phase,))[phase]
    results = root / WP1 / "results/04_electronic_structure"
    rows = [r for r in table_rows(results / "electronic_structure_summary.csv") if r["structure"] == phase]
    if len(rows) != 1:
        raise ValueError(f"Missing/duplicate accepted summary for {phase}")
    summary = rows[0]
    parent = accepted.parse_poscar(sources["parent"] / "POSCAR")
    fu = accepted.formula_units(parent)
    if fu != EXPECTED_FU[phase]:
        raise ValueError(f"Unexpected composition/cell size for {phase}")
    parent_incar = accepted.read_incar(sources["parent"] / "INCAR")
    if parent_incar.get("GGA") != "PE" or int(parent_incar["IVDW"]) != 12:
        raise ValueError("Stage 03 method lineage differs from accepted PBE+D3(BJ)")
    parent_out = accepted.scan_outcar(sources["parent"] / "OUTCAR")
    if not parent_out["normal_termination"] or not parent_out["ediff_reached"]:
        raise ValueError("Stage 03 parent is incomplete")
    eigs, evidence = {}, {}
    for branch in ("dos", "band"):
        structure = accepted.parse_poscar(sources[branch] / "POSCAR")
        if not accepted.structures_equivalent(parent, structure):
            raise ValueError(f"Stage 03/04 structure or site-order mismatch: {phase}/{branch}")
        eigs[branch], evidence[branch], _ = validate_branch(sources[branch], structure, branch)
    mesh, shift = accepted.parse_uniform_kpoints(sources["dos"] / "KPOINTS")
    if mesh != accepted.PHASE_META[phase]["mesh"] or shift != (0., 0., 0.):
        raise ValueError("DOS branch is not the accepted dense Gamma mesh")
    for branch, column in (("dos", "pbe_fundamental_gap_eV"), ("band", "band_path_gap_eV")):
        edges = accepted.band_edges(eigs[branch])
        close_value(edges["fundamental_gap"], float(summary[column]), f"{phase}/{column}")
        if branch == "dos":
            for field, column_name in (("vbm", "uniform_VBM_energy_eV"), ("cbm", "uniform_CBM_energy_eV")):
                close_value(edges[field], float(summary[column_name]), f"{phase}/{field}")
    path_definition = accepted.parse_path_kpoints(sources["band"] / "KPOINTS")
    labels = tuple((a, b) for _, a, _, b in path_definition.segments)
    if path_definition.points_per_segment != 40 or labels != accepted.PHASE_META[phase]["segments"]:
        raise ValueError("Actual labelled path differs from accepted topology/sample count")
    endpoints = [(a, b) for a, _, b, _ in path_definition.segments]
    path = prepare_path(parent.lattice, endpoints, eigs["band"].kpoints, path_definition.points_per_segment)
    # Reuse the accepted segmentation/distance check on intended (possibly unwrapped) samples.
    intended_eig = accepted.EigenvalData(eigs['band'].nelect, eigs['band'].nkpoints, eigs['band'].nbands,
        path.intended_fractional.tolist(), eigs['band'].weights, eigs['band'].energies, eigs['band'].occupancies)
    descriptors = accepted.validate_path_sampling(path_definition, intended_eig, accepted.reciprocal_lattice(parent.lattice))
    np.testing.assert_allclose(path.distances, [r['distance'] for r in descriptors], atol=1e-10, rtol=0)

    run = Vasprun(sources['dos'] / 'vasprun.xml', parse_dos=True, parse_eigen=False,
                  parse_projected_eigen=False, parse_potcar_file=False, exception_on_bad_xml=True)
    if run.parameters['ISPIN'] != 1 or run.parameters['ICHARG'] != 11 or run.dos_has_errors:
        raise ValueError("Unexpected spin, fixed-density mode or malformed DOS XML")
    complete = run.complete_dos
    xml_structure = complete.structure
    if [s.specie.symbol for s in xml_structure] != parent.site_elements:
        raise ValueError("DOS XML site order differs from POSCAR")
    np.testing.assert_allclose(xml_structure.lattice.matrix, parent.lattice, atol=1e-7, rtol=0)
    delta = xml_structure.frac_coords - parent.frac_coords
    np.testing.assert_allclose(delta - np.rint(delta), 0, atol=1e-7, rtol=0)
    if set(complete.densities) != {Spin.up}:
        raise ValueError("Expected the single physical ISPIN=1 DOS channel")
    site_channels = []
    for site in xml_structure:
        orbitals = complete.pdos[site]
        if any(set(densities) != {Spin.up} for densities in orbitals.values()):
            raise ValueError("Unexpected projected spin channels")
        site_channels.append({orbital.name: densities[Spin.up] for orbital, densities in orbitals.items()})
    curves = aggregate_projections(parent.site_elements, site_channels, GROUPS)
    curves['TDOS'] = np.asarray(complete.densities[Spin.up])
    for element in accepted.ELEMENT_ORDER:
        curves[element] = sum(curves[f'{element}-{orbital}'] for orbital in ('s', 'p', 'd'))
        for orbital, projection in complete.get_element_spd_dos(element).items():
            np.testing.assert_allclose(curves[f'{element}-{orbital.name}'], projection.densities[Spin.up], atol=1e-10, rtol=1e-10)
    raw_dos = accepted.parse_doscar(sources['dos'] / 'DOSCAR', parent.site_elements)
    close_value(raw_dos.efermi, float(summary['dos_efermi_eV']), f'{phase}/reference')
    close_value(raw_dos.efermi, complete.efermi, f'{phase}/XML reference', atol=1e-3)
    close_value(raw_dos.efermi, evidence['dos']['efermi'], f'{phase}/OUTCAR reference', atol=1e-3)
    np.testing.assert_allclose(complete.energies, raw_dos.energies, atol=1.1e-3, rtol=0)
    # DOSCAR uses four significant figures for densities; XML retains more.
    np.testing.assert_allclose(curves['TDOS'], raw_dos.total_dos, atol=1e-3, rtol=5.1e-4)
    for (element, orbital), density in raw_dos.element_orbitals.items():
        # XML site components have four decimal places: propagate their rounding
        # bounds through the sum, rather than using a cell-size-independent floor.
        rounding_bound = parent.site_elements.count(element) * len(GROUPS[orbital]) * 0.00005 + 1e-6
        np.testing.assert_allclose(curves[f'{element}-{orbital}'], density, atol=rounding_bound, rtol=5.1e-4)
    zero = raw_dos.efermi if reference == 'fermi' else float(accepted.band_edges(eigs['dos'])['vbm'])
    energy, validated_curves = prepare_dos(complete.energies, curves, reference_ev=zero)
    bands = shift_energy(eigs['band'].energies, zero)
    primitive = SpacegroupAnalyzer(xml_structure, symprec=1e-5).find_primitive()
    multiplicity = len(xml_structure) // len(primitive)
    if multiplicity != {'spinel': 4, 'alpha1': 3}.get(phase, 1):
        raise ValueError(f"Unexpected primitive-cell multiplicity: {phase}/{multiplicity}")
    # Conservative upper coverage: every sampled k has computed bands up to this energy.
    dos_eigen = np.asarray(eigs['dos'].energies) - zero
    dos_coverage = (max(energy[0], np.max(dos_eigen[:, 0])), min(energy[-1], np.min(dos_eigen[:, -1])))
    band_coverage = (float(np.max(bands[:, 0])), float(np.min(bands[:, -1])))
    validation = {'phase': LABELS[phase], 'formula_units': fu, 'spin': 'ISPIN=1 (single channel)',
        'reference_mode': reference, 'reference_eV': zero, 'reference_source': str((sources['dos'] / ('DOSCAR' if reference == 'fermi' else 'EIGENVAL')).relative_to(root)),
        'dos_kpoints': eigs['dos'].nkpoints, 'band_kpoints': eigs['band'].nkpoints,
        'bands': eigs['band'].nbands, 'points_per_segment': 40, 'primitive_cells': multiplicity,
        'uniform_gap_eV': float(summary['pbe_fundamental_gap_eV']), 'path_gap_eV': float(summary['band_path_gap_eV']),
        'classification': summary['sampled_gap_classification'], 'validation': 'PASS'}
    return PhaseData(phase, sources, fu, zero, reference, energy, validated_curves, bands, path, labels, validation, dos_coverage, band_coverage)


def load_all(root: Path, phases: Sequence[str] = PHASES, reference: str = 'fermi') -> list[PhaseData]:
    """Load phases sequentially; retain compact arrays rather than parser objects."""
    source_mapping(root, phases)
    return [load_phase(root, phase, reference) for phase in phases]


def energy_window(data: Sequence[PhaseData], family: str, requested=None) -> tuple[float, float]:
    """Use a common window lying inside all calculated energy/band coverage."""
    coverage = [d.dos_coverage if family == 'dos' else d.band_coverage for d in data]
    lower, upper = max(c[0] for c in coverage), min(c[1] for c in coverage)
    window = tuple(requested) if requested is not None else (max(-8., float(np.ceil(lower))), min(4., float(np.floor(upper))))
    if len(window) != 2 or not np.isfinite(window).all() or not lower <= window[0] < window[1] <= upper:
        raise ValueError(f"Requested window {window} exceeds shared calculated coverage {(lower, upper)}")
    return window


def prepared_dos(data: PhaseData, normalization='formula_unit'):
    """Scale all already aligned raw DOS curves equally; never shift a second time."""
    if normalization not in ('formula_unit', 'cell'):
        raise ValueError('Normalization must be formula_unit or cell')
    return prepare_dos(data.dos_energy, data.curves, reference_ev=0,
                       divisor=data.formula_units if normalization == 'formula_unit' else 1)
