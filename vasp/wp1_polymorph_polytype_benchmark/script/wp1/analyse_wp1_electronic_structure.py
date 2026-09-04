#!/usr/bin/env python3
"""Validate and summarize WP1 Stage 04 VASP electronic-structure outputs.

The analysis is deliberately dependency-free and numerical only. It reads the
existing fixed-geometry VASP outputs, derives sampled band edges from EIGENVAL,
integrates element/orbital DOSCAR projections near those edges, and writes the
compact Stage 04 Markdown and CSV deliverables. It never invokes VASP or plots.
"""

from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import reduce
from pathlib import Path
from typing import Iterable, Sequence


BENCHMARK_ROOT = Path(__file__).resolve().parents[2]
CALC_ROOT = BENCHMARK_ROOT / "calculation" / "04_electronic_structure"
STAGE03_ROOT = BENCHMARK_ROOT / "calculation" / "03_static_scf"
RESULT_ROOT = BENCHMARK_ROOT / "results" / "04_electronic_structure"

SUMMARY_CSV = RESULT_ROOT / "electronic_structure_summary.csv"
EDGES_CSV = RESULT_ROOT / "band_edges.csv"
PDOS_CSV = RESULT_ROOT / "pdos_band_edge_character.csv"
REPORT_MD = RESULT_ROOT / "ELECTRONIC_STRUCTURE_ANALYSIS.md"

PHASE_ORDER = ("spinel", "alpha1", "beta", "IIa_prime", "IIb")
ELEMENT_ORDER = ("Zn", "In", "S")
ORBITAL_ORDER = ("s", "p", "d")
WINDOWS_EV = (0.20, 0.50)
EDGE_TOL_EV = 1.0e-5
DIRECT_TOL_EV = 0.010
DIRECT_EXACT_TOL_EV = 0.001
COORDINATION_CUTOFF_ANG = 3.0

PHASE_META = {
    "spinel": {
        "display": "Spinel",
        "space_group": "Fd-3m (No. 227)",
        "path_type": "cF2",
        "mesh": (6, 6, 6),
        "segments": (
            ("GAMMA", "X"),
            ("X", "U"),
            ("K", "GAMMA"),
            ("GAMMA", "L"),
            ("L", "W"),
            ("W", "X"),
        ),
        "path": "GAMMA-X-U | K-GAMMA-L-W-X",
    },
    "alpha1": {
        "display": "Alpha1",
        "space_group": "R3m (No. 160)",
        "path_type": "hR1",
        "mesh": (18, 18, 6),
        "segments": (
            ("GAMMA", "T"),
            ("T", "H_2"),
            ("H_0", "L"),
            ("L", "GAMMA"),
            ("GAMMA", "S_0"),
            ("S_2", "F"),
            ("F", "GAMMA"),
        ),
        "path": "GAMMA-T-H_2 | H_0-L-GAMMA-S_0 | S_2-F-GAMMA",
    },
    "beta": {
        "display": "Beta",
        "space_group": "P3m1 (No. 156)",
        "path_type": "hP2",
        "mesh": (18, 18, 6),
        "segments": (
            ("GAMMA", "M"),
            ("M", "K"),
            ("K", "GAMMA"),
            ("GAMMA", "A"),
            ("A", "L"),
            ("L", "H"),
            ("H", "A"),
            ("L", "M"),
            ("H", "K"),
        ),
        "path": "GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K",
    },
    "IIa_prime": {
        "display": "IIa-prime",
        "space_group": "P-3m1 (No. 164)",
        "path_type": "hP2",
        "mesh": (18, 18, 6),
        "segments": (
            ("GAMMA", "M"),
            ("M", "K"),
            ("K", "GAMMA"),
            ("GAMMA", "A"),
            ("A", "L"),
            ("L", "H"),
            ("H", "A"),
            ("L", "M"),
            ("H", "K"),
        ),
        "path": "GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K",
    },
    "IIb": {
        "display": "IIb",
        "space_group": "P6_3mc (No. 186)",
        "path_type": "hP2",
        "mesh": (18, 18, 6),
        "segments": (
            ("GAMMA", "M"),
            ("M", "K"),
            ("K", "GAMMA"),
            ("GAMMA", "A"),
            ("A", "L"),
            ("L", "H"),
            ("H", "A"),
            ("L", "M"),
            ("H", "K"),
        ),
        "path": "GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K",
    },
}

ATOMIC_MASS = {"Zn": 65.38, "In": 114.818, "S": 32.06}


@dataclass
class PoscarData:
    lattice: list[tuple[float, float, float]]
    species: list[str]
    counts: list[int]
    site_elements: list[str]
    frac_coords: list[tuple[float, float, float]]


@dataclass
class EigenvalData:
    nelect: int
    nkpoints: int
    nbands: int
    kpoints: list[tuple[float, float, float]]
    weights: list[float]
    energies: list[list[float]]
    occupancies: list[list[float]]


@dataclass
class DoscarData:
    nions: int
    nedos: int
    efermi: float
    energies: list[float]
    total_dos: list[float]
    integrated_dos: list[float]
    element_orbitals: dict[tuple[str, str], list[float]]
    site_orbitals: list[dict[str, list[float]]]


@dataclass
class PathData:
    points_per_segment: int
    segments: list[
        tuple[
            tuple[float, float, float],
            str,
            tuple[float, float, float],
            str,
        ]
    ]


def fail(message: str) -> None:
    raise RuntimeError(message)


def vector_sub(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vector_add(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vector_scale(value: float, vector: Sequence[float]) -> tuple[float, float, float]:
    return (value * vector[0], value * vector[1], value * vector[2])


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(vector: Sequence[float]) -> float:
    return math.sqrt(dot(vector, vector))


def frac_to_cart(
    frac: Sequence[float], lattice: Sequence[Sequence[float]]
) -> tuple[float, float, float]:
    return (
        sum(frac[i] * lattice[i][0] for i in range(3)),
        sum(frac[i] * lattice[i][1] for i in range(3)),
        sum(frac[i] * lattice[i][2] for i in range(3)),
    )


def reciprocal_lattice(
    lattice: Sequence[Sequence[float]],
) -> list[tuple[float, float, float]]:
    volume = dot(lattice[0], cross(lattice[1], lattice[2]))
    if abs(volume) < 1.0e-12:
        fail("Singular real-space lattice")
    factor = 2.0 * math.pi / volume
    return [
        vector_scale(factor, cross(lattice[1], lattice[2])),
        vector_scale(factor, cross(lattice[2], lattice[0])),
        vector_scale(factor, cross(lattice[0], lattice[1])),
    ]


def periodic_fractional_delta(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    delta = vector_sub(a, b)
    return tuple(value - round(value) for value in delta)  # type: ignore[return-value]


def read_incar(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].split("!", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip().upper()] = value.strip()
    return values


def incar_int(values: dict[str, str], key: str) -> int:
    return int(float(values[key]))


def incar_float(values: dict[str, str], key: str) -> float:
    return float(values[key].replace("D", "E").replace("d", "e"))


def parse_poscar(path: Path) -> PoscarData:
    lines = path.read_text().splitlines()
    if len(lines) < 8:
        fail(f"Incomplete POSCAR: {path}")
    scale_parts = lines[1].split()
    if len(scale_parts) != 1:
        fail(f"Unsupported three-factor POSCAR scale: {path}")
    scale = float(scale_parts[0])
    if scale <= 0:
        fail(f"Unsupported non-positive POSCAR scale: {path}")
    lattice = []
    for line in lines[2:5]:
        values = [float(value) * scale for value in line.split()[:3]]
        if len(values) != 3:
            fail(f"Malformed lattice vector: {path}")
        lattice.append(tuple(values))
    species = lines[5].split()
    try:
        counts = [int(value) for value in lines[6].split()]
    except ValueError as exc:
        raise RuntimeError(f"VASP 5 species labels required in {path}") from exc
    if len(species) != len(counts):
        fail(f"Species/count mismatch: {path}")
    coord_line = 7
    if lines[coord_line].strip().lower().startswith("s"):
        coord_line += 1
    coordinate_mode = lines[coord_line].strip().lower()
    coord_line += 1
    if not coordinate_mode.startswith("d"):
        fail(f"Only Direct coordinates are supported: {path}")
    site_elements = [element for element, count in zip(species, counts) for _ in range(count)]
    coordinates = []
    for line in lines[coord_line : coord_line + len(site_elements)]:
        values = tuple(float(value) for value in line.split()[:3])
        if len(values) != 3:
            fail(f"Malformed coordinate in {path}")
        coordinates.append(values)
    if len(coordinates) != len(site_elements):
        fail(f"Coordinate count mismatch: {path}")
    return PoscarData(lattice, species, counts, site_elements, coordinates)


def structures_equivalent(a: PoscarData, b: PoscarData, tolerance: float = 1.0e-10) -> bool:
    if a.site_elements != b.site_elements:
        return False
    for row_a, row_b in zip(a.lattice, b.lattice):
        if max(abs(x - y) for x, y in zip(row_a, row_b)) > tolerance:
            return False
    for coord_a, coord_b in zip(a.frac_coords, b.frac_coords):
        if max(abs(value) for value in periodic_fractional_delta(coord_a, coord_b)) > tolerance:
            return False
    return True


def files_equal(a: Path, b: Path, chunk_size: int = 1024 * 1024) -> bool:
    if a.stat().st_size != b.stat().st_size:
        return False
    with a.open("rb") as left, b.open("rb") as right:
        while True:
            left_chunk = left.read(chunk_size)
            right_chunk = right.read(chunk_size)
            if left_chunk != right_chunk:
                return False
            if not left_chunk:
                return True


def normalized_formula(structure: PoscarData) -> str:
    counts = Counter(structure.site_elements)
    divisor = reduce(math.gcd, counts.values())
    reduced = {element: count // divisor for element, count in counts.items()}
    if reduced != {"Zn": 1, "In": 2, "S": 4}:
        fail(f"Unexpected reduced composition: {reduced}")
    return "ZnIn2S4"


def formula_units(structure: PoscarData) -> int:
    counts = Counter(structure.site_elements)
    normalized_formula(structure)
    return counts["Zn"]


def cell_volume(structure: PoscarData) -> float:
    return abs(dot(structure.lattice[0], cross(structure.lattice[1], structure.lattice[2])))


def density_g_cm3(structure: PoscarData) -> float:
    mass_amu = sum(ATOMIC_MASS[element] for element in structure.site_elements)
    return mass_amu * 1.66053906660 / cell_volume(structure)


def parse_eigenval(path: Path) -> EigenvalData:
    lines = path.read_text().splitlines()
    if len(lines) < 7:
        fail(f"Incomplete EIGENVAL: {path}")
    nelect_raw, nkpoints_raw, nbands_raw = lines[5].split()[:3]
    nelect = int(round(float(nelect_raw)))
    nkpoints = int(nkpoints_raw)
    nbands = int(nbands_raw)
    index = 6
    kpoints: list[tuple[float, float, float]] = []
    weights: list[float] = []
    energies: list[list[float]] = []
    occupancies: list[list[float]] = []
    for _ in range(nkpoints):
        while index < len(lines) and not lines[index].strip():
            index += 1
        if index >= len(lines):
            fail(f"Missing k-point block in {path}")
        k_values = [float(value) for value in lines[index].split()[:4]]
        index += 1
        if len(k_values) != 4:
            fail(f"Malformed k-point in {path}")
        kpoints.append(tuple(k_values[:3]))
        weights.append(k_values[3])
        k_energies: list[float] = []
        k_occupancies: list[float] = []
        for expected_band in range(1, nbands + 1):
            fields = lines[index].split()
            index += 1
            if int(fields[0]) != expected_band or len(fields) < 3:
                fail(f"Malformed band block in {path}")
            k_energies.append(float(fields[1]))
            k_occupancies.append(float(fields[2]))
        energies.append(k_energies)
        occupancies.append(k_occupancies)
    return EigenvalData(nelect, nkpoints, nbands, kpoints, weights, energies, occupancies)


def parse_uniform_kpoints(path: Path) -> tuple[tuple[int, int, int], tuple[float, float, float]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 5 or lines[1] != "0" or not lines[2].lower().startswith("g"):
        fail(f"Expected automatic Gamma-centred KPOINTS: {path}")
    mesh = tuple(int(value) for value in lines[3].split()[:3])
    shift = tuple(float(value) for value in lines[4].split()[:3])
    return mesh, shift  # type: ignore[return-value]


def parse_path_kpoints(path: Path) -> PathData:
    lines = path.read_text().splitlines()
    if len(lines) < 6:
        fail(f"Incomplete Line-mode KPOINTS: {path}")
    points_per_segment = int(lines[1].strip())
    if not lines[2].strip().lower().startswith("line"):
        fail(f"Expected Line-mode KPOINTS: {path}")
    if not lines[3].strip().lower().startswith("rec"):
        fail(f"Expected reciprocal Line-mode coordinates: {path}")
    entries = []
    for line in lines[4:]:
        if not line.strip():
            continue
        coordinate_text, separator, label = line.partition("!")
        if not separator or not label.strip():
            fail(f"Every path endpoint must be labelled: {path}")
        coordinate = tuple(float(value) for value in coordinate_text.split()[:3])
        if len(coordinate) != 3:
            fail(f"Malformed path endpoint: {path}")
        entries.append((coordinate, label.strip()))
    if len(entries) % 2:
        fail(f"Unpaired path endpoint: {path}")
    segments = []
    for index in range(0, len(entries), 2):
        start, end = entries[index], entries[index + 1]
        segments.append((start[0], start[1], end[0], end[1]))
    return PathData(points_per_segment, segments)


def parse_procar_header(path: Path) -> tuple[int, int, int]:
    with path.open(errors="replace") as handle:
        header = "".join(next(handle) for _ in range(3))
    match = re.search(
        r"# of k-points:\s*(\d+)\s+# of bands:\s*(\d+)\s+# of ions:\s*(\d+)",
        header,
    )
    if not match:
        fail(f"Cannot parse PROCAR header: {path}")
    return tuple(int(value) for value in match.groups())  # type: ignore[return-value]


def validate_procar_edge_bands(
    path: Path, eigenval: EigenvalData, edge_data: dict[str, object]
) -> None:
    """Independently match edge energies/occupations in PROCAR and EIGENVAL."""
    valence_band = int(edge_data["valence_index"]) + 1
    conduction_band = int(edge_data["conduction_index"]) + 1
    kpoint_pattern = re.compile(r"^\s*k-point\s+(\d+)\s*:")
    number_pattern = re.compile(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
    )
    band_pattern = re.compile(
        r"^\s*band\s+(\d+)\s+# energy\s+([-+0-9.Ee]+)\s+# occ\.\s+([-+0-9.Ee]+)"
    )
    current_kpoint = None
    matched = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            kpoint_match = kpoint_pattern.match(line)
            if kpoint_match:
                current_kpoint = int(kpoint_match.group(1)) - 1
                if not 0 <= current_kpoint < eigenval.nkpoints:
                    fail(f"PROCAR k-point index is out of range: {path}")
                coordinate_text = line.split(":", 1)[1].split("weight", 1)[0]
                coordinate_values = number_pattern.findall(coordinate_text)
                if len(coordinate_values) != 3:
                    fail(f"Cannot parse PROCAR k-point coordinates: {path}")
                coordinate = tuple(float(value) for value in coordinate_values)
                if max(
                    abs(a - b)
                    for a, b in zip(coordinate, eigenval.kpoints[current_kpoint])
                ) > 6.0e-7:
                    fail(
                        f"PROCAR/EIGENVAL k-point mismatch in {path}: "
                        f"k={current_kpoint + 1}, PROCAR={coordinate}, "
                        f"EIGENVAL={eigenval.kpoints[current_kpoint]}"
                    )
                continue
            band_match = band_pattern.match(line)
            if not band_match or current_kpoint is None:
                continue
            band_number = int(band_match.group(1))
            if band_number not in {valence_band, conduction_band}:
                continue
            energy = float(band_match.group(2))
            occupancy = float(band_match.group(3))
            expected_energy = eigenval.energies[current_kpoint][band_number - 1]
            if abs(energy - expected_energy) > 2.0e-6:
                fail(
                    f"PROCAR/EIGENVAL edge-energy mismatch in {path}: "
                    f"k={current_kpoint + 1}, band={band_number}, PROCAR={energy:.8f}, "
                    f"EIGENVAL={expected_energy:.8f}, "
                    f"delta={abs(energy - expected_energy):.3e} eV"
                )
            if band_number == valence_band and occupancy < 1.99:
                fail(f"PROCAR valence-edge band is not fully occupied: {path}")
            if band_number == conduction_band and occupancy > 0.01:
                fail(f"PROCAR conduction-edge band is occupied: {path}")
            matched += 1
    if matched != 2 * eigenval.nkpoints:
        fail(
            f"PROCAR edge-band coverage mismatch in {path}: "
            f"expected {2 * eigenval.nkpoints}, found {matched}"
        )


def scan_outcar(path: Path) -> dict[str, float | int | bool]:
    result: dict[str, float | int | bool] = {
        "normal_termination": False,
        "ediff_reached": False,
        "fatal_marker": False,
    }
    with path.open(errors="replace") as handle:
        for line in handle:
            if "General timing and accounting informations" in line:
                result["normal_termination"] = True
            if "aborting loop because EDIFF is reached" in line:
                result["ediff_reached"] = True
            if any(marker in line for marker in ("VERY BAD NEWS", "ZBRENT: fatal error", "internal error")):
                result["fatal_marker"] = True
            match = re.search(r"NKPTS\s*=\s*(\d+).*NBANDS=\s*(\d+)", line)
            if match:
                result["nkpoints"] = int(match.group(1))
                result["nbands"] = int(match.group(2))
            match = re.search(r"NELECT\s*=\s*([-+0-9.Ee]+)", line)
            if match:
                result["nelect"] = float(match.group(1))
            match = re.search(r"E-fermi\s*:\s*([-+0-9.Ee]+)", line)
            if match:
                result["efermi"] = float(match.group(1))
    return result


def scan_vasprun(path: Path) -> dict[str, float | bool]:
    efermi = None
    with path.open(errors="replace") as handle:
        for line in handle:
            match = re.search(r'<i name="efermi">\s*([-+0-9.Ee]+)', line)
            if match:
                efermi = float(match.group(1))
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size - 4096))
        complete = b"</modeling>" in handle.read()
    return {"complete": complete, "efermi": efermi if efermi is not None else math.nan}


def final_electronic_iteration(path: Path) -> int:
    iteration = None
    pattern = re.compile(r"^\s*(?:DAV|RMM):\s*(\d+)")
    for line in path.read_text(errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            iteration = int(match.group(1))
    if iteration is None:
        fail(f"No electronic iteration found in {path}")
    return iteration


def parse_doscar(path: Path, site_elements: Sequence[str]) -> DoscarData:
    with path.open(errors="replace") as handle:
        header = [next(handle) for _ in range(6)]
        nions = int(header[0].split()[0])
        values = header[5].split()
        if len(values) < 5:
            fail(f"Malformed DOSCAR header: {path}")
        nedos = int(values[2])
        efermi = float(values[3])
        energies = []
        total_dos = []
        integrated_dos = []
        for _ in range(nedos):
            fields = handle.readline().split()
            if len(fields) != 3:
                fail(f"Expected non-spin-polarized total DOS in {path}")
            energies.append(float(fields[0]))
            total_dos.append(float(fields[1]))
            integrated_dos.append(float(fields[2]))
        if nions != len(site_elements):
            fail(f"DOSCAR/POSCAR atom-count mismatch: {path}")
        element_orbitals = {
            (element, orbital): [0.0] * nedos
            for element in ELEMENT_ORDER
            for orbital in ORBITAL_ORDER
        }
        site_orbitals = []
        for site_index, element in enumerate(site_elements):
            block_header = handle.readline().split()
            if len(block_header) < 5 or int(block_header[2]) != nedos:
                fail(f"Malformed projected-DOS block {site_index + 1} in {path}")
            site_data = {orbital: [] for orbital in ORBITAL_ORDER}
            for energy_index in range(nedos):
                fields = [float(value) for value in handle.readline().split()]
                if len(fields) != 10:
                    fail(
                        "Expected LORBIT=11 non-spin s/py/pz/px/dxy/dyz/dz2/dxz/dx2-y2 "
                        f"columns in {path}; found {len(fields) - 1} projected columns"
                    )
                if abs(fields[0] - energies[energy_index]) > 1.1e-3:
                    fail(f"Projected-DOS energy-grid mismatch in {path}")
                orbital_values = {
                    "s": fields[1],
                    "p": sum(fields[2:5]),
                    "d": sum(fields[5:10]),
                }
                for orbital, value in orbital_values.items():
                    site_data[orbital].append(value)
                    element_orbitals[(element, orbital)][energy_index] += value
            site_orbitals.append(site_data)
    return DoscarData(
        nions,
        nedos,
        efermi,
        energies,
        total_dos,
        integrated_dos,
        element_orbitals,
        site_orbitals,
    )


def integrate_window(
    x: Sequence[float], y: Sequence[float], lower: float, upper: float
) -> float:
    if upper <= lower:
        fail("Integration window must have positive width")
    area = 0.0
    for index in range(len(x) - 1):
        left_x, right_x = x[index], x[index + 1]
        if right_x <= lower or left_x >= upper:
            continue
        start = max(left_x, lower)
        stop = min(right_x, upper)
        if stop <= start:
            continue
        span = right_x - left_x
        start_y = y[index] + (y[index + 1] - y[index]) * (start - left_x) / span
        stop_y = y[index] + (y[index + 1] - y[index]) * (stop - left_x) / span
        area += 0.5 * (start_y + stop_y) * (stop - start)
    return area


def interpolate(x: Sequence[float], y: Sequence[float], target: float) -> float:
    if target <= x[0]:
        return y[0]
    if target >= x[-1]:
        return y[-1]
    for index in range(len(x) - 1):
        if x[index] <= target <= x[index + 1]:
            fraction = (target - x[index]) / (x[index + 1] - x[index])
            return y[index] + fraction * (y[index + 1] - y[index])
    fail("Interpolation target is outside the grid")


def band_edges(eigenval: EigenvalData) -> dict[str, object]:
    if eigenval.nelect % 2:
        fail("Odd NELECT is incompatible with the verified ISPIN=1 setup")
    occupied_bands = eigenval.nelect // 2
    if occupied_bands >= eigenval.nbands:
        fail("No unoccupied band is available in EIGENVAL")
    valence_index = occupied_bands - 1
    conduction_index = occupied_bands
    valence = [bands[valence_index] for bands in eigenval.energies]
    conduction = [bands[conduction_index] for bands in eigenval.energies]
    if min(row[valence_index] for row in eigenval.occupancies) < 0.99:
        fail("Nominal valence-edge band is not fully occupied")
    if max(row[conduction_index] for row in eigenval.occupancies) > 0.01:
        fail("Nominal conduction-edge band has non-negligible occupation")
    vbm = max(valence)
    cbm = min(conduction)
    fundamental_gap = cbm - vbm
    direct_gaps = [cb - vb for vb, cb in zip(valence, conduction)]
    minimum_direct_gap = min(direct_gaps)
    delta = max(0.0, minimum_direct_gap - fundamental_gap)
    vbm_indices = [index for index, value in enumerate(valence) if abs(value - vbm) <= EDGE_TOL_EV]
    cbm_indices = [index for index, value in enumerate(conduction) if abs(value - cbm) <= EDGE_TOL_EV]
    direct_indices = [
        index
        for index, value in enumerate(direct_gaps)
        if abs(value - minimum_direct_gap) <= EDGE_TOL_EV
    ]
    common_edge_indices = sorted(set(vbm_indices).intersection(cbm_indices))
    if fundamental_gap <= 0:
        classification = "METALLIC / ZERO-GAP ON THE SAMPLED UNIFORM MESH"
    elif common_edge_indices and delta <= DIRECT_EXACT_TOL_EV:
        classification = "DIRECT ON THE SAMPLED UNIFORM MESH"
    elif delta <= DIRECT_TOL_EV:
        classification = "AMBIGUOUS / NEAR-DEGENERATE"
    else:
        classification = "INDIRECT ON THE SAMPLED UNIFORM MESH"
    return {
        "occupied_bands": occupied_bands,
        "valence_index": valence_index,
        "conduction_index": conduction_index,
        "valence": valence,
        "conduction": conduction,
        "vbm": vbm,
        "cbm": cbm,
        "fundamental_gap": fundamental_gap,
        "minimum_direct_gap": minimum_direct_gap,
        "delta_direct": delta,
        "direct_gaps": direct_gaps,
        "vbm_indices": vbm_indices,
        "cbm_indices": cbm_indices,
        "direct_indices": direct_indices,
        "common_edge_indices": common_edge_indices,
        "classification": classification,
    }


def validate_path_sampling(
    path_data: PathData, eigenval: EigenvalData, reciprocal: Sequence[Sequence[float]]
) -> list[dict[str, object]]:
    expected_nkpoints = path_data.points_per_segment * len(path_data.segments)
    if eigenval.nkpoints != expected_nkpoints:
        fail(
            f"Line-mode count mismatch: expected {expected_nkpoints}, "
            f"found {eigenval.nkpoints}"
        )
    descriptors = []
    cumulative_offset = 0.0
    for segment_index, (start, start_label, end, end_label) in enumerate(path_data.segments):
        segment_vector = vector_sub(end, start)
        segment_length = norm(frac_to_cart(segment_vector, reciprocal))
        for point_index in range(path_data.points_per_segment):
            fraction = point_index / (path_data.points_per_segment - 1)
            expected = vector_add(start, vector_scale(fraction, segment_vector))
            global_index = segment_index * path_data.points_per_segment + point_index
            actual = eigenval.kpoints[global_index]
            mismatch = norm(frac_to_cart(vector_sub(actual, expected), reciprocal))
            if mismatch > 2.0e-5:
                fail(
                    f"Line-mode interpolation mismatch at k-point {global_index + 1}: "
                    f"{mismatch:.3e} Angstrom^-1"
                )
            if point_index == 0:
                location = start_label
            elif point_index == path_data.points_per_segment - 1:
                location = end_label
            else:
                near = ""
                if fraction <= 0.05:
                    near = f", near {start_label}"
                elif fraction >= 0.95:
                    near = f", near {end_label}"
                location = f"{start_label}-{end_label} (t={fraction:.3f}{near})"
            descriptors.append(
                {
                    "index": global_index,
                    "segment_index": segment_index,
                    "segment": f"{start_label}-{end_label}",
                    "fraction": fraction,
                    "location": location,
                    "distance": cumulative_offset + fraction * segment_length,
                }
            )
        cumulative_offset += segment_length
    return descriptors


def unique_path_locations(indices: Iterable[int], descriptors: Sequence[dict[str, object]]) -> list[str]:
    locations = []
    for index in indices:
        location = str(descriptors[index]["location"])
        if location not in locations:
            locations.append(location)
    return locations


def compact_path_locations(
    indices: Iterable[int], descriptors: Sequence[dict[str, object]]
) -> list[str]:
    """Summarize tied path extrema without dumping every point of a flat segment."""
    selected = sorted(set(indices))
    by_segment: defaultdict[int, list[int]] = defaultdict(list)
    for index in selected:
        by_segment[int(descriptors[index]["segment_index"])].append(index)
    ranges: list[str] = []
    points: list[str] = []
    covered_labels: set[str] = set()
    for segment_index in sorted(by_segment):
        segment_indices = by_segment[segment_index]
        segment_all = [
            index
            for index, descriptor in enumerate(descriptors)
            if int(descriptor["segment_index"]) == segment_index
        ]
        segment = str(descriptors[segment_indices[0]]["segment"])
        start_label, end_label = segment.split("-", 1)
        fractions = [float(descriptors[index]["fraction"]) for index in segment_indices]
        contiguous = all(
            right == left + 1 for left, right in zip(segment_indices, segment_indices[1:])
        )
        if segment_indices == segment_all:
            ranges.append(f"{segment} (entire segment)")
            covered_labels.update((start_label, end_label))
        elif len(segment_indices) > 1 and contiguous:
            start_fraction, end_fraction = min(fractions), max(fractions)
            near = ""
            if end_fraction <= 0.15:
                near = f", near {start_label}"
            elif start_fraction >= 0.85:
                near = f", near {end_label}"
            ranges.append(
                f"{segment} (t={start_fraction:.3f}-{end_fraction:.3f}{near})"
            )
            if start_fraction == 0.0:
                covered_labels.add(start_label)
            if end_fraction == 1.0:
                covered_labels.add(end_label)
        else:
            for index in segment_indices:
                location = str(descriptors[index]["location"])
                if location not in points:
                    points.append(location)
    compact = ranges + [point for point in points if point not in covered_labels]
    return list(dict.fromkeys(compact))


def high_symmetry_label(
    kpoint: Sequence[float],
    path_data: PathData,
    reciprocal: Sequence[Sequence[float]],
) -> str | None:
    exact_labels = []
    for start, start_label, end, end_label in path_data.segments:
        for coordinate, label in ((start, start_label), (end, end_label)):
            exact_delta = vector_sub(kpoint, coordinate)
            if (
                norm(frac_to_cart(exact_delta, reciprocal)) <= 2.0e-4
                and label not in exact_labels
            ):
                exact_labels.append(label)
    return "/".join(exact_labels) if exact_labels else None


def format_kpoint(kpoint: Sequence[float]) -> str:
    return f"({kpoint[0]:.6f}, {kpoint[1]:.6f}, {kpoint[2]:.6f})"


def uniform_locations(
    indices: Sequence[int],
    eigenval: EigenvalData,
    path_data: PathData,
    reciprocal: Sequence[Sequence[float]],
) -> list[str]:
    locations = []
    for index in indices:
        label = high_symmetry_label(eigenval.kpoints[index], path_data, reciprocal)
        location = label if label else format_kpoint(eigenval.kpoints[index])
        if location not in locations:
            locations.append(location)
    return locations


def element_orbital_integrals(
    doscar: DoscarData, vbm: float, cbm: float
) -> dict[tuple[str, float, str, str], tuple[float, float]]:
    raw: dict[tuple[str, float, str, str], float] = {}
    for edge, origin, direction in (("VBM", vbm, -1), ("CBM", cbm, 1)):
        for window in WINDOWS_EV:
            lower, upper = (
                (origin - window, origin) if direction < 0 else (origin, origin + window)
            )
            total = 0.0
            for element in ELEMENT_ORDER:
                for orbital in ORBITAL_ORDER:
                    value = integrate_window(
                        doscar.energies,
                        doscar.element_orbitals[(element, orbital)],
                        lower,
                        upper,
                    )
                    raw[(edge, window, element, orbital)] = value
                    total += value
            if total <= 0:
                fail(f"No projected weight in {edge} {window:.2f} eV window")
            for element in ELEMENT_ORDER:
                for orbital in ORBITAL_ORDER:
                    key = (edge, window, element, orbital)
                    raw[key] = max(0.0, raw[key])
    result = {}
    for edge in ("VBM", "CBM"):
        for window in WINDOWS_EV:
            total = sum(
                raw[(edge, window, element, orbital)]
                for element in ELEMENT_ORDER
                for orbital in ORBITAL_ORDER
            )
            for element in ELEMENT_ORDER:
                for orbital in ORBITAL_ORDER:
                    value = raw[(edge, window, element, orbital)]
                    result[(edge, window, element, orbital)] = (value, 100.0 * value / total)
    return result


def dominant_character(
    integrals: dict[tuple[str, float, str, str], tuple[float, float]],
    edge: str,
    window: float = 0.50,
    limit: int = 3,
) -> str:
    entries = [
        (integrals[(edge, window, element, orbital)][1], f"{element}-{orbital}")
        for element in ELEMENT_ORDER
        for orbital in ORBITAL_ORDER
    ]
    entries.sort(reverse=True)
    return ", ".join(f"{label} {percent:.1f}%" for percent, label in entries[:limit])


def coordination_numbers(structure: PoscarData) -> list[int | None]:
    coordination: list[int | None] = []
    sulfur_indices = [
        index for index, element in enumerate(structure.site_elements) if element == "S"
    ]
    for site_index, element in enumerate(structure.site_elements):
        if element not in {"Zn", "In"}:
            coordination.append(None)
            continue
        count = 0
        origin = structure.frac_coords[site_index]
        for sulfur_index in sulfur_indices:
            sulfur = structure.frac_coords[sulfur_index]
            for tx in (-1, 0, 1):
                for ty in (-1, 0, 1):
                    for tz in (-1, 0, 1):
                        delta = (
                            sulfur[0] + tx - origin[0],
                            sulfur[1] + ty - origin[1],
                            sulfur[2] + tz - origin[2],
                        )
                        distance = norm(frac_to_cart(delta, structure.lattice))
                        if 1.0e-8 < distance <= COORDINATION_CUTOFF_ANG:
                            count += 1
        coordination.append(count)
    return coordination


def site_edge_diagnostics(
    doscar: DoscarData,
    structure: PoscarData,
    vbm: float,
    cbm: float,
) -> dict[str, list[tuple[str, float]]]:
    coordination = coordination_numbers(structure)
    element_counter: defaultdict[str, int] = defaultdict(int)
    labels = []
    for index, element in enumerate(structure.site_elements):
        element_counter[element] += 1
        suffix = f", S{coordination[index]}" if coordination[index] is not None else ""
        labels.append(f"{element}{element_counter[element]}{suffix}")
    result = {}
    for edge, lower, upper in (
        ("VBM", vbm - 0.50, vbm),
        ("CBM", cbm, cbm + 0.50),
    ):
        values = []
        for label, site_data in zip(labels, doscar.site_orbitals):
            weight = sum(
                integrate_window(doscar.energies, site_data[orbital], lower, upper)
                for orbital in ORBITAL_ORDER
            )
            values.append((label, weight))
        total = sum(value for _, value in values)
        result[edge] = [(label, 100.0 * value / total) for label, value in values]
    return result


def dos_diagnostics(
    doscar: DoscarData, vbm: float, cbm: float, n_formula_units: int
) -> dict[str, float | bool]:
    threshold = max(1.0e-5, max(doscar.total_dos) * 1.0e-5)
    midpoint = 0.5 * (vbm + cbm)
    occupied_points = [
        (energy, density)
        for energy, density in zip(doscar.energies, doscar.total_dos)
        if energy < midpoint and density > threshold
    ]
    empty_points = [
        (energy, density)
        for energy, density in zip(doscar.energies, doscar.total_dos)
        if energy > midpoint and density > threshold
    ]
    if not occupied_points or not empty_points:
        fail("Cannot locate DOS band onsets")
    valence_onset = max(energy for energy, _ in occupied_points)
    conduction_onset = min(energy for energy, _ in empty_points)
    valence_floor = min(energy for energy, _ in occupied_points)
    interior = [
        density
        for energy, density in zip(doscar.energies, doscar.total_dos)
        if valence_onset < energy < conduction_onset
    ]
    max_gap_dos = max(interior, default=0.0)
    return {
        "threshold": threshold,
        "valence_onset": valence_onset,
        "conduction_onset": conduction_onset,
        "max_gap_dos": max_gap_dos,
        "clean_gap": max_gap_dos <= threshold,
        "valence_width": valence_onset - valence_floor,
        "conduction_0_05_states_per_fu": integrate_window(
            doscar.energies, doscar.total_dos, cbm, cbm + 0.50
        )
        / n_formula_units,
        "conduction_05_10_states_per_fu": integrate_window(
            doscar.energies, doscar.total_dos, cbm + 0.50, cbm + 1.00
        )
        / n_formula_units,
        "integrated_at_fermi": interpolate(
            doscar.energies, doscar.integrated_dos, doscar.efermi
        ),
    }


def segment_dispersions(
    path_data: PathData, edge_data: dict[str, object]
) -> dict[str, tuple[float, float]]:
    valence = edge_data["valence"]
    conduction = edge_data["conduction"]
    assert isinstance(valence, list) and isinstance(conduction, list)
    result = {}
    count = path_data.points_per_segment
    for segment_index, (_, start_label, _, end_label) in enumerate(path_data.segments):
        start = segment_index * count
        stop = start + count
        valence_values = valence[start:stop]
        conduction_values = conduction[start:stop]
        result[f"{start_label}-{end_label}"] = (
            max(valence_values) - min(valence_values),
            max(conduction_values) - min(conduction_values),
        )
    return result


def pearson(values_a: Sequence[float], values_b: Sequence[float]) -> float:
    mean_a = sum(values_a) / len(values_a)
    mean_b = sum(values_b) / len(values_b)
    covariance = sum((a - mean_a) * (b - mean_b) for a, b in zip(values_a, values_b))
    scale_a = math.sqrt(sum((a - mean_a) ** 2 for a in values_a))
    scale_b = math.sqrt(sum((b - mean_b) ** 2 for b in values_b))
    return covariance / (scale_a * scale_b) if scale_a and scale_b else math.nan


def analyse_phase(phase: str) -> dict[str, object]:
    meta = PHASE_META[phase]
    dos_dir = CALC_ROOT / phase / "dos_pdos"
    band_dir = CALC_ROOT / phase / "band"
    stage03_dir = STAGE03_ROOT / phase
    required = ("INCAR", "KPOINTS", "POSCAR", "CONTCAR", "OUTCAR", "vasprun.xml", "DOSCAR", "EIGENVAL", "PROCAR", "CHGCAR")
    for branch_dir in (dos_dir, band_dir):
        missing = [name for name in required if not (branch_dir / name).is_file()]
        if missing:
            fail(f"{phase}/{branch_dir.name} missing: {', '.join(missing)}")

    dos_incar = read_incar(dos_dir / "INCAR")
    band_incar = read_incar(band_dir / "INCAR")
    shared_expected = {
        "ICHARG": 11,
        "ENCUT": 500,
        "ISPIN": 1,
        "IBRION": -1,
        "NSW": 0,
        "LORBIT": 11,
    }
    for branch_name, values in (("dos_pdos", dos_incar), ("band", band_incar)):
        for key, expected in shared_expected.items():
            if incar_int(values, key) != expected:
                fail(f"{phase}/{branch_name}: {key} != {expected}")
        if values.get("GGA", "").upper() != "PE" or incar_int(values, "IVDW") != 12:
            fail(f"{phase}/{branch_name}: expected PBE with IVDW=12 geometry provenance")
        if incar_float(values, "EDIFF") > 1.0e-7:
            fail(f"{phase}/{branch_name}: EDIFF is looser than 1e-7 eV")
    if incar_int(dos_incar, "ISMEAR") != -5 or incar_int(dos_incar, "ISYM") != 2:
        fail(f"{phase}/dos_pdos: unexpected occupation or symmetry setup")
    if incar_int(dos_incar, "NEDOS") != 3000:
        fail(f"{phase}/dos_pdos: expected NEDOS=3000")
    if incar_int(band_incar, "ISMEAR") != 0 or incar_int(band_incar, "ISYM") != 0:
        fail(f"{phase}/band: unexpected occupation or symmetry setup")

    mesh, shift = parse_uniform_kpoints(dos_dir / "KPOINTS")
    if mesh != meta["mesh"] or any(abs(value) > 1.0e-12 for value in shift):
        fail(f"{phase}: unexpected DOS mesh or shift")
    path_data = parse_path_kpoints(band_dir / "KPOINTS")
    actual_segments = tuple((start_label, end_label) for _, start_label, _, end_label in path_data.segments)
    if path_data.points_per_segment != 40 or actual_segments != meta["segments"]:
        fail(f"{phase}: Line-mode path does not match validated {meta['path_type']} topology")

    dos_poscar = parse_poscar(dos_dir / "POSCAR")
    band_poscar = parse_poscar(band_dir / "POSCAR")
    stage03_poscar = parse_poscar(stage03_dir / "POSCAR")
    normalized_formula(dos_poscar)
    if not files_equal(dos_dir / "POSCAR", stage03_dir / "POSCAR"):
        fail(f"{phase}: DOS POSCAR is not byte-identical to Stage 03")
    if not files_equal(band_dir / "POSCAR", stage03_dir / "POSCAR"):
        fail(f"{phase}: band POSCAR is not byte-identical to Stage 03")
    if not files_equal(dos_dir / "CHGCAR", stage03_dir / "CHGCAR"):
        fail(f"{phase}: DOS CHGCAR is not byte-identical to Stage 03")
    if not files_equal(band_dir / "CHGCAR", stage03_dir / "CHGCAR"):
        fail(f"{phase}: band CHGCAR is not byte-identical to Stage 03")
    if not structures_equivalent(dos_poscar, band_poscar):
        fail(f"{phase}: branch structures differ")
    for branch_dir, poscar in ((dos_dir, dos_poscar), (band_dir, band_poscar)):
        if not structures_equivalent(poscar, parse_poscar(branch_dir / "CONTCAR")):
            fail(f"{phase}/{branch_dir.name}: geometry changed")

    dos_eigenval = parse_eigenval(dos_dir / "EIGENVAL")
    band_eigenval = parse_eigenval(band_dir / "EIGENVAL")
    if dos_eigenval.nelect != band_eigenval.nelect:
        fail(f"{phase}: branch electron counts differ")
    uniform_edge = band_edges(dos_eigenval)
    path_edge = band_edges(band_eigenval)
    reciprocal = reciprocal_lattice(dos_poscar.lattice)
    path_descriptors = validate_path_sampling(path_data, band_eigenval, reciprocal)

    branch_validation = {}
    for branch_name, branch_dir, values, eigenval, poscar, edge_data in (
        ("dos_pdos", dos_dir, dos_incar, dos_eigenval, dos_poscar, uniform_edge),
        ("band", band_dir, band_incar, band_eigenval, band_poscar, path_edge),
    ):
        outcar = scan_outcar(branch_dir / "OUTCAR")
        vasprun = scan_vasprun(branch_dir / "vasprun.xml")
        iteration = final_electronic_iteration(branch_dir / "OSZICAR")
        procar_shape = parse_procar_header(branch_dir / "PROCAR")
        if not outcar["normal_termination"] or not outcar["ediff_reached"] or outcar["fatal_marker"]:
            fail(f"{phase}/{branch_name}: OUTCAR termination/convergence check failed")
        if not vasprun["complete"]:
            fail(f"{phase}/{branch_name}: vasprun.xml is incomplete")
        if iteration >= incar_int(values, "NELM"):
            fail(f"{phase}/{branch_name}: electronic iteration reached NELM")
        expected_shape = (eigenval.nkpoints, eigenval.nbands, len(poscar.site_elements))
        if procar_shape != expected_shape:
            fail(f"{phase}/{branch_name}: PROCAR shape mismatch {procar_shape} != {expected_shape}")
        validate_procar_edge_bands(branch_dir / "PROCAR", eigenval, edge_data)
        if int(outcar.get("nkpoints", -1)) != eigenval.nkpoints or int(outcar.get("nbands", -1)) != eigenval.nbands:
            fail(f"{phase}/{branch_name}: OUTCAR/EIGENVAL shape mismatch")
        if abs(float(outcar.get("nelect", math.nan)) - eigenval.nelect) > 1.0e-6:
            fail(f"{phase}/{branch_name}: OUTCAR/EIGENVAL NELECT mismatch")
        if abs(float(outcar.get("efermi", math.nan)) - float(vasprun["efermi"])) > 1.0e-3:
            fail(f"{phase}/{branch_name}: OUTCAR/vasprun Fermi-level mismatch")
        branch_validation[branch_name] = {
            "iteration": iteration,
            "efermi": float(vasprun["efermi"]),
            "nkpoints": eigenval.nkpoints,
            "nbands": eigenval.nbands,
            "nions": len(poscar.site_elements),
        }

    doscar = parse_doscar(dos_dir / "DOSCAR", dos_poscar.site_elements)
    if abs(doscar.efermi - float(branch_validation["dos_pdos"]["efermi"])) > 1.0e-3:
        fail(f"{phase}: DOSCAR/vasprun Fermi-level mismatch")
    if doscar.nedos != 3000 or doscar.nions != len(dos_poscar.site_elements):
        fail(f"{phase}: DOSCAR dimensions are inconsistent")
    integrated_at_fermi = interpolate(doscar.energies, doscar.integrated_dos, doscar.efermi)
    if abs(integrated_at_fermi - dos_eigenval.nelect) > 0.5:
        fail(
            f"{phase}: integrated DOS at E-fermi ({integrated_at_fermi:.3f}) "
            f"does not recover NELECT ({dos_eigenval.nelect})"
        )

    vbm = float(uniform_edge["vbm"])
    cbm = float(uniform_edge["cbm"])
    pdos_integrals = element_orbital_integrals(doscar, vbm, cbm)
    dos_info = dos_diagnostics(doscar, vbm, cbm, formula_units(dos_poscar))
    site_info = site_edge_diagnostics(doscar, dos_poscar, vbm, cbm)
    dispersions = segment_dispersions(path_data, path_edge)

    return {
        "phase": phase,
        "meta": meta,
        "structure": dos_poscar,
        "formula_units": formula_units(dos_poscar),
        "density": density_g_cm3(dos_poscar),
        "mesh": mesh,
        "dos_eigenval": dos_eigenval,
        "band_eigenval": band_eigenval,
        "uniform_edge": uniform_edge,
        "path_edge": path_edge,
        "path_data": path_data,
        "path_descriptors": path_descriptors,
        "reciprocal": reciprocal,
        "branch_validation": branch_validation,
        "doscar": doscar,
        "dos_info": dos_info,
        "pdos_integrals": pdos_integrals,
        "site_info": site_info,
        "dispersions": dispersions,
        "acceptance": (
            "PASS WITH CAVEAT"
            if uniform_edge["classification"] == "AMBIGUOUS / NEAR-DEGENERATE"
            else "PASS"
        ),
    }


def write_summary_csv(results: Sequence[dict[str, object]]) -> None:
    fields = [
        "structure",
        "space_group",
        "acceptance",
        "electronic_structure_method",
        "dos_run_mode",
        "dos_mesh",
        "uniform_irreducible_kpoints",
        "dos_efermi_eV",
        "density_g_cm3",
        "pbe_fundamental_gap_eV",
        "pbe_min_direct_gap_eV",
        "direct_minus_fundamental_meV",
        "sampled_gap_classification",
        "uniform_VBM_energy_eV",
        "uniform_CBM_energy_eV",
        "uniform_VBM_kx",
        "uniform_VBM_ky",
        "uniform_VBM_kz",
        "uniform_CBM_kx",
        "uniform_CBM_ky",
        "uniform_CBM_kz",
        "band_path_gap_eV",
        "band_path_minus_uniform_gap_meV",
        "band_path_VBM_location",
        "band_path_CBM_location",
        "band_path_type",
        "band_path",
        "verdict",
    ]
    with SUMMARY_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            uniform = result["uniform_edge"]
            path_edge = result["path_edge"]
            dos_eigenval = result["dos_eigenval"]
            band_eigenval = result["band_eigenval"]
            path_data = result["path_data"]
            reciprocal = result["reciprocal"]
            descriptors = result["path_descriptors"]
            meta = result["meta"]
            assert isinstance(uniform, dict) and isinstance(path_edge, dict)
            assert isinstance(dos_eigenval, EigenvalData) and isinstance(band_eigenval, EigenvalData)
            assert isinstance(path_data, PathData) and isinstance(meta, dict)
            vbm_index = int(uniform["vbm_indices"][0])
            cbm_index = int(uniform["cbm_indices"][0])
            path_vbm_locations = compact_path_locations(path_edge["vbm_indices"], descriptors)
            path_cbm_locations = compact_path_locations(path_edge["cbm_indices"], descriptors)
            writer.writerow(
                {
                    "structure": result["phase"],
                    "space_group": meta["space_group"],
                    "acceptance": result["acceptance"],
                    "electronic_structure_method": "PBE on PBE+D3(BJ)-relaxed geometry",
                    "dos_run_mode": "fixed-density non-SCF (ICHARG=11)",
                    "dos_mesh": "x".join(str(value) for value in result["mesh"]),
                    "uniform_irreducible_kpoints": dos_eigenval.nkpoints,
                    "dos_efermi_eV": f"{result['doscar'].efermi:.6f}",
                    "density_g_cm3": f"{result['density']:.6f}",
                    "pbe_fundamental_gap_eV": f"{uniform['fundamental_gap']:.6f}",
                    "pbe_min_direct_gap_eV": f"{uniform['minimum_direct_gap']:.6f}",
                    "direct_minus_fundamental_meV": f"{1000.0 * uniform['delta_direct']:.3f}",
                    "sampled_gap_classification": uniform["classification"],
                    "uniform_VBM_energy_eV": f"{uniform['vbm']:.6f}",
                    "uniform_CBM_energy_eV": f"{uniform['cbm']:.6f}",
                    "uniform_VBM_kx": f"{dos_eigenval.kpoints[vbm_index][0]:.8f}",
                    "uniform_VBM_ky": f"{dos_eigenval.kpoints[vbm_index][1]:.8f}",
                    "uniform_VBM_kz": f"{dos_eigenval.kpoints[vbm_index][2]:.8f}",
                    "uniform_CBM_kx": f"{dos_eigenval.kpoints[cbm_index][0]:.8f}",
                    "uniform_CBM_ky": f"{dos_eigenval.kpoints[cbm_index][1]:.8f}",
                    "uniform_CBM_kz": f"{dos_eigenval.kpoints[cbm_index][2]:.8f}",
                    "band_path_gap_eV": f"{path_edge['fundamental_gap']:.6f}",
                    "band_path_minus_uniform_gap_meV": (
                        f"{1000.0 * (path_edge['fundamental_gap'] - uniform['fundamental_gap']):.3f}"
                    ),
                    "band_path_VBM_location": "; ".join(path_vbm_locations),
                    "band_path_CBM_location": "; ".join(path_cbm_locations),
                    "band_path_type": meta["path_type"],
                    "band_path": meta["path"],
                    "verdict": result["acceptance"],
                }
            )


def write_band_edges_csv(results: Sequence[dict[str, object]]) -> None:
    fields = [
        "structure",
        "edge",
        "energy_eV",
        "kx",
        "ky",
        "kz",
        "source",
        "high_symmetry_label_or_segment",
        "notes",
    ]
    with EDGES_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            for source, eigen_key, edge_key in (
                ("dense_uniform_mesh", "dos_eigenval", "uniform_edge"),
                ("high_symmetry_path", "band_eigenval", "path_edge"),
            ):
                eigenval = result[eigen_key]
                edge = result[edge_key]
                assert isinstance(eigenval, EigenvalData) and isinstance(edge, dict)
                if source == "dense_uniform_mesh":
                    location_for = lambda index: high_symmetry_label(  # noqa: E731
                        eigenval.kpoints[index], result["path_data"], result["reciprocal"]
                    ) or "unlabelled uniform-mesh point"
                else:
                    location_for = lambda index: result["path_descriptors"][index]["location"]  # noqa: E731
                for edge_name, value_key, index_key in (
                    ("VBM", "vbm", "vbm_indices"),
                    ("CBM", "cbm", "cbm_indices"),
                ):
                    seen = set()
                    for index in edge[index_key]:
                        kpoint = eigenval.kpoints[index]
                        record_key = (tuple(round(value, 8) for value in kpoint), str(location_for(index)))
                        if record_key in seen:
                            continue
                        seen.add(record_key)
                        notes = (
                            f"irreducible-mesh representative; weight={eigenval.weights[index]:.8f}"
                            if source == "dense_uniform_mesh"
                            else (
                                f"path index={index + 1}; cumulative distance="
                                f"{result['path_descriptors'][index]['distance']:.8f} Angstrom^-1"
                            )
                        )
                        writer.writerow(
                            {
                                "structure": result["phase"],
                                "edge": edge_name,
                                "energy_eV": f"{edge[value_key]:.6f}",
                                "kx": f"{kpoint[0]:.8f}",
                                "ky": f"{kpoint[1]:.8f}",
                                "kz": f"{kpoint[2]:.8f}",
                                "source": source,
                                "high_symmetry_label_or_segment": location_for(index),
                                "notes": notes,
                            }
                        )
                if source == "dense_uniform_mesh":
                    seen_direct = set()
                    for index in edge["direct_indices"]:
                        kpoint = eigenval.kpoints[index]
                        record_key = tuple(round(value, 8) for value in kpoint)
                        if record_key in seen_direct:
                            continue
                        seen_direct.add(record_key)
                        writer.writerow(
                            {
                                "structure": result["phase"],
                                "edge": "minimum_direct_gap",
                                "energy_eV": f"{edge['minimum_direct_gap']:.6f}",
                                "kx": f"{kpoint[0]:.8f}",
                                "ky": f"{kpoint[1]:.8f}",
                                "kz": f"{kpoint[2]:.8f}",
                                "source": source,
                                "high_symmetry_label_or_segment": location_for(index),
                                "notes": "energy_eV is E_CBM(k)-E_VBM(k)",
                            }
                        )


def write_pdos_csv(results: Sequence[dict[str, object]]) -> None:
    fields = [
        "structure",
        "edge",
        "window_eV",
        "element",
        "orbital",
        "integrated_weight",
        "normalized_percent",
    ]
    with PDOS_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            integrals = result["pdos_integrals"]
            assert isinstance(integrals, dict)
            for edge in ("VBM", "CBM"):
                for window in WINDOWS_EV:
                    for element in ELEMENT_ORDER:
                        for orbital in ORBITAL_ORDER:
                            value, percent = integrals[(edge, window, element, orbital)]
                            writer.writerow(
                                {
                                    "structure": result["phase"],
                                    "edge": edge,
                                    "window_eV": f"{window:.2f}",
                                    "element": element,
                                    "orbital": orbital,
                                    "integrated_weight": f"{value:.8f}",
                                    "normalized_percent": f"{percent:.5f}",
                                }
                            )


def compact_uniform_location(result: dict[str, object], edge_name: str) -> str:
    edge = result["uniform_edge"]
    eigenval = result["dos_eigenval"]
    assert isinstance(edge, dict) and isinstance(eigenval, EigenvalData)
    locations = uniform_locations(
        edge[f"{edge_name.lower()}_indices"],
        eigenval,
        result["path_data"],
        result["reciprocal"],
    )
    return "; ".join(locations)


def dispersion_sentence(result: dict[str, object]) -> str:
    phase = str(result["phase"])
    dispersions = result["dispersions"]
    assert isinstance(dispersions, dict)
    if phase in {"beta", "IIa_prime", "IIb"}:
        in_plane_v, in_plane_c = dispersions["GAMMA-M"]
        out_plane_v, out_plane_c = dispersions["GAMMA-A"]
        return (
            f"On GAMMA-M, the top-valence/bottom-conduction band ranges are "
            f"{in_plane_v:.3f}/{in_plane_c:.3f} eV; on GAMMA-A they are "
            f"{out_plane_v:.3f}/{out_plane_c:.3f} eV. These path-wise ranges "
            "indicate anisotropy but are not effective masses."
        )
    if phase == "alpha1":
        v1, c1 = dispersions["GAMMA-T"]
        v2, c2 = dispersions["L-GAMMA"]
        return (
            f"The top-valence/bottom-conduction ranges are {v1:.3f}/{c1:.3f} eV "
            f"on GAMMA-T and {v2:.3f}/{c2:.3f} eV on L-GAMMA; the rhombohedral "
            "path geometry prevents treating this as a simple Cartesian mass comparison."
        )
    v1, c1 = dispersions["GAMMA-X"]
    v2, c2 = dispersions["GAMMA-L"]
    return (
        f"The top-valence/bottom-conduction ranges are {v1:.3f}/{c1:.3f} eV "
        f"on GAMMA-X and {v2:.3f}/{c2:.3f} eV on GAMMA-L."
    )


def site_note(result: dict[str, object]) -> str:
    site_info = result["site_info"]
    assert isinstance(site_info, dict)
    notes = []
    for edge in ("VBM", "CBM"):
        entries = site_info[edge]
        by_element: defaultdict[str, list[tuple[str, float]]] = defaultdict(list)
        for label, percent in entries:
            match = re.match(r"[A-Za-z]+", label)
            if not match:
                fail(f"Cannot parse site label {label}")
            by_element[match.group(0)].append((label, percent))
        for element, element_entries in by_element.items():
            total = sum(percent for _, percent in element_entries)
            if total < 10.0 or len(element_entries) < 2:
                continue
            top_label, top_percent = max(element_entries, key=lambda item: item[1])
            same_element_average = total / len(element_entries)
            if top_percent >= 20.0 or (
                top_percent >= 5.0 and top_percent >= 2.5 * same_element_average
            ):
                notes.append(
                    f"{edge}: {top_label} carries {top_percent:.1f}% of total projected edge weight"
                )
    if not notes:
        return "No individual same-element site dominates the 0.50 eV edge windows."
    return "; ".join(notes) + "."


def write_report(results: Sequence[dict[str, object]]) -> None:
    gaps = [float(result["uniform_edge"]["fundamental_gap"]) for result in results]
    lowest = min(results, key=lambda result: float(result["uniform_edge"]["fundamental_gap"]))
    highest = max(results, key=lambda result: float(result["uniform_edge"]["fundamental_gap"]))
    classifications = Counter(str(result["uniform_edge"]["classification"]) for result in results)
    vbm_top = [dominant_character(result["pdos_integrals"], "VBM", limit=1).split()[0] for result in results]
    acceptance_counts = Counter(str(result["acceptance"]) for result in results)
    lines = [
        "# WP1 ZnIn2S4 PBE Electronic-Structure Analysis",
        "",
        "## Executive Summary",
        "",
        (
            "All 10 Stage 04 calculations (five dense DOS/PDOS branches and five "
            "high-symmetry band branches) passed termination, electronic-convergence, "
            "fixed-geometry, structure/charge-density-lineage, dimensional, and parseability checks."
        ),
        "",
        (
            f"The dense-mesh sampled PBE fundamental gaps span {min(gaps):.3f}-{max(gaps):.3f} eV. "
            f"{PHASE_META[str(lowest['phase'])]['display']} has the smallest sampled gap "
            f"({lowest['uniform_edge']['fundamental_gap']:.3f} eV), and "
            f"{PHASE_META[str(highest['phase'])]['display']} has the largest "
            f"({highest['uniform_edge']['fundamental_gap']:.3f} eV)."
        ),
        "",
        (
            f"Sampled classifications comprise {classifications.get('DIRECT ON THE SAMPLED UNIFORM MESH', 0)} direct, "
            f"{classifications.get('INDIRECT ON THE SAMPLED UNIFORM MESH', 0)} indirect, and "
            f"{classifications.get('AMBIGUOUS / NEAR-DEGENERATE', 0)} near-degenerate cases. "
            "S-p is the dominant VBM projection in every phase. The layered CBMs are strongly "
            "hybridized S-p/In-s/S-s states whose first-place channel is window-dependent; the "
            "spinel CBM instead has substantial In-s/S-s/Zn-s character."
        ),
        "",
        (
            f"Phase-level acceptance is {acceptance_counts.get('PASS', 0)} PASS, "
            f"{acceptance_counts.get('PASS WITH CAVEAT', 0)} PASS WITH CAVEAT, and 0 REQUIRES "
            "RECALCULATION. No Stage 04 recalculation is required before HSE06. Stage 04 can be frozen, "
            "with the finite uniform-mesh extremum and PBE gap-magnitude limitations carried "
            "forward explicitly into Stage 05."
        ),
        "",
        (
            "All values below are **PBE electronic structures on PBE+D3(BJ)-relaxed geometries**. "
            "D3(BJ) determined the geometry and is not an electronic band-gap correction."
        ),
        "",
        "## Calculation Validation",
        "",
        (
            "Both branches use the accepted Stage 03 POSCAR and a byte-identical Stage 03 CHGCAR. "
            "`NSW = 0` and `IBRION = -1` freeze the geometry; output CONTCAR coordinates differ "
            "from POSCAR only at numerical round-off. Both DOS/PDOS and band branches are "
            "fixed-density non-self-consistent evaluations (`ICHARG = 11`). This is valid here "
            "because the charge-density/structure lineage is exact, but it must not be described "
            "as a new dense-mesh self-consistent charge calculation."
        ),
        "",
        "| Phase | Branch | k sampling | Irreducible/path k points | Final electronic iteration | Validation |",
        "| --- | --- | --- | ---: | ---: | --- |",
    ]
    for result in results:
        meta = result["meta"]
        validation = result["branch_validation"]
        lines.append(
            f"| {meta['display']} | DOS/PDOS | Gamma {' x '.join(str(value) for value in result['mesh'])} | "
            f"{validation['dos_pdos']['nkpoints']} | {validation['dos_pdos']['iteration']} | PASS |"
        )
        lines.append(
            f"| {meta['display']} | Band | {meta['path_type']} Line-mode, 40/segment | "
            f"{validation['band']['nkpoints']} | {validation['band']['iteration']} | PASS |"
        )
    path_differences = [
        1000.0
        * (
            float(result["path_edge"]["fundamental_gap"])
            - float(result["uniform_edge"]["fundamental_gap"])
        )
        for result in results
    ]
    lines.extend(
        [
            "",
            (
                "All INCARs resolve to PBE (`GGA = PE`), 500 eV, `ISPIN = 1`, `LORBIT = 11`, "
                "and `EDIFF <= 1e-7` eV. DOS/PDOS uses the stated unshifted dense Gamma mesh, "
                "`ISMEAR = -5`, `NEDOS = 3000`, and `ISYM = 2`; band calculations use explicit "
                "reciprocal-coordinate Line-mode, `ISMEAR = 0`, and `ISYM = 0`. Every OUTCAR "
                "contains both the EDIFF and normal timing markers, every vasprun.xml closes "
                "normally, and every EIGENVAL/PROCAR/DOSCAR dimension agrees with the structure. "
                "Edge energies and occupied/unoccupied assignments independently agree between "
                "EIGENVAL and PROCAR at every sampled k point."
            ),
            "",
            "## Band-Gap Summary",
            "",
            "| Phase | Uniform-mesh fundamental gap (eV) | Minimum direct gap (eV) | Delta direct-indirect (meV) | Classification | Sampled VBM | Sampled CBM |",
            "| --- | ---: | ---: | ---: | --- | --- | --- |",
        ]
    )
    for result in results:
        edge = result["uniform_edge"]
        lines.append(
            f"| {result['meta']['display']} | {edge['fundamental_gap']:.3f} | "
            f"{edge['minimum_direct_gap']:.3f} | {1000.0 * edge['delta_direct']:.1f} | "
            f"{edge['classification']} | {compact_uniform_location(result, 'VBM')} | "
            f"{compact_uniform_location(result, 'CBM')} |"
        )
    lines.extend(
        [
            "",
            (
                f"Classification uses the dense uniform mesh. A direct-minus-fundamental difference "
                f"of <= {1000.0 * DIRECT_TOL_EV:.0f} meV is treated as near-degenerate unless the "
                "same sampled k point explicitly hosts both extrema. The CSV retains raw VASP "
                "eigenvalue references; absolute eigenvalues from separate cells are not vacuum-aligned "
                "band positions and must not be compared as such."
            ),
            "",
            (
                "The path-minus-uniform gap differences (spinel, alpha1, beta, IIa-prime, IIb) are "
                + ", ".join(f"{value:+.1f}" for value in path_differences)
                + " meV. The small negative differences for spinel and IIa-prime mean that the finer "
                "Line-mode sampling intersects a slightly higher valence point than the uniform mesh; "
                "they are retained explicitly rather than treated as corrupt data."
            ),
            "",
            "## High-Symmetry Band Structures",
            "",
        ]
    )
    for result in results:
        edge = result["path_edge"]
        descriptors = result["path_descriptors"]
        vbm_locations = "; ".join(compact_path_locations(edge["vbm_indices"], descriptors))
        cbm_locations = "; ".join(compact_path_locations(edge["cbm_indices"], descriptors))
        lines.extend(
            [
                f"### {result['meta']['display']}",
                "",
                f"Path: `{result['meta']['path']}` ({result['meta']['path_type']}; 40 points per segment).",
                "",
                (
                    f"The path-only gap is {edge['fundamental_gap']:.3f} eV, with the path VBM at "
                    f"{vbm_locations} and path CBM at {cbm_locations}. {dispersion_sentence(result)}"
                ),
                "",
            ]
        )
    lines.extend(
        [
            "## DOS / PDOS",
            "",
            (
                "DOS onsets are used only as a qualitative check on the eigenvalue-derived gaps. "
                "Band-edge character is the integrated LORBIT=11 PAW projection over 0.50 eV edge "
                "windows; the machine-readable table also includes 0.20 eV windows. Percentages are "
                "normalized within each projected window and are comparative descriptors, not exact "
                "chemical populations."
            ),
            "",
            "| Phase | E-fermi (eV) | DOS valence/conduction onset (eV) | Gross occupied DOS width (eV) | Conduction TDOS states/f.u. (0-0.5 / 0.5-1.0 eV) | Dominant VBM character | Dominant CBM character | Notes |",
            "| --- | ---: | --- | ---: | --- | --- | --- | --- |",
        ]
    )
    for result in results:
        dos = result["dos_info"]
        lines.append(
            f"| {result['meta']['display']} | {result['doscar'].efermi:.3f} | "
            f"{dos['valence_onset']:.3f} / {dos['conduction_onset']:.3f} | "
            f"{dos['valence_width']:.2f} | {dos['conduction_0_05_states_per_fu']:.2f} / "
            f"{dos['conduction_05_10_states_per_fu']:.2f} | "
            f"{dominant_character(result['pdos_integrals'], 'VBM')} | "
            f"{dominant_character(result['pdos_integrals'], 'CBM')} | "
            f"{'Clean sampled DOS gap' if dos['clean_gap'] else 'Residual DOS in nominal gap'}; "
            f"{site_note(result)} |"
        )
    lines.extend(
        [
            "",
            (
                "The DOS Fermi/reference energies agree between OUTCAR, vasprun.xml, and DOSCAR "
                "within 1 meV. Integrated DOS at the Fermi level recovers the actual electron count, "
                "supporting the occupied-band assignment `NELECT / 2` for these non-spin-polarized runs."
            ),
            "",
            (
                "The narrower 0.20 eV integrations preserve the S-p-dominated VBM assignment. "
                "They also preserve the mixed S-p/In-s/S-s layered CBM and In-s/S-s/Zn-s spinel "
                "CBM descriptions, so the qualitative edge-character conclusions are not artifacts "
                "of choosing the 0.50 eV reporting window."
            ),
            "",
            "## Cross-Phase Electronic Trends",
            "",
        ]
    )
    ordered = sorted(results, key=lambda result: float(result["uniform_edge"]["fundamental_gap"]))
    ordering_text = " < ".join(
        f"{result['meta']['display']} ({result['uniform_edge']['fundamental_gap']:.3f} eV)"
        for result in ordered
    )
    lines.append(f"The sampled PBE gap ordering is {ordering_text}.")
    lines.append("")
    if len(set(vbm_top)) == 1:
        lines.append(f"VBM character is broadly conserved: {vbm_top[0]} is dominant in all five phases.")
    else:
        lines.append("VBM character is similar but not identical across the phase set: " + ", ".join(vbm_top) + ".")
    lines.append("")
    lines.append(
        "The layered CBMs consistently mix S-p, In-s, and S-s weight; S-p and In-s are close "
        "enough that the leading label should not be interpreted as a pure orbital assignment. "
        "Spinel is distinct in retaining a much larger Zn-s contribution near the CBM."
    )
    lines.append("")
    phase_lookup = {str(result["phase"]): result for result in results}
    beta_gap = float(phase_lookup["beta"]["uniform_edge"]["fundamental_gap"])
    iib_gap = float(phase_lookup["IIb"]["uniform_edge"]["fundamental_gap"])
    alpha_gap = float(phase_lookup["alpha1"]["uniform_edge"]["fundamental_gap"])
    lines.append(
        f"Among the layered stacking variants, beta and IIb differ by {abs(beta_gap - iib_gap):.3f} eV, "
        f"while alpha1 differs from beta by {abs(alpha_gap - beta_gap):.3f} eV. Similar dominant edge "
        "channels support a common local electronic motif, but the nonzero gap and dispersion changes "
        "show that stacking is not electronically invisible."
    )
    lines.append("")
    lines.append(
        "IIa-prime does not introduce a new leading edge-orbital species, but it is distinct through "
        "its off-GAMMA sampled VBM, 42.0 meV direct-minus-fundamental separation, only 0.170 eV "
        "bottom-conduction-band range on GAMMA-A, and concentration of 34.6% of the 0.50 eV VBM "
        "projection on S8. These are topology-associated signatures, not proof of a single causal mechanism."
    )
    lines.append("")
    densities = [float(result["density"]) for result in results]
    correlation = pearson(densities, gaps)
    layered_results = [result for result in results if result["phase"] != "spinel"]
    layered_correlation = pearson(
        [float(result["density"]) for result in layered_results],
        [float(result["uniform_edge"]["fundamental_gap"]) for result in layered_results],
    )
    lines.append(
        f"The all-phase gap-density Pearson coefficient is {correlation:.2f}, but it is dominated by "
        f"the dense, wide-gap spinel. Within the four layered phases it is {layered_correlation:.2f}, "
        "so density does not provide a robust standalone explanation of their gap ordering. The sample "
        "is too small for a causal claim."
    )
    lines.extend(
        [
            "",
            "Spinel is retained as a qualitatively distinct three-dimensional reference. IIa-prime is "
            "treated separately because its reconstructed layer topology can change site weighting and "
            "dispersion even where the leading element/orbital labels remain the same.",
            "",
            "## Literature Context",
            "",
            (
                "Lee et al. focused the accurate electronic-structure comparison of the layered "
                "polytypes on hybrid-functional results and proposed the revised IIa-prime structure. "
                "Accordingly, no direct numerical literature match is claimed for these present PBE gaps. "
                "The Stage 04 results are used to define hypotheses for HSE06, not as substitutes for the "
                "literature HSE06 values. See [Lee et al., Chemistry of Materials 31, 9148-9155 (2019)]"
                "(https://doi.org/10.1021/acs.chemmater.9b03539)."
            ),
            "",
            "## Direct vs Indirect Behaviour",
            "",
            (
                "The primary fundamental gap is `min_k E_CBM(k) - max_k E_VBM(k)` on the dense uniform "
                "mesh. The minimum sampled direct gap is `min_k [E_CBM(k) - E_VBM(k)]`. Their difference "
                "is reported explicitly. High-symmetry path extrema are analysed separately and never "
                "substituted for the full-mesh sampled extrema."
            ),
            "",
            (
                "These are finite-k sampled assignments, not mathematical proofs of the continuous-zone "
                "global extrema. Symmetry reduction means each reported uniform coordinate is an "
                "irreducible representative; equivalent full-zone points are implicit. HSE06 should "
                "recheck the important direct/indirect assignments rather than assume a rigid scissor shift."
            ),
            "",
            "## Stage 05 HSE06 Validation Targets",
            "",
            "| Phase | PBE gap (eV) | PBE sampled classification | Sampled VBM | Sampled CBM | Dominant VBM | Dominant CBM | Main HSE06 question |",
            "| --- | ---: | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for result in results:
        phase = str(result["phase"])
        classification = str(result["uniform_edge"]["classification"])
        if phase == "spinel":
            question = "Does the three-dimensional reference remain electronically distinct and retain its extremum assignment?"
        elif phase == "IIa_prime":
            question = "Does the revised layer topology retain its gap, edge locations, and site-weighting signature?"
        elif classification == "AMBIGUOUS / NEAR-DEGENERATE":
            question = "Does HSE06 resolve the near-degenerate direct/indirect assignment without moving an extremum off the sampled point?"
        else:
            question = "Does HSE06 preserve the sampled direct/indirect class and layered-phase gap ordering?"
        lines.append(
            f"| {result['meta']['display']} | {result['uniform_edge']['fundamental_gap']:.3f} | "
            f"{classification} | {compact_uniform_location(result, 'VBM')} | "
            f"{compact_uniform_location(result, 'CBM')} | "
            f"{dominant_character(result['pdos_integrals'], 'VBM', limit=1)} | "
            f"{dominant_character(result['pdos_integrals'], 'CBM', limit=3)} | {question} |"
        )
    lines.extend(
        [
            "",
            "Across the phase set, HSE06 should test gap ordering, sampled direct/indirect character, "
            "edge-location stability, and whether the leading orbital-projection picture survives. "
            "SOC remains deferred; Stage 04 supplies no specific result that changes that project decision.",
            "",
            "## Acceptance Decision",
            "",
            "| Phase | Decision | Basis |",
            "| --- | --- | --- |",
        ]
    )
    for result in results:
        lines.append(
            f"| {result['meta']['display']} | {result['acceptance']} | Complete, converged, frozen, "
            "lineage-consistent, parseable, semiconducting, and internally cross-checked. |"
        )
    lines.extend(
        [
            "",
            "**Stage 04 freeze decision: PASS.** All five phases can be frozen for the stated PBE scope.",
            "",
            "## Recommended Next Step",
            "",
            "Proceed to the already planned Stage 05 HSE06 validation. No new PBE calculation is required "
            "before HSE06. Preserve the sampled-language and fixed-density provenance when comparing the "
            "two stages, and do not interpret HSE06 as a guaranteed rigid correction to every PBE band edge.",
            "",
        ]
    )
    REPORT_MD.write_text("\n".join(lines))


def print_summary(results: Sequence[dict[str, object]]) -> None:
    print("Validated 10/10 Stage 04 calculations")
    for result in results:
        edge = result["uniform_edge"]
        print(
            f"{result['phase']}: gap={edge['fundamental_gap']:.6f} eV; "
            f"direct={edge['minimum_direct_gap']:.6f} eV; "
            f"delta={1000.0 * edge['delta_direct']:.3f} meV; "
            f"class={edge['classification']}"
        )
        print(f"  VBM {compact_uniform_location(result, 'VBM')} | {dominant_character(result['pdos_integrals'], 'VBM')}")
        print(f"  CBM {compact_uniform_location(result, 'CBM')} | {dominant_character(result['pdos_integrals'], 'CBM')}")
        print(f"  path gap={result['path_edge']['fundamental_gap']:.6f} eV")
        print(f"  {dispersion_sentence(result)}")
        print(f"  {site_note(result)}")
    print(f"Wrote {REPORT_MD}")
    print(f"Wrote {SUMMARY_CSV}")
    print(f"Wrote {EDGES_CSV}")
    print(f"Wrote {PDOS_CSV}")


def main() -> None:
    if not CALC_ROOT.is_dir() or not RESULT_ROOT.is_dir():
        fail("Expected Stage 04 calculation/result directories do not exist")
    results = [analyse_phase(phase) for phase in PHASE_ORDER]
    write_summary_csv(results)
    write_band_edges_csv(results)
    write_pdos_csv(results)
    write_report(results)
    print_summary(results)


if __name__ == "__main__":
    main()
