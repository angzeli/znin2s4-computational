#!/usr/bin/env python3
"""Reconstruct ordered alpha1-ZnIn2S4 from the disordered IIIa COD CIF.

The reconstruction is occupancy-only: it preserves the experimental parent
lattice and fractional coordinates, exhaustively enumerates all minimum-cell
Zn/In assignments, removes parent-symmetry-equivalent assignments, and selects
alpha1 only when exactly one ordering orbit is identified as R3m by spglib.

Runtime dependencies: NumPy and spglib.
"""

from __future__ import annotations

import itertools
import math
import re
import shlex
from collections import Counter
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
    import spglib
except ImportError as exc:  # pragma: no cover - dependency failure path
    raise SystemExit(
        "generate_alpha1.py requires NumPy and spglib; install them in an "
        "isolated environment before running this script"
    ) from exc


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_CIF = PROJECT_ROOT / "structure" / "raw" / "IIIa_cod_1525846.cif"
OUTPUT_CIF = (
    PROJECT_ROOT / "structure" / "raw" / "alpha1_reconstructed_from_IIIa.cif"
)

POSITION_TOLERANCE = 5.0e-5
SYMMETRY_TOLERANCES_ANGSTROM = (1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2)
ATOMIC_NUMBERS = {"S": 16, "Zn": 30, "In": 49}
ELEMENTS_BY_NUMBER = {value: key for key, value in ATOMIC_NUMBERS.items()}


class ReconstructionError(RuntimeError):
    """Raised when the parent or an enumerated structure fails validation."""


def parse_number(token: str) -> float:
    """Parse a CIF number, dropping a trailing standard uncertainty."""

    return float(re.sub(r"\([^)]*\)$", "", token))


def scalar_value(lines: list[str], *keys: str) -> str:
    for raw_line in lines:
        fields = shlex.split(raw_line.strip())
        if fields and fields[0] in keys and len(fields) > 1:
            return fields[1]
    raise ReconstructionError(f"Missing required CIF scalar: {' or '.join(keys)}")


def loop_rows(lines: list[str], target_header: str) -> list[dict[str, str]]:
    """Read one simple whitespace-delimited CIF loop."""

    for index, raw_line in enumerate(lines):
        if raw_line.strip() != "loop_":
            continue
        cursor = index + 1
        headers: list[str] = []
        while cursor < len(lines) and lines[cursor].strip().startswith("_"):
            headers.append(lines[cursor].strip())
            cursor += 1
        if target_header not in headers:
            continue

        rows: list[dict[str, str]] = []
        while cursor < len(lines):
            stripped = lines[cursor].strip()
            if (
                not stripped
                or stripped == "loop_"
                or stripped.startswith("_")
                or stripped.startswith(";")
            ):
                break
            fields = shlex.split(stripped)
            if len(fields) < len(headers):
                raise ReconstructionError(
                    f"Malformed CIF loop row near line {cursor + 1}"
                )
            rows.append(dict(zip(headers, fields, strict=True)))
            cursor += 1
        return rows
    raise ReconstructionError(f"Missing required CIF loop: {target_header}")


def element_from_type_symbol(type_symbol: str) -> str:
    match = re.match(r"[A-Z][a-z]?", type_symbol)
    if match is None:
        raise ReconstructionError(f"Cannot parse element from {type_symbol!r}")
    return match.group(0)


def positions_close(first: np.ndarray, second: np.ndarray) -> bool:
    delta = np.abs(first - second) % 1.0
    delta = np.minimum(delta, 1.0 - delta)
    return bool(np.all(delta < POSITION_TOLERANCE))


def apply_cif_operation(operation: str, position: np.ndarray) -> np.ndarray:
    x, y, z = (float(value) for value in position)
    namespace = {"x": x, "y": y, "z": z}
    values = [
        eval(expression, {"__builtins__": {}}, namespace)  # noqa: S307
        for expression in operation.split(",")
    ]
    return np.mod(np.asarray(values, dtype=float), 1.0)


def lattice_vectors(cell: dict[str, float]) -> np.ndarray:
    a, b, c = (cell[name] for name in ("a", "b", "c"))
    alpha, beta, gamma = (
        math.radians(cell[name]) for name in ("alpha", "beta", "gamma")
    )
    vector_a = np.array([a, 0.0, 0.0])
    vector_b = np.array([b * math.cos(gamma), b * math.sin(gamma), 0.0])
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(
        gamma
    )
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return np.array([vector_a, vector_b, [cx, cy, cz]], dtype=float)


def parse_parent(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise ReconstructionError(f"Parent CIF does not exist: {path}")

    lines = path.read_text(encoding="utf-8").splitlines()
    cell = {
        "a": parse_number(scalar_value(lines, "_cell_length_a")),
        "b": parse_number(scalar_value(lines, "_cell_length_b")),
        "c": parse_number(scalar_value(lines, "_cell_length_c")),
        "alpha": parse_number(scalar_value(lines, "_cell_angle_alpha")),
        "beta": parse_number(scalar_value(lines, "_cell_angle_beta")),
        "gamma": parse_number(scalar_value(lines, "_cell_angle_gamma")),
    }
    operations = [
        row["_symmetry_equiv_pos_as_xyz"]
        for row in loop_rows(lines, "_symmetry_equiv_pos_as_xyz")
    ]
    atom_rows = loop_rows(lines, "_atom_site_label")

    asymmetric_sites: list[dict[str, object]] = []
    expanded_sites: list[dict[str, object]] = []
    for row in atom_rows:
        element = element_from_type_symbol(row["_atom_site_type_symbol"])
        occupancy = parse_number(row["_atom_site_occupancy"])
        position = np.array(
            [
                parse_number(row["_atom_site_fract_x"]),
                parse_number(row["_atom_site_fract_y"]),
                parse_number(row["_atom_site_fract_z"]),
            ],
            dtype=float,
        )
        label = row["_atom_site_label"]
        unique_positions: list[np.ndarray] = []
        for operation in operations:
            transformed = apply_cif_operation(operation, position)
            if not any(
                positions_close(transformed, previous)
                for previous in unique_positions
            ):
                unique_positions.append(transformed)

        asymmetric_sites.append(
            {
                "label": label,
                "element": element,
                "occupancy": occupancy,
                "position": np.mod(position, 1.0),
                "multiplicity": len(unique_positions),
            }
        )
        for transformed in unique_positions:
            existing = next(
                (
                    site
                    for site in expanded_sites
                    if positions_close(transformed, site["position"])
                ),
                None,
            )
            if existing is None:
                expanded_sites.append(
                    {
                        "position": transformed,
                        "composition": {element: occupancy},
                        "source_labels": {label},
                    }
                )
            else:
                composition = existing["composition"]
                if element in composition:
                    raise ReconstructionError(
                        f"Duplicate expanded {element} site near {transformed}"
                    )
                composition[element] = occupancy
                existing["source_labels"].add(label)

    return {
        "formula": scalar_value(lines, "_chemical_formula_sum"),
        "space_group_number": int(
            parse_number(scalar_value(lines, "_space_group_IT_number"))
        ),
        "space_group_name": scalar_value(
            lines, "_symmetry_space_group_name_H-M"
        ),
        "cell": cell,
        "lattice": lattice_vectors(cell),
        "operations": operations,
        "asymmetric_sites": asymmetric_sites,
        "expanded_sites": expanded_sites,
    }


def weighted_composition(sites: list[dict[str, object]]) -> Counter[str]:
    composition: Counter[str] = Counter()
    for site in sites:
        for element, occupancy in site["composition"].items():
            composition[element] += occupancy
    return composition


def validate_parent(parent: dict[str, object]) -> list[dict[str, object]]:
    if parent["space_group_number"] != 166:
        raise ReconstructionError(
            f"Expected IIIa parent space group 166, found "
            f"{parent['space_group_number']}"
        )
    if str(parent["space_group_name"]).replace(" ", "") not in {
        "R-3m:H",
        "R-3m",
    }:
        raise ReconstructionError(
            f"Unexpected IIIa space-group metadata: {parent['space_group_name']}"
        )

    cell = parent["cell"]
    expected_cell = {
        "a": 3.8728,
        "b": 3.8728,
        "c": 37.0664,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 120.0,
    }
    for key, expected in expected_cell.items():
        if not math.isclose(cell[key], expected, abs_tol=1.0e-4):
            raise ReconstructionError(
                f"Unexpected parent {key}: {cell[key]} instead of {expected}"
            )

    sites = parent["expanded_sites"]
    if len(sites) != 21:
        raise ReconstructionError(
            f"Expected 21 expanded parent sites, found {len(sites)}"
        )
    composition = weighted_composition(sites)
    if composition != Counter({"S": 12.0, "In": 6.0, "Zn": 3.0}):
        raise ReconstructionError(f"Unexpected parent composition: {composition}")

    mixed_sites = [
        site for site in sites if set(site["composition"]) == {"In", "Zn"}
    ]
    if len(mixed_sites) != 6:
        raise ReconstructionError(
            f"Expected six mixed Zn/In sites, found {len(mixed_sites)}"
        )
    for site in mixed_sites:
        if not math.isclose(site["composition"]["Zn"], 0.5, abs_tol=1.0e-12):
            raise ReconstructionError("Mixed-site Zn occupancy is not 0.5")
        if not math.isclose(site["composition"]["In"], 0.5, abs_tol=1.0e-12):
            raise ReconstructionError("Mixed-site In occupancy is not 0.5")
        if not math.isclose(sum(site["composition"].values()), 1.0):
            raise ReconstructionError("A mixed site has a vacancy or over-occupancy")

    mixed_site_ids = {id(site) for site in mixed_sites}
    unexpected_partial = [
        site
        for site in sites
        if id(site) not in mixed_site_ids
        and not math.isclose(sum(site["composition"].values()), 1.0)
    ]
    if unexpected_partial:
        raise ReconstructionError("Unexpected partial occupancy outside mixed sites")

    zinc_count = sum(site["composition"]["Zn"] for site in mixed_sites)
    indium_count = sum(site["composition"]["In"] for site in mixed_sites)
    if not (zinc_count.is_integer() and indium_count.is_integer()):
        raise ReconstructionError("Parent cell cannot support integer Zn/In ordering")
    if (int(zinc_count), int(indium_count)) != (3, 3):
        raise ReconstructionError(
            f"Expected three Zn and three In substitutions, found "
            f"{zinc_count:g} and {indium_count:g}"
        )

    return sorted(
        mixed_sites,
        key=lambda site: tuple(float(value) for value in site["position"][[2, 0, 1]]),
    )


def mixed_site_permutations(
    mixed_sites: list[dict[str, object]], operations: list[str]
) -> list[tuple[int, ...]]:
    permutations: list[tuple[int, ...]] = []
    for operation in operations:
        permutation: list[int] = []
        for site in mixed_sites:
            transformed = apply_cif_operation(operation, site["position"])
            matches = [
                index
                for index, target in enumerate(mixed_sites)
                if positions_close(transformed, target["position"])
            ]
            if len(matches) != 1:
                raise ReconstructionError(
                    "A parent symmetry operation does not uniquely permute mixed sites"
                )
            permutation.append(matches[0])
        permutation_tuple = tuple(permutation)
        if permutation_tuple not in permutations:
            permutations.append(permutation_tuple)

    identity = tuple(range(len(mixed_sites)))
    if identity not in permutations:
        raise ReconstructionError("Parent permutation group lacks the identity")
    for first in permutations:
        for second in permutations:
            composed = tuple(first[second[index]] for index in identity)
            if composed not in permutations:
                raise ReconstructionError("Parent mixed-site permutations are not closed")
    return permutations


def enumerate_ordering_orbits(
    site_count: int,
    zinc_count: int,
    permutations: list[tuple[int, ...]],
) -> tuple[list[frozenset[int]], list[dict[str, object]]]:
    raw_combinations = [
        frozenset(indices)
        for indices in itertools.combinations(range(site_count), zinc_count)
    ]
    unseen = set(raw_combinations)
    orbits: list[dict[str, object]] = []
    while unseen:
        seed = min(unseen, key=lambda item: tuple(sorted(item)))
        members = {
            frozenset(permutation[index] for index in seed)
            for permutation in permutations
        }
        representative = min(members, key=lambda item: tuple(sorted(item)))
        orbits.append(
            {
                "representative": representative,
                "members": members,
                "degeneracy": len(members),
            }
        )
        unseen -= members

    covered = set().union(*(orbit["members"] for orbit in orbits))
    if covered != set(raw_combinations):
        raise ReconstructionError("Ordering orbits do not cover all raw combinations")
    orbits.sort(key=lambda orbit: tuple(sorted(orbit["representative"])))
    return raw_combinations, orbits


def build_candidate(
    parent: dict[str, object],
    mixed_sites: list[dict[str, object]],
    zinc_indices: frozenset[int],
) -> dict[str, object]:
    mixed_site_ids = {id(site) for site in mixed_sites}
    fixed_sites = [
        site
        for site in parent["expanded_sites"]
        if id(site) not in mixed_site_ids
    ]
    candidate_sites: list[dict[str, object]] = []
    for site in fixed_sites:
        if len(site["composition"]) != 1:
            raise ReconstructionError("A supposedly fixed site remains disordered")
        element = next(iter(site["composition"]))
        candidate_sites.append(
            {
                "position": np.array(site["position"], copy=True),
                "element": element,
                "origin": "fixed",
            }
        )
    for index, site in enumerate(mixed_sites):
        candidate_sites.append(
            {
                "position": np.array(site["position"], copy=True),
                "element": "Zn" if index in zinc_indices else "In",
                "origin": "mixed",
            }
        )

    positions = np.array([site["position"] for site in candidate_sites])
    numbers = np.array(
        [ATOMIC_NUMBERS[site["element"]] for site in candidate_sites], dtype=int
    )
    return {
        "sites": candidate_sites,
        "spglib_cell": (parent["lattice"], positions, numbers),
    }


def classify_candidate(candidate: dict[str, object]) -> dict[str, object]:
    datasets = []
    for symmetry_tolerance in SYMMETRY_TOLERANCES_ANGSTROM:
        dataset = spglib.get_symmetry_dataset(
            candidate["spglib_cell"], symprec=symmetry_tolerance
        )
        if dataset is None:
            raise ReconstructionError(
                f"spglib failed at symprec={symmetry_tolerance:g} angstrom"
            )
        datasets.append(dataset)

    fingerprints = {
        (dataset.international.replace(" ", ""), int(dataset.number))
        for dataset in datasets
    }
    if len(fingerprints) != 1:
        raise ReconstructionError(
            f"Candidate symmetry is tolerance-sensitive: {sorted(fingerprints)}"
        )
    name, number = next(iter(fingerprints))
    return {
        "space_group_name": name,
        "space_group_number": number,
        "dataset": datasets[2],
        "fingerprints": fingerprints,
    }


def Cartesian_displacement(
    first: np.ndarray, second: np.ndarray, lattice: np.ndarray
) -> float:
    delta = second - first
    best = math.inf
    for shift in itertools.product((-1, 0, 1), repeat=3):
        vector = (delta + np.asarray(shift)) @ lattice
        best = min(best, float(np.linalg.norm(vector)))
    return best


def minimum_pair_distance(candidate: dict[str, object], lattice: np.ndarray) -> float:
    sites = candidate["sites"]
    distances = [
        Cartesian_displacement(first["position"], second["position"], lattice)
        for index, first in enumerate(sites)
        for second in sites[index + 1 :]
    ]
    return min(distances)


def sulfur_neighbour_count(
    site: dict[str, object],
    candidate: dict[str, object],
    lattice: np.ndarray,
    cutoff: float = 2.90,
) -> int:
    count = 0
    for sulfur in candidate["sites"]:
        if sulfur["element"] != "S":
            continue
        for shift in itertools.product((-1, 0, 1), repeat=3):
            delta = sulfur["position"] + np.asarray(shift) - site["position"]
            if float(np.linalg.norm(delta @ lattice)) < cutoff:
                count += 1
    return count


def matrix_rank(vectors: list[tuple[int, int, int]]) -> int:
    nonzero = [vector for vector in vectors if any(vector)]
    if not nonzero:
        return 0
    return int(np.linalg.matrix_rank(np.asarray(nonzero, dtype=float), tol=1.0e-9))


def periodic_bond_components(
    candidate: dict[str, object], lattice: np.ndarray
) -> list[tuple[int, int]]:
    sites = candidate["sites"]
    edges: list[list[tuple[int, tuple[int, int, int]]]] = [
        [] for _ in sites
    ]
    for cation_index, cation in enumerate(sites):
        if cation["element"] not in {"Zn", "In"}:
            continue
        for sulfur_index, sulfur in enumerate(sites):
            if sulfur["element"] != "S":
                continue
            for shift in itertools.product((-1, 0, 1), repeat=3):
                delta = sulfur["position"] + np.asarray(shift) - cation["position"]
                if float(np.linalg.norm(delta @ lattice)) < 2.90:
                    edges[cation_index].append((sulfur_index, shift))
                    edges[sulfur_index].append(
                        (cation_index, tuple(-value for value in shift))
                    )

    visited: set[int] = set()
    components: list[tuple[int, int]] = []
    for root in range(len(sites)):
        if root in visited:
            continue
        assigned = {root: (0, 0, 0)}
        stack = [root]
        cycles: list[tuple[int, int, int]] = []
        while stack:
            current = stack.pop()
            visited.add(current)
            for target, shift in edges[current]:
                expected = tuple(
                    assigned[current][axis] + shift[axis] for axis in range(3)
                )
                if target not in assigned:
                    assigned[target] = expected
                    stack.append(target)
                else:
                    cycles.append(
                        tuple(
                            expected[axis] - assigned[target][axis]
                            for axis in range(3)
                        )
                    )
        components.append((len(assigned), matrix_rank(cycles)))
    return sorted(components)


def vertical_plane_gaps(
    candidate: dict[str, object], c_length: float
) -> list[float]:
    planes = sorted(
        {
            round(float(site["position"][2] % 1.0), 10)
            for site in candidate["sites"]
        }
    )
    return sorted(
        [
            ((planes[(index + 1) % len(planes)] - planes[index]) % 1.0)
            * c_length
            for index in range(len(planes))
        ],
        reverse=True,
    )


def validate_alpha1(
    candidate: dict[str, object],
    classification: dict[str, object],
    parent: dict[str, object],
) -> dict[str, object]:
    sites = candidate["sites"]
    counts = Counter(site["element"] for site in sites)
    if counts != Counter({"S": 12, "In": 6, "Zn": 3}):
        raise ReconstructionError(f"Incorrect alpha1 composition: {counts}")
    if len(sites) != 21:
        raise ReconstructionError(f"Expected 21 alpha1 atoms, found {len(sites)}")
    if classification["space_group_number"] != 160:
        raise ReconstructionError(
            f"Selected alpha1 is not R3m: {classification['space_group_name']} "
            f"({classification['space_group_number']})"
        )

    minimum_distance = minimum_pair_distance(candidate, parent["lattice"])
    if minimum_distance < 2.0:
        raise ReconstructionError(
            f"Unphysically short alpha1 contact: {minimum_distance:.6f} angstrom"
        )

    mixed_coordination = [
        sulfur_neighbour_count(site, candidate, parent["lattice"])
        for site in sites
        if site["origin"] == "mixed"
    ]
    fixed_indium_coordination = [
        sulfur_neighbour_count(site, candidate, parent["lattice"])
        for site in sites
        if site["origin"] == "fixed" and site["element"] == "In"
    ]
    if mixed_coordination != [4] * 6:
        raise ReconstructionError(
            f"Unexpected tetrahedral-site coordination: {mixed_coordination}"
        )
    if fixed_indium_coordination != [6] * 3:
        raise ReconstructionError(
            f"Unexpected octahedral In coordination: {fixed_indium_coordination}"
        )

    components = periodic_bond_components(candidate, parent["lattice"])
    if components != [(7, 2), (7, 2), (7, 2)]:
        raise ReconstructionError(
            f"Alpha1 no longer has three layered septuple networks: {components}"
        )
    large_gaps = [
        gap
        for gap in vertical_plane_gaps(candidate, parent["cell"]["c"])
        if gap > 2.8
    ]
    if len(large_gaps) != 3:
        raise ReconstructionError(
            f"Expected three van der Waals plane gaps, found {large_gaps}"
        )

    return {
        "composition": counts,
        "minimum_distance": minimum_distance,
        "components": components,
        "large_gaps": large_gaps,
    }


def operation_expression(rotation: np.ndarray, translation: np.ndarray) -> str:
    variables = ("x", "y", "z")
    coordinate_expressions: list[str] = []
    for row, offset in zip(rotation, translation, strict=True):
        expression = ""
        for coefficient, variable in zip(row, variables, strict=True):
            coefficient = int(coefficient)
            if coefficient == 0:
                continue
            if abs(coefficient) != 1:
                raise ReconstructionError(
                    f"Unsupported symmetry coefficient: {coefficient}"
                )
            term = variable if coefficient > 0 else f"-{variable}"
            if expression and coefficient > 0:
                expression += "+"
            expression += term

        wrapped_offset = float(offset) % 1.0
        if math.isclose(wrapped_offset, 1.0, abs_tol=1.0e-8):
            wrapped_offset = 0.0
        fraction = Fraction(wrapped_offset).limit_denominator(24)
        if not math.isclose(float(fraction), wrapped_offset, abs_tol=1.0e-7):
            raise ReconstructionError(
                f"Cannot express symmetry translation exactly: {wrapped_offset}"
            )
        if fraction:
            translation_text = (
                str(fraction.numerator)
                if fraction.denominator == 1
                else f"{fraction.numerator}/{fraction.denominator}"
            )
            if expression:
                expression += "+"
            expression += translation_text
        coordinate_expressions.append(expression or "0")
    return ",".join(coordinate_expressions)


def asymmetric_representatives(
    candidate: dict[str, object], dataset: object
) -> list[dict[str, object]]:
    groups: dict[int, list[int]] = {}
    for atom_index, representative in enumerate(dataset.equivalent_atoms):
        groups.setdefault(int(representative), []).append(atom_index)

    sites = candidate["sites"]
    representatives: list[dict[str, object]] = []
    for group in groups.values():
        representative_index = group[0]
        elements = {sites[index]["element"] for index in group}
        if len(elements) != 1:
            raise ReconstructionError("spglib mixed unlike elements in one orbit")
        representatives.append(
            {
                "element": next(iter(elements)),
                "position": sites[representative_index]["position"],
                "multiplicity": len(group),
                "wyckoff": str(dataset.wyckoffs[representative_index]),
            }
        )

    element_order = {"Zn": 0, "In": 1, "S": 2}
    representatives.sort(
        key=lambda site: (
            element_order[site["element"]],
            float(site["position"][2]),
            float(site["position"][0]),
            float(site["position"][1]),
        )
    )
    return representatives


def validate_asymmetric_model(
    representatives: list[dict[str, object]],
    operations: list[tuple[np.ndarray, np.ndarray]],
    candidate: dict[str, object],
) -> None:
    expanded: list[tuple[str, np.ndarray]] = []
    for representative in representatives:
        unique_positions: list[np.ndarray] = []
        for rotation, translation in operations:
            transformed = np.mod(
                rotation @ representative["position"] + translation, 1.0
            )
            if not any(
                positions_close(transformed, previous)
                for previous in unique_positions
            ):
                unique_positions.append(transformed)
        if len(unique_positions) != representative["multiplicity"]:
            raise ReconstructionError(
                "Serialized asymmetric-site multiplicity does not match spglib"
            )
        expanded.extend(
            (representative["element"], position) for position in unique_positions
        )

    if len(expanded) != len(candidate["sites"]):
        raise ReconstructionError("Serialized CIF would expand to the wrong atom count")
    for site in candidate["sites"]:
        matches = [
            index
            for index, (element, position) in enumerate(expanded)
            if element == site["element"]
            and positions_close(position, site["position"])
        ]
        if len(matches) != 1:
            raise ReconstructionError(
                "Serialized CIF does not reproduce the selected alpha1 candidate"
            )


def cif_text(
    parent: dict[str, object],
    candidate: dict[str, object],
    classification: dict[str, object],
    raw_count: int,
    distinct_count: int,
) -> str:
    dataset = classification["dataset"]
    operations = list(zip(dataset.rotations, dataset.translations, strict=True))
    representatives = asymmetric_representatives(candidate, dataset)
    validate_asymmetric_model(representatives, operations, candidate)

    cell = parent["cell"]
    volume = abs(float(np.linalg.det(parent["lattice"])))
    lines = [
        "data_alpha1_ZnIn2S4_reconstructed",
        "_audit_creation_method",
        ";",
        "Ordered alpha1-ZnIn2S4 reconstructed from experimental disordered",
        "IIIa-ZnIn2S4 parent COD 1525846. The experimental lattice and all",
        "fractional coordinates were preserved; only Zn/In occupations were",
        f"resolved by exhaustive symmetry-aware enumeration of {raw_count} raw",
        f"minimum-cell assignments ({distinct_count} parent-symmetry orbits).",
        "The unique R3m ordering orbit was selected using spglib. This is not",
        "an experimentally refined alpha1 CIF and no DFT relaxation was applied.",
        ";",
        "_chemical_formula_structural       'ZnIn2S4'",
        "_chemical_formula_sum              'In6 S12 Zn3'",
        "_cell_formula_units_Z              3",
        "_space_group_name_H-M_alt          'R 3 m :H'",
        "_space_group_IT_number             160",
        f"_cell_length_a                     {cell['a']:.6f}",
        f"_cell_length_b                     {cell['b']:.6f}",
        f"_cell_length_c                     {cell['c']:.6f}",
        f"_cell_angle_alpha                  {cell['alpha']:.6f}",
        f"_cell_angle_beta                   {cell['beta']:.6f}",
        f"_cell_angle_gamma                  {cell['gamma']:.6f}",
        f"_cell_volume                       {volume:.6f}",
        "",
        "loop_",
        "_space_group_symop_id",
        "_space_group_symop_operation_xyz",
    ]
    for operation_id, (rotation, translation) in enumerate(operations, start=1):
        expression = operation_expression(rotation, translation)
        lines.append(f"{operation_id} '{expression}'")

    lines.extend(
        [
            "",
            "loop_",
            "_atom_site_label",
            "_atom_site_type_symbol",
            "_atom_site_symmetry_multiplicity",
            "_atom_site_Wyckoff_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_occupancy",
        ]
    )
    label_counts: Counter[str] = Counter()
    for representative in representatives:
        element = representative["element"]
        label_counts[element] += 1
        label = f"{element}{label_counts[element]}"
        x, y, z = representative["position"]
        lines.append(
            f"{label:<4} {element:<2} {representative['multiplicity']:>2} "
            f"{representative['wyckoff']} "
            f"{x:.10f} {y:.10f} {z:.10f} 1.000000"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    if OUTPUT_CIF.exists():
        raise ReconstructionError(f"Refusing to overwrite existing file: {OUTPUT_CIF}")

    parent = parse_parent(INPUT_CIF)
    mixed_sites = validate_parent(parent)
    permutations = mixed_site_permutations(mixed_sites, parent["operations"])
    raw_combinations, orbits = enumerate_ordering_orbits(
        site_count=len(mixed_sites),
        zinc_count=3,
        permutations=permutations,
    )

    evaluated: list[dict[str, object]] = []
    for orbit_index, orbit in enumerate(orbits, start=1):
        candidate = build_candidate(
            parent, mixed_sites, orbit["representative"]
        )
        classification = classify_candidate(candidate)
        evaluated.append(
            {
                "index": orbit_index,
                "orbit": orbit,
                "candidate": candidate,
                "classification": classification,
            }
        )

    alpha1_matches = [
        item
        for item in evaluated
        if item["classification"]["space_group_number"] == 160
        and item["classification"]["space_group_name"] == "R3m"
    ]
    if len(alpha1_matches) != 1:
        raise ReconstructionError(
            f"Expected exactly one R3m ordering orbit, found {len(alpha1_matches)}"
        )

    alpha1 = alpha1_matches[0]
    validation = validate_alpha1(
        alpha1["candidate"], alpha1["classification"], parent
    )
    output_text = cif_text(
        parent,
        alpha1["candidate"],
        alpha1["classification"],
        raw_count=len(raw_combinations),
        distinct_count=len(orbits),
    )
    with OUTPUT_CIF.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(output_text)

    print(f"Parent mixed sites: {len(mixed_sites)} x (Zn0.5 In0.5)")
    print(f"Raw occupational combinations: {len(raw_combinations)}")
    print(f"Parent-symmetry-distinct configurations: {len(orbits)}")
    for item in evaluated:
        representative = tuple(sorted(item["orbit"]["representative"]))
        classification = item["classification"]
        print(
            f"  candidate {item['index']}: Zn mixed-site indices "
            f"{representative}; degeneracy {item['orbit']['degeneracy']}; "
            f"{classification['space_group_name']} "
            f"({classification['space_group_number']})"
        )
    print(
        f"alpha1: candidate {alpha1['index']} (unique R3m orbit), "
        f"ZnIn2S4, {len(alpha1['candidate']['sites'])} atoms"
    )
    print(
        f"Minimum contact: {validation['minimum_distance']:.6f} angstrom; "
        f"layered components: {validation['components']}"
    )
    print(f"Wrote: {OUTPUT_CIF}")


if __name__ == "__main__":
    main()
