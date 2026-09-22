"""Geometry-only diagnostics for unrelaxed WP2 seeds; distances are in angstrom.

Bond graphs use phase/species first shells validated against the WP1 audit.
Formal charges and broken bonds are screening diagnostics, never energies or
proof of physical polarity. No calculator, optimization, or job API is used.
"""

from collections import Counter

import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

CHARGES = {"Zn": 2, "In": 3, "S": -2}
# Resolve spinel's distinct planes separated by about 0.04 A; 0.05 A merges them.
PLANE_TOL = 0.01
VACUUM = 20.0
FACE_DEPTH = 3.6


def composition_record(structure):
    """Return exact integer formula-unit bookkeeping, including nonstoichiometry."""
    counts = structure.composition.get_el_amt_dict()
    units = counts.get("Zn", 0)
    stoichiometric = units > 0 and counts == {"Zn": units, "In": 2 * units, "S": 4 * units}
    return {
        "composition": structure.composition.formula,
        "atom_count": len(structure),
        "n_formula_units_if_integer": int(units) if stoichiometric and units.is_integer() else None,
        "stoichiometric": stoichiometric,
        "nominal_charge": sum(CHARGES[s] * n for s, n in counts.items()),
    }


def coordination_reference(parent, phase):
    """Locate each site's first distance gap >=0.5 A after 4 or 6 S neighbors.

    The 0.5 A separation is validated explicitly, not assumed to be a bond
    cutoff. A species-wide cutoff lies halfway between the largest included
    and smallest excluded distance in that phase. Abort on ambiguous shells.
    """
    criteria = {}
    for element in ("Zn", "In"):
        shells = []
        for i, site in enumerate(parent):
            if site.specie.symbol != element:
                continue
            distances = sorted(n.nn_distance for n in parent.get_neighbors(site, 5.0)
                               if n.specie.symbol == "S")
            choices = [n for n in (4, 6) if distances[n] - distances[n - 1] >= 0.5]
            if not choices:
                raise ValueError(f"{phase} {element}{i}: no resolved first shell")
            count = choices[0]
            shells.append(dict(parent_site_index=i, coordination=count,
                               minimum_A=float(distances[0]), maximum_A=float(distances[count - 1]),
                               next_distance_A=float(distances[count])))
        observed = sorted(s["coordination"] for s in shells)
        expected = ([4] * len(shells) if element == "Zn" else
                    [6] * len(shells) if phase == "spinel" else
                    sorted([4, 6] * (len(shells) // 2)))
        if observed != expected:
            raise ValueError(f"{phase} {element}: first shells disagree with WP1")
        last = max(s["maximum_A"] for s in shells)
        following = min(s["next_distance_A"] for s in shells)
        if following - last < 0.5:
            raise ValueError("Species-wide shell boundary is ambiguous")
        criteria[element] = dict(cutoff_A=(last + following) / 2,
                                 included_max_A=last, excluded_min_A=following, sites=shells)
    return criteria


def bond_graph(structure, criteria):
    """Directed cation-S neighbors with periodic image offsets; no z bonds in vacuum."""
    adjacency = [[] for _ in structure]
    for i, site in enumerate(structure):
        symbol = site.specie.symbol
        if symbol not in criteria:
            continue
        for neighbor in structure.get_neighbors(site, criteria[symbol]["cutoff_A"]):
            if neighbor.specie.symbol != "S":
                continue
            j = int(neighbor.index)
            offset = np.rint(neighbor.image).astype(int)
            adjacency[i].append((j, offset))
            adjacency[j].append((i, -offset))
    return adjacency


def connectivity(structure, adjacency):
    """Resolve components and periodic rank from translation cycles of the bond graph.

    A rank-0 component is finite even if it contains several atoms. Rank-1
    edge ribbons and rank-2 septuple sheets can be legitimate vdW-separated
    components; graph disconnection alone is not a rejection criterion.
    """
    seen = set()
    components = []
    for root in range(len(structure)):
        if root in seen:
            continue
        offsets = {root: np.zeros(3, dtype=int)}
        queue = [root]
        cycles = []
        for i in queue:
            seen.add(i)
            for j, image in adjacency[i]:
                proposed = offsets[i] + image
                if j not in offsets:
                    offsets[j] = proposed
                    queue.append(j)
                else:
                    cycles.append(proposed - offsets[j])
        nonzero = [x for x in cycles if np.any(x)]
        rank = int(np.linalg.matrix_rank(np.array(nonzero, dtype=float))) if nonzero else 0
        components.append({"site_indices": sorted(offsets), "atom_count": len(offsets),
                           "composition": dict(Counter(structure[i].specie.symbol for i in offsets)),
                           "periodic_rank": rank,
                           "periodic_cycle_vectors": sorted({tuple(map(int, x)) for x in nonzero})})
    return components


def cut_diagnostics(generator, shift, criteria):
    """Measure the bracketing gap and count bulk cation-S bonds crossing the cut.

    Count each bond from its cation endpoint only, including periodic images.
    Fractional shifts refer to the oriented unit cell, before centering.
    """
    bulk = generator.oriented_unit_cell
    z = bulk.frac_coords[:, 2] % 1
    height = float(abs(np.dot(bulk.lattice.matrix[2], generator._normal)))
    below = float(np.min((shift - z) % 1))
    above = float(np.min((z - shift) % 1))
    crosses = Counter({"Zn": 0, "In": 0})
    for i, neighbors in enumerate(bond_graph(bulk, criteria)):
        element = bulk[i].specie.symbol
        if element not in criteria:
            continue
        for j, image in neighbors:
            lo, hi = sorted((float(bulk.frac_coords[i, 2]), float(bulk.frac_coords[j, 2] + image[2])))
            crosses[element] += sum(lo + 1e-8 < shift + k < hi - 1e-8
                                    for k in range(int(np.floor(lo)) - 1, int(np.ceil(hi)) + 1))
    return {"fractional_shift": float(shift), "nearest_gap_width_A": (below + above) * height,
            "crossing_ZnS": crosses["Zn"], "crossing_InS": crosses["In"],
            "no_first_shell_bond_crosses": not any(crosses.values()),
            "oriented_bulk_planes": plane_sequence(bulk, z * height)}


def present_slab(raw, generator):
    """Rigidly rotate to a||x, normal||z, then replace only the vacuum cell vector.

    Cartesian atomic positions undergo one rigid rotation and translation.
    In-plane periodic wrapping is a representation change, never relaxation.
    """
    a, b, _ = raw.lattice.matrix
    x = a / np.linalg.norm(a)
    z = np.cross(a, b)
    z /= np.linalg.norm(z)
    y = np.cross(z, x)
    rotation = np.column_stack((x, y, z))
    coords = raw.cart_coords @ rotation
    thickness = float(np.ptp(coords[:, 2]))
    shift_z = VACUUM / 2 - float(coords[:, 2].min())
    coords[:, 2] += shift_z
    lattice = np.array([a @ rotation, b @ rotation, [0, 0, thickness + VACUUM]])
    lattice[np.abs(lattice) < 1e-12] = 0
    slab = Structure(Lattice(lattice), raw.species, coords, coords_are_cartesian=True,
                     site_properties=raw.site_properties)
    fractional = slab.frac_coords.copy()
    wraps = np.floor(fractional[:, :2]).astype(int)
    fractional[:, :2] -= wraps
    slab = Structure(slab.lattice, slab.species, fractional, site_properties=slab.site_properties)
    transformation = {
        "oriented_unit_cell_matrix": np.asarray(generator.slab_scale_factor).tolist(),
        "oriented_unit_cell_lattice_A": generator.oriented_unit_cell.lattice.matrix.tolist(),
        "oriented_repeat_count": len(raw) // len(generator.oriented_unit_cell),
        "raw_lattice_A": raw.lattice.matrix.tolist(),
        "row_cartesian_rotation": rotation.tolist(), "translation_z_A": shift_z,
        "inplane_wrap_integers": wraps.tolist(),
        "final_lattice_A": lattice.tolist(),
        "parent_site_indices": list(map(int, slab.site_properties["parent_index"])),
        "literal_parent_basis_preserved": False,
    }
    return slab, transformation


def plane_sequence(structure, heights=None):
    """Sorted normal planes clustered within 0.01 A, retaining composition/charge."""
    if heights is None:
        heights = structure.cart_coords[:, 2]
    groups = []
    for i in np.argsort(heights, kind="stable"):
        if not groups or heights[i] - heights[groups[-1][0]] > PLANE_TOL:
            groups.append([])
        groups[-1].append(int(i))
    origin = min(heights)
    return [dict(z_A=float(np.mean(np.asarray(heights)[g]) - origin),
                 species=dict(sorted(Counter(structure[i].specie.symbol for i in g).items())),
                 formal_charge_e=sum(CHARGES[structure[i].specie.symbol] for i in g),
                 site_indices=g) for g in groups]


def faces_equivalent(structure):
    """Require a whole-slab symmetry operation that exchanges its exposed faces."""
    for operation in SpacegroupAnalyzer(structure, symprec=1e-3).get_symmetry_operations():
        if np.allclose(operation.rotation_matrix[2], [0, 0, -1], atol=1e-7):
            mapped = operation.operate_multi(structure.frac_coords) % 1
            expected = 1 - structure.frac_coords[:, 2]
            if np.allclose(mapped[:, 2], expected, atol=1e-5):
                return True
    return False


def face_records(structure, coordination, deficits):
    """Face fingerprints use outer 3.6 A; all lost-coordination sites use their nearest face.

    The finite-depth formal charge is partition-dependent and is reported only
    as a diagnostic. Both exposed faces always keep separate stable IDs.
    """
    z = structure.cart_coords[:, 2]
    midpoint = (z.min() + z.max()) / 2
    faces = {}
    for side in ("upper", "lower"):
        depth = z.max() - z if side == "upper" else z - z.min()
        indices = np.flatnonzero(depth <= FACE_DEPTH)
        nearest = z >= midpoint if side == "upper" else z < midpoint
        under = Counter(structure[i].specie.symbol for i in range(len(structure))
                        if nearest[i] and deficits[i] > 0)
        planes = plane_sequence(structure, depth)
        layers = [dict(depth_A=round(p["z_A"], 5), species=p["species"],
                       coordination=sorted((structure[i].specie.symbol, int(coordination[i]))
                                           for i in p["site_indices"]))
                  for p in planes if p["z_A"] <= FACE_DEPTH]
        faces[side] = dict(surface_species_fingerprint=layers,
                           outermost_species=planes[0]["species"],
                           face_composition=dict(Counter(structure[i].specie.symbol for i in indices)),
                           undercoord_Zn=under["Zn"], undercoord_In=under["In"], undercoord_S=under["S"],
                           formal_surface_charge_diagnostic=sum(CHARGES[structure[i].specie.symbol]
                                                                for i in indices))
    return faces


def diagnose_slab(slab, parent, criteria, basal):
    """Evaluate inherited coordination, topology, formal charge, and geometric validity."""
    parent_cn = np.array([len(n) for n in bond_graph(parent, criteria)])
    graph = bond_graph(slab, criteria)
    cn = np.array([len(n) for n in graph])
    indices = slab.site_properties["parent_index"]
    deficits = parent_cn[indices] - cn
    components = connectivity(slab, graph)
    z = slab.cart_coords[:, 2]
    thickness = float(np.ptp(z))
    height = float(slab.lattice.matrix[2, 2])
    neighbors = slab.get_neighbor_list(5)
    minimum = float(min(neighbors[3]))
    area = float(np.linalg.norm(np.cross(*slab.lattice.matrix[:2])))
    equivalent = faces_equivalent(slab)
    charge = np.array([CHARGES[s.specie.symbol] for s in slab])
    dipole = float(np.dot(charge, z - (z.min() + z.max()) / 2))
    density = dipole / area
    comp = composition_record(slab)
    # Transparent triage thresholds, not physical polar/nonpolar assignments.
    if not comp["stoichiometric"] or comp["nominal_charge"] != 0 or abs(density) > 0.1:
        polarity = "HIGH_CONCERN"
    elif not equivalent or abs(density) > 1e-3:
        polarity = "WARNING"
    else:
        polarity = "LOW_CONCERN"
    errors = []
    if not np.isfinite(slab.lattice.matrix).all() or slab.volume <= 0 or not np.isfinite(slab.cart_coords).all():
        errors.append("nonfinite or invalid cell/coordinates")
    parent_min = float(min(parent.get_neighbor_list(5)[3]))
    if minimum < 1.8 or minimum < parent_min - 1e-5:
        errors.append("overlap or contact shorter than inherited parent minimum")
    if abs(height - thickness - VACUUM) > 1e-8 or abs(z.min() - VACUUM / 2) > 1e-8:
        errors.append("incorrect vacuum or centering")
    if not np.allclose(slab.lattice.matrix[:2, 2], 0, atol=1e-10) or not np.allclose(slab.lattice.matrix[2, :2], 0):
        errors.append("incorrect normal convention")
    if any(image[2] != 0 for nn in graph for _, image in nn):
        errors.append("unexpected bond across vacuum")
    if np.any(deficits < 0):
        errors.append("coordination exceeds inherited parent")
    finite = [c for c in components if c["periodic_rank"] == 0]
    if finite:
        errors.append("detached finite fragment(s)")
    complete = (all(c["composition"] == {"Zn": 1, "In": 2, "S": 4}
                    and c["periodic_rank"] == 2 for c in components) and not np.any(deficits))
    if basal and not complete:
        errors.append("incomplete or bond-cut basal layered unit")
    broken = {e: int(sum(deficits[i] for i, site in enumerate(slab) if site.specie.symbol == e))
              for e in ("Zn", "In")}
    return dict(**comp, cell_normal_length_A=height, slab_thickness_A=thickness,
                vacuum_A=height - thickness, surface_area_A2=area, minimum_distance_A=minimum,
                faces_equivalent=equivalent, faces=face_records(slab, cn, deficits),
                coordination=dict(parent_counts=parent_cn.tolist(), slab_counts=cn.tolist(),
                                  deficits=deficits.tolist(), broken_bonds=broken),
                components=components,
                connectivity_status="FINITE_FRAGMENTS" if finite else
                "CONNECTED" if len(components) == 1 else "PERIODIC_VDW_COMPONENTS",
                complete_layer_count=len(components) if basal and complete else None,
                polarity_screen=polarity, formal_dipole_eA=dipole, formal_dipole_per_area_e_per_A=density,
                charged_plane_sequence=plane_sequence(slab), validation_errors=errors)


def slabs_match(first, second, first_diagnostics, second_diagnostics):
    """Strict composition/lattice/layer/environment prefilter, then periodic structure matching.

    Slab-pair matching permits normal reversal. It never determines whether
    the upper and lower faces of one slab are chemically equivalent.
    """
    if first.composition != second.composition:
        return False
    if not np.allclose(first.lattice.matrix[:2] @ first.lattice.matrix[:2].T,
                       second.lattice.matrix[:2] @ second.lattice.matrix[:2].T, atol=1e-5, rtol=0):
        return False
    if abs(first_diagnostics["slab_thickness_A"] - second_diagnostics["slab_thickness_A"]) > 1e-4:
        return False
    def environment(structure, data):
        return sorted((s.specie.symbol, cn) for s, cn in zip(structure, data["coordination"]["slab_counts"]))
    if environment(first, first_diagnostics) != environment(second, second_diagnostics):
        return False
    a, b = plane_sequence(first), plane_sequence(second)
    def layers_match(left, right):
        return (len(left) == len(right) and all(x["species"] == y["species"] for x, y in zip(left, right))
                and np.allclose([x["z_A"] for x in left], [x["z_A"] for x in right], atol=1e-4, rtol=0))
    reverse = [dict(p, z_A=second_diagnostics["slab_thickness_A"] - p["z_A"]) for p in reversed(b)]
    if not (layers_match(a, b) or layers_match(a, reverse)):
        return False
    matcher = StructureMatcher(ltol=1e-5, stol=0.003, angle_tol=0.01,
                               primitive_cell=False, scale=False, attempt_supercell=False)
    return bool(matcher.fit(first, second, skip_structure_reduction=True))
