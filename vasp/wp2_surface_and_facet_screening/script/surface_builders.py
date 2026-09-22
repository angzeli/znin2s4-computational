"""Pure vacuum and complete-repeat slab construction from Stage 0 provenance.

No files or jobs are created. Thickness construction makes a new bulk-derived
seed; it refuses a relaxed or otherwise changed input instead of silently
discarding its displacements. Scientific approval is the caller's separate gate.
"""

import copy
import json
import math

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator

from surface_model_utils import coordination_reference, diagnose_slab, present_slab


def slab_geometry(slab: Structure) -> dict:
    """Return atom-envelope geometry for an unwrapped, normal-aligned slab.

    The atom envelope must lie inside one cell with its uniquely largest normal
    gap across the cell boundary. A split or ambiguous envelope is rejected;
    unwrapping needs an explicit decision. In-plane vectors may be oblique.
    """
    lattice = np.asarray(slab.lattice.matrix, dtype=float)
    coords = np.asarray(slab.cart_coords, dtype=float)
    if not len(slab) or not np.isfinite(lattice).all() or not np.isfinite(coords).all():
        raise ValueError("A finite, nonempty slab is required")
    normal = np.cross(lattice[0], lattice[1])
    area = float(np.linalg.norm(normal))
    if area <= 0:
        raise ValueError("The in-plane cell is degenerate")
    normal /= area
    height = float(np.dot(lattice[2], normal))
    if height <= 0 or not np.allclose(lattice[2], height * normal, atol=1e-9, rtol=0):
        raise ValueError("The third lattice vector must follow the slab normal")
    z = coords @ normal
    lower, upper = float(z.min()), float(z.max())
    if lower < -1e-9 or upper > height + 1e-9 or height - (upper - lower) <= 0:
        raise ValueError("The slab must be unwrapped inside a cell with atom-free vacuum")
    vacuum = height - (upper - lower)
    if len(z) > 1 and vacuum <= float(np.max(np.diff(np.sort(z)))) + 1e-9:
        raise ValueError("The boundary vacuum is not uniquely resolved; unwrap or review the slab envelope")
    return dict(normal=normal.tolist(), heights_A=z.tolist(), lower_A=lower,
                upper_A=upper, slab_thickness_A=upper - lower,
                cell_normal_length_A=height, vacuum_A=vacuum,
                surface_area_A2=area)


def change_vacuum(slab: Structure, vacuum_A: float, *, face_context=None):
    """Return a centered copy with the requested actual atom-free separation (A).

    Only the normal cell vector and a common Cartesian translation change.
    In-plane vectors, species ordering, site properties, and internal Cartesian
    pair vectors are retained to floating-point precision, without wrapping.
    """
    if not np.isfinite(vacuum_A) or vacuum_A <= 0:
        raise ValueError("vacuum_A must be finite and positive")
    geometry = slab_geometry(slab)
    normal = np.asarray(geometry["normal"])
    lattice = slab.lattice.matrix.copy()
    lattice[2] = normal * (geometry["slab_thickness_A"] + vacuum_A)
    translation = normal * (vacuum_A / 2 - geometry["lower_A"])
    result = Structure(Lattice(lattice), slab.species, slab.cart_coords + translation,
                       coords_are_cartesian=True, site_properties=copy.deepcopy(slab.site_properties))
    metadata = dict(builder="rigid_vacuum_change", target_vacuum_A=float(vacuum_A),
                    **slab_geometry(result), translation_cartesian_A=translation.tolist(),
                    face_context=copy.deepcopy(face_context), internal_geometry_preserved=True)
    return result, metadata


def _regenerate(parent, metadata, repeats):
    """Build exactly the specified number of the recorded oriented unit cell."""
    settings = copy.deepcopy(metadata["construction_settings"])
    recorded = np.asarray(metadata["transformation"]["oriented_unit_cell_lattice_A"])
    normal = np.cross(recorded[0], recorded[1])
    normal /= np.linalg.norm(normal)
    height = abs(float(np.dot(recorded[2], normal)))
    settings["min_slab_size"] = height * (repeats - 0.5)
    generator = SlabGenerator(parent, tuple(metadata["hkl"]), **settings)
    if not np.allclose(generator.oriented_unit_cell.lattice.matrix, recorded, atol=1e-8, rtol=0):
        raise ValueError("Supplied parent does not reproduce the recorded oriented lattice")
    if not np.array_equal(generator.slab_scale_factor,
                          metadata["transformation"]["oriented_unit_cell_matrix"]):
        raise ValueError("Oriented-cell transformation changed")
    raw = generator.get_slab(shift=float(metadata["raw_cut_shift"]))
    result, transform = present_slab(raw, generator)
    if transform["oriented_repeat_count"] != repeats:
        raise ValueError("Slab generator did not retain the requested integer repeats")
    return result, transform, height


def _matches_seed(supplied, regenerated):
    """Compare ordering and inherited geometry, allowing a normal vacuum translation."""
    if supplied.species != regenerated.species or not np.allclose(
            supplied.lattice.matrix[:2], regenerated.lattice.matrix[:2], atol=1e-8, rtol=0):
        return False
    slab_geometry(supplied)
    left, right = supplied.cart_coords.copy(), regenerated.cart_coords.copy()
    left[:, 2] -= left[:, 2].min()
    right[:, 2] -= right[:, 2].min()
    # In-plane wrapping of the same atoms is a representation change.
    difference = (left - right) @ np.linalg.inv(regenerated.lattice.matrix)
    difference[:, :2] -= np.rint(difference[:, :2])
    return bool(np.allclose(difference @ regenerated.lattice.matrix, 0, atol=1e-7, rtol=0))


def extend_thickness(parent: Structure, seed: Structure, seed_metadata: dict, *,
                     n_layers=None, oriented_repeats=None, target_thickness_A=None,
                     vacuum_A=None):
    """Rebuild one complete-repeat thickness seed while retaining its termination.

    Supply one target: complete basal layers, integer OUC repeats, or an edge
    physical-thickness lower bound in A. Beta has one layer per OUC; IIb has
    two. Physical edge targets round upward to complete repeats, never trim.
    The supplied parent and seed must regenerate the recorded Stage 0 model.
    No relaxed displacements are copied or silently replaced.
    """
    targets = (n_layers, oriented_repeats, target_thickness_A)
    if sum(value is not None for value in targets) != 1:
        raise ValueError("Specify exactly one thickness target")
    phase = seed_metadata["phase"]
    if phase not in ("beta", "IIb"):
        raise ValueError("Only existing beta and IIb families are supported")
    basal = tuple(seed_metadata["hkl"]) == (0, 0, 1)
    source_parent = parent.copy()
    source_parent.add_site_property("parent_index", list(range(len(parent))))
    original_repeats = int(seed_metadata["transformation"]["oriented_repeat_count"])
    baseline, _, height = _regenerate(source_parent, seed_metadata, original_repeats)
    if not _matches_seed(seed, baseline):
        raise ValueError("Input is not the recorded bulk-derived seed; relaxed geometry cannot be extended implicitly")
    layers_per_repeat = 1 if phase == "beta" else 2
    if n_layers is not None:
        if (not basal or isinstance(n_layers, bool) or int(n_layers) != n_layers
                or n_layers <= 0 or n_layers % layers_per_repeat):
            raise ValueError("Basal layer target must contain complete parent stacking repeats")
        repeats = int(n_layers) // layers_per_repeat
    elif oriented_repeats is not None:
        if (isinstance(oriented_repeats, bool) or int(oriented_repeats) != oriented_repeats
                or oriented_repeats <= 0):
            raise ValueError("oriented_repeats must be a positive integer")
        repeats = int(oriented_repeats)
    else:
        if basal or not np.isfinite(target_thickness_A) or target_thickness_A <= 0:
            raise ValueError("Physical thickness targets are positive edge/prismatic lower bounds")
        unit, _, _ = _regenerate(source_parent, seed_metadata, 1)
        repeats = max(original_repeats, 1 + math.ceil(
            (target_thickness_A - np.ptp(unit.cart_coords[:, 2])) / height - 1e-10))
    if repeats < original_repeats:
        raise ValueError("Thickness extension cannot remove existing complete repeats")
    result, transform, _ = _regenerate(source_parent, seed_metadata, repeats)
    criteria = coordination_reference(source_parent, phase)
    diagnostics = diagnose_slab(result, source_parent, criteria, basal)
    if diagnostics["validation_errors"]:
        raise ValueError("Rebuilt slab is invalid: " + "; ".join(diagnostics["validation_errors"]))
    for side in ("lower", "upper"):
        fingerprint = json.loads(json.dumps(diagnostics["faces"][side]["surface_species_fingerprint"]))
        if fingerprint != \
                seed_metadata["diagnostics"]["faces"][side]["surface_species_fingerprint"]:
            raise ValueError(f"{side} termination changed during extension")
    if basal and diagnostics["complete_layer_count"] != repeats * layers_per_repeat:
        raise ValueError("Rebuilt slab does not contain the expected complete layers")
    context = dict(surface_id=seed_metadata["slab_id"], face_ids={
        side: seed_metadata[f"{side}_face_id"] for side in ("lower", "upper")})
    selected_vacuum = slab_geometry(seed)["vacuum_A"] if vacuum_A is None else vacuum_A
    result, vacuum_metadata = change_vacuum(result, selected_vacuum, face_context=context)
    transform["final_lattice_A"] = result.lattice.matrix.tolist()
    transform["translation_z_A"] += vacuum_metadata["translation_cartesian_A"][2]
    metadata = dict(builder="complete_parent_oriented_repeats", phase=phase,
                    parent_surface_id=seed_metadata["slab_id"], face_context=context,
                    face_identity_scope="same termination family; new thickness model",
                    termination_id=seed_metadata["termination_id"],
                    raw_cut_shift=seed_metadata["raw_cut_shift"],
                    parent_sha256=seed_metadata["parent_sha256"],
                    requested_n_layers=n_layers, requested_oriented_repeats=oriented_repeats,
                    target_thickness_A=target_thickness_A, oriented_repeat_count=repeats,
                    n_layers=diagnostics["complete_layer_count"],
                    layers_per_oriented_repeat=layers_per_repeat if basal else None,
                    transformation=transform, **slab_geometry(result),
                    reconstruction_status="UNRELAXED_BULK_DERIVED_SEED",
                    electronic_bulk_like_character="UNVERIFIED")
    return result, metadata
