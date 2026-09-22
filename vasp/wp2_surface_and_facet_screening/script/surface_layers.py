"""Geometry-based surface groups and explicitly requested relaxation masks.

Layer labels are spatial regions, not proof of bulk-like electronic states.
Thin slabs keep missing/ambiguous groups rather than acquiring an invented
interior. Thresholds are explicit proposal parameters in angstrom.
"""

import copy
import hashlib
import json

import numpy as np

from surface_builders import slab_geometry
from surface_model_utils import FACE_DEPTH, PLANE_TOL, plane_sequence

GROUP_LABELS = ("lower_surface", "lower_subsurface", "interior",
                "upper_subsurface", "upper_surface")


def group_surface_layers(slab, *, surface_depth_A=FACE_DEPTH, subsurface_depth_A=FACE_DEPTH,
                         min_interior_span_A=FACE_DEPTH, face_context=None):
    """Partition whole atomic planes by distance from each exposed face.

    Surface/subsurface depths are geometric proposal parameters. When bands
    overlap, a plane is assigned only to its nearer face; an exactly central
    ambiguous plane stays unassigned. A candidate interior needs at least
    three planes and the requested physical span beyond both face regions.
    A narrow middle is left unassigned rather than relabeled as bulk-like.
    geometry_id binds evidence to the ordered species, cell, and Cartesian
    coordinates exactly as supplied; site properties are excluded.
    """
    parameters = dict(surface_depth_A=surface_depth_A, subsurface_depth_A=subsurface_depth_A,
                      min_interior_span_A=min_interior_span_A, plane_tolerance_A=PLANE_TOL)
    if any(not np.isfinite(value) or value <= 0 for value in parameters.values()):
        raise ValueError("Grouping depths and tolerances must be finite and positive")
    geometry = slab_geometry(slab)
    identity = dict(species=[site.species_string for site in slab],
                    lattice_A=slab.lattice.matrix.tolist(), cartesian_positions_A=slab.cart_coords.tolist())
    geometry_id = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":"),
                                            allow_nan=False).encode()).hexdigest()
    heights = np.asarray(geometry["heights_A"])
    planes = plane_sequence(slab, heights)
    thickness = geometry["slab_thickness_A"]
    groups = {label: [] for label in GROUP_LABELS}
    labels = [None] * len(slab)
    interior_planes = []
    ambiguous = []
    outer_depth = surface_depth_A + subsurface_depth_A
    for plane in planes:
        lower, upper = plane["z_A"], thickness - plane["z_A"]
        nearest = min(lower, upper)
        if nearest > outer_depth:
            interior_planes.append(plane)
            continue
        if abs(lower - upper) <= PLANE_TOL:
            ambiguous.extend(plane["site_indices"])
            continue
        side = "lower" if lower < upper else "upper"
        region = "surface" if nearest <= surface_depth_A else "subsurface"
        label = f"{side}_{region}"
        groups[label].extend(plane["site_indices"])
        for index in plane["site_indices"]:
            labels[index] = label
    span = (interior_planes[-1]["z_A"] - interior_planes[0]["z_A"]
            if interior_planes else 0.0)
    interior_available = len(interior_planes) >= 3 and span >= min_interior_span_A
    if interior_available:
        for plane in interior_planes:
            groups["interior"].extend(plane["site_indices"])
            for index in plane["site_indices"]:
                labels[index] = "interior"
    unassigned = [index for index, label in enumerate(labels) if label is None]
    status = ("AMBIGUOUS" if ambiguous else "GEOMETRY_GROUPED" if interior_available
              else "NO_BULK_LIKE_INTERIOR")
    context = copy.deepcopy(face_context or {})
    return dict(groups=groups, site_labels=labels, unassigned_site_indices=unassigned,
                status=status, interior_available=interior_available,
                interior_interpretation="geometric candidate; bulk-like electronic character unverified",
                interior_span_A=span if interior_available else None,
                interior_plane_count=len(interior_planes), thickness_A=thickness,
                parameters=parameters, surface_id=context.get("surface_id"),
                face_ids=context.get("face_ids", {}), atom_count=len(slab), geometry_id=geometry_id)


def propose_relaxed_region(slab, *, mode, explicitly_requested=False,
                           interior_evidence=None, grouping=None, face_context=None):
    """Return a requested Selective Dynamics mask and its evidence record.

    True means free along that lattice-coordinate direction. Fixed-region
    proposals require geometry with material free on both faces and explicit
    corroborating bulk-like evidence for the same surface, geometry_id, and sites.
    The evidence must name an established criterion and reference; this routine
    does not establish bulk-like character from geometry alone.
    """
    if not explicitly_requested:
        raise ValueError("A relaxation-mask proposal must be explicitly requested")
    if mode not in ("full", "central_fixed"):
        raise ValueError("mode must be full or central_fixed")
    mask = np.ones((len(slab), 3), dtype=bool)
    if mode == "full":
        return mask.tolist(), dict(mode=mode, fixed_site_indices=[], face_context=copy.deepcopy(face_context))
    current = group_surface_layers(slab, face_context=face_context,
                                   **({key: value for key, value in grouping["parameters"].items()
                                       if key != "plane_tolerance_A"} if grouping else {}))
    if grouping is not None and (grouping["groups"] != current["groups"]
                                 or grouping["site_labels"] != current["site_labels"]
                                 or grouping.get("geometry_id") != current["geometry_id"]):
        raise ValueError("Supplied grouping no longer matches the current geometry")
    if (current["status"] != "GEOMETRY_GROUPED" or not current["interior_available"]
            or current["unassigned_site_indices"]):
        raise ValueError("No unambiguous, sufficiently thick central region is available")
    if any(not current["groups"][label] for label in GROUP_LABELS if label != "interior"):
        raise ValueError("A fixed region would leave insufficient surface/subsurface material free")
    evidence = interior_evidence or {}
    if (evidence.get("bulk_like_established") is not True or not evidence.get("reference")
            or not evidence.get("criterion") or not current["surface_id"]
            or evidence.get("surface_id") != current["surface_id"]
            or evidence.get("geometry_id") != current["geometry_id"]
            or sorted(evidence.get("site_indices", [])) != sorted(current["groups"]["interior"])):
        raise ValueError("Matching, explicit bulk-like evidence for the proposed interior is required")
    fixed = current["groups"]["interior"]
    mask[fixed] = False
    return mask.tolist(), dict(mode=mode, fixed_site_indices=fixed, grouping=current,
                              interior_evidence=copy.deepcopy(evidence),
                              scientific_status="PROPOSED_REQUIRES_HUMAN_REVIEW")
