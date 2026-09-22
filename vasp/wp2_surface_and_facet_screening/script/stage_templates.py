"""Pure, explicitly gated scientific configurations for later WP2 calculations.

These functions return dictionaries/structures; they never write input files,
inspect a runtime, invoke the frozen launcher, or interact with CMW. Approvals
are supplied review records, not a gate inferred from numerical results.
"""

import hashlib
import math

import numpy as np
from pymatgen.io.vasp import Poscar

from generate_surface_candidates import ACCEPTED
from surface_energetics import validate_paw_identities
from surface_layers import group_surface_layers

# Selected beta(001) 4L recipe; other thicknesses and tighter tolerances still
# require their own reviewed iteration budget (the accepted 6L retry used 300).
BASE = dict(GGA="PE", IVDW=12, PREC="Accurate", LASPH=True, LREAL=False,
            ISYM=0, ISPIN=1, EDIFF=1e-6, NELM=120, NELMIN=2, ALGO="Normal", AMIN=0.01,
            ISMEAR=0, SIGMA=0.05, ISTART=0, ICHARG=2,
            LCHARG=True, LWAVE=False, LVHAR=True)
ENERGY_CONVENTION = "energy(sigma->0)"


def approved_context(approval):
    """Validate and preserve the explicitly reviewed surface and Stage 0 face IDs."""
    if not isinstance(approval, dict) or approval.get("approved") is not True:
        raise ValueError("An explicitly approved surface record is required")
    for key in ("surface_id", "review_reference"):
        if not isinstance(approval.get(key), str) or not approval[key].strip():
            raise ValueError(f"Approval requires a nonempty {key}")
    faces = approval.get("face_ids")
    if not isinstance(faces, dict) or set(faces) != {"upper", "lower"}:
        raise ValueError("Both existing Stage 0 face IDs must be supplied")
    if any(not isinstance(value, str) or not value.strip() for value in faces.values()):
        raise ValueError("Face IDs must be nonempty strings")
    if faces["upper"] == faces["lower"]:
        raise ValueError("The two reaction-facing sides retain distinct face IDs")
    return dict(approval, face_ids=dict(faces))


def _numerical_settings(numerical, *, slab, selected=False):
    """Validate explicit ENCUT and Gamma mesh, with a Stage 02 selection gate."""
    encut = numerical.get("encut_eV")
    mesh = numerical.get("k_mesh")
    if isinstance(encut, bool) or not isinstance(encut, (int, float)) or not math.isfinite(encut) or encut <= 0:
        raise ValueError("An explicit positive finite ENCUT in eV is required")
    if not isinstance(mesh, (list, tuple)) or len(mesh) != 3 or any(type(n) is not int or n < 1 for n in mesh):
        raise ValueError("k_mesh must contain three positive integers")
    if slab and mesh[2] != 1:
        raise ValueError("Slab Gamma meshes must have N3=1")
    if selected and (not isinstance(numerical.get("selected_by_review"), str)
                     or not numerical["selected_by_review"].strip()):
        raise ValueError("Later stages require the explicit Stage 02 numerical-selection review")
    return float(encut), list(mesh)


def _magmom(values, atom_count):
    """Validate explicit collinear initial moments in current POSCAR site order."""
    if type(atom_count) is not int or atom_count < 1:
        raise ValueError("An explicit atom count is required for MAGMOM")
    array = np.asarray(values, dtype=float)
    if array.shape != (atom_count,) or not np.all(np.isfinite(array)):
        raise ValueError("MAGMOM must have one finite collinear initialization per atom")
    return array.tolist()


def surface_template(approval, numerical, *, mode="static", dipole_correction=True,
                     magmom=None, atom_count=None, smearing=None, retain_wavecar=False,
                     paw_identities=None):
    """Return a self-consistent PBE+D3(BJ) slab configuration with supplied sampling.

    The input slab must be centered along the third lattice vector, normal to
    the surface, as produced by the WP2 builders. DIPOL is fractional (0.5,0.5,0.5).
    Stage 02 may use this common primitive before settings are converged; later
    stages use the stricter wrappers below. No universal magnetic pattern is made.
    """
    context = approved_context(approval)
    encut, mesh = _numerical_settings(numerical, slab=True)
    if mode not in ("static", "relaxation"):
        raise ValueError("mode must be static or relaxation")
    if type(dipole_correction) is not bool or type(retain_wavecar) is not bool:
        raise ValueError("Dipole and WAVECAR choices must be explicit booleans")
    incar = dict(BASE, ENCUT=encut, ISIF=2, LDIPOL=dipole_correction,
                 IBRION=2 if mode == "relaxation" else -1,
                 NSW=100 if mode == "relaxation" else 0, LWAVE=retain_wavecar)
    if dipole_correction:
        incar.update(IDIPOL=3, DIPOL=[0.5, 0.5, 0.5])
    if mode == "relaxation":
        incar["EDIFFG"] = -0.02
    else:
        incar["LORBIT"] = 11
    if magmom is not None:
        incar.update(ISPIN=2, MAGMOM=_magmom(magmom, atom_count))
    smearing_note = "Working starting convention; electronic character remains to be established."
    if smearing is not None:
        if not isinstance(smearing.get("reason"), str) or not smearing["reason"].strip():
            raise ValueError("A smearing sensitivity change requires an explicit scientific reason")
        method, sigma = smearing.get("ismear"), smearing.get("sigma_eV")
        if type(method) is not int or method not in (-5, -1, 0, 1, 2):
            raise ValueError("Unsupported explicit smearing method")
        if isinstance(sigma, bool) or not isinstance(sigma, (int, float)) or not math.isfinite(sigma) or sigma <= 0:
            raise ValueError("SIGMA must be a finite positive value in eV")
        incar.update(ISMEAR=method, SIGMA=float(sigma))
        smearing_note = smearing["reason"]
    identities = None if paw_identities is None else validate_paw_identities(paw_identities)
    return dict(context=context, incar=incar, k_mesh=mesh, k_centering="Gamma",
                encut_eV=encut, xc="PBE", dispersion="D3(BJ)", vasp_version="6.6.1",
                formula_unit="ZnIn2S4", paw_identities=identities,
                energy_convention=ENERGY_CONVENTION, smearing_note=smearing_note,
                numerical_selection_review=numerical.get("selected_by_review"),
                integration_compatibility_id=numerical.get("integration_compatibility_id"),
                integration_review_reference=numerical.get("integration_review_reference"),
                requires_human_review=True)


def relaxation_template(approval, numerical, *, tier="candidate_screening",
                        validation_check=None, reason=None, **options):
    """Select screening or final force targets, with explicit shortlist checks."""
    context = approved_context(approval)
    _numerical_settings(numerical, slab=True, selected=True)
    if tier not in ("candidate_screening", "shortlist_validation", "final_relaxation"):
        raise ValueError("Unknown relaxation tier")
    if tier != "candidate_screening" and context.get("shortlisted") is not True:
        raise ValueError("Validation and final relaxation require explicit shortlisting")
    if tier == "shortlist_validation":
        if validation_check not in ("reconstruction", "spin", "relaxed_region"):
            raise ValueError("Select exactly one shortlist validation question")
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError("A shortlist validation check needs an explicit scientific reason")
    elif validation_check is not None or reason is not None:
        raise ValueError("Optional validation questions belong to shortlist_validation")
    result = surface_template(context, numerical, mode="relaxation", **options)
    if tier == "final_relaxation":
        result["incar"].update(EDIFF=1e-7, EDIFFG=-0.01)
    result.update(relaxation_tier=tier, validation_check=validation_check, validation_reason=reason)
    return result


def surface_static_template(approval, numerical, *, paw_identities,
                            band_edge_charge_reason=None, **options):
    """Prepare the final self-consistent static; retain WAVECAR only by justification."""
    context = approved_context(approval)
    _numerical_settings(numerical, slab=True, selected=True)
    if not isinstance(context.get("final_relaxation_review"), str) or not context["final_relaxation_review"].strip():
        raise ValueError("Final surface static requires the accepted final-relaxation review")
    if context.get("parent_phase") not in ACCEPTED:
        raise ValueError("An explicit parent phase is required for matched surface statics")
    for key in ("integration_compatibility_id", "integration_review_reference"):
        if not isinstance(numerical.get(key), str) or not numerical[key].strip():
            raise ValueError("A final surface static needs a reviewed bulk/slab integration convention")
    retain = band_edge_charge_reason is not None
    if retain and (context.get("shortlisted") is not True
                   or not isinstance(band_edge_charge_reason, str) or not band_edge_charge_reason.strip()):
        raise ValueError("Band-edge WAVECAR retention requires a shortlisted model and explicit reason")
    result = surface_template(context, numerical, mode="static", retain_wavecar=retain,
                              paw_identities=paw_identities, **options)
    result["incar"]["EDIFF"] = 1e-7
    result["band_edge_charge_reason"] = band_edge_charge_reason
    return result


def bulk_static_template(parent_phase, parent_poscar_bytes, parent_metadata, *, numerical,
                         paw_identities, matching_surface, integration_review_reference,
                         purpose="final_reference", convergence_case_review=None):
    """Prepare a VASP 6.6.1 static on the exact accepted WP1 parent geometry.

    Accepted Stage 0 parent bytes/metadata are explicit inputs. Historical bulk
    energies are never consumed. An explicit convergence_reference purpose permits
    one matched reference before numerical settings have been selected. The
    common integration review ID is required
    because a bulk three-dimensional mesh cannot equal a slab mesh in general.
    """
    if parent_phase not in ACCEPTED:
        raise ValueError("An explicit accepted parent phase is required")
    if not isinstance(parent_poscar_bytes, bytes):
        raise ValueError("Supply exact accepted parent POSCAR bytes")
    expected = ACCEPTED[parent_phase]["sha256"]
    if (hashlib.sha256(parent_poscar_bytes).hexdigest() != expected
            or parent_metadata.get("phase") != parent_phase
            or parent_metadata.get("source_sha256") != expected
            or parent_metadata.get("source_role") != "accepted relaxed WP1 parent"):
        raise ValueError("Parent geometry is not the accepted, unchanged WP1 parent")
    if purpose not in ("final_reference", "convergence_reference"):
        raise ValueError("Choose final_reference or one convergence_reference")
    if purpose == "convergence_reference":
        if not isinstance(convergence_case_review, str) or not convergence_case_review.strip():
            raise ValueError("A single convergence reference needs an explicit case review")
    elif convergence_case_review is not None:
        raise ValueError("A convergence-case review belongs to convergence_reference")
    context = approved_context(matching_surface.get("context"))
    if purpose == "final_reference" and not context.get("final_relaxation_review"):
        raise ValueError("A final bulk reference must identify the accepted final surface model")
    encut, mesh = _numerical_settings(numerical, slab=False, selected=purpose == "final_reference")
    identities = validate_paw_identities(paw_identities)
    if not isinstance(integration_review_reference, str) or not integration_review_reference.strip():
        raise ValueError("Explicit bulk/slab reciprocal-space compatibility review is required")
    compatibility_id = numerical.get("integration_compatibility_id")
    if not isinstance(compatibility_id, str) or not compatibility_id.strip():
        raise ValueError("An explicit bulk/slab integration compatibility identity is required")
    if (matching_surface.get("context", {}).get("parent_phase") != parent_phase
            or matching_surface.get("xc") != "PBE"
            or matching_surface.get("dispersion") != "D3(BJ)"
            or matching_surface.get("vasp_version") != "6.6.1"
            or matching_surface.get("encut_eV") != encut
            or matching_surface.get("paw_identities") != identities
            or matching_surface.get("energy_convention") != ENERGY_CONVENTION
            or matching_surface.get("integration_compatibility_id") != compatibility_id):
        raise ValueError("Bulk settings do not match the identified parent surface static")
    incar = dict(BASE, ENCUT=encut, IBRION=-1, NSW=0, LDIPOL=False)
    surface_incar = matching_surface["incar"]
    incar.update(ISMEAR=surface_incar["ISMEAR"], SIGMA=surface_incar["SIGMA"],
                 EDIFF=surface_incar["EDIFF"])
    return dict(structure=Poscar.from_str(parent_poscar_bytes.decode()).structure,
                parent_phase=parent_phase, parent_geometry_sha256=expected,
                incar=incar, k_mesh=mesh, k_centering="Gamma", encut_eV=encut,
                xc="PBE", dispersion="D3(BJ)", vasp_version="6.6.1",
                paw_identities=identities, formula_unit="ZnIn2S4",
                energy_convention=ENERGY_CONVENTION,
                integration_compatibility_id=compatibility_id,
                integration_review_reference=integration_review_reference,
                matched_surface_id=matching_surface["context"]["surface_id"],
                numerical_selection_review=numerical.get("selected_by_review"),
                purpose=purpose, convergence_case_review=convergence_case_review,
                requires_human_review=True)


def reconstruction_probe(structure, approval, *, supercell, surface_indices,
                         amplitude_A, reason, seed=0):
    """Replicate a shortlisted slab and perturb only explicitly selected surface sites.

    Indices refer to the supplied cell and are corroborated by the reusable
    geometric surface grouping. Displacements are in-plane, so normal center,
    slab thickness and actual atom-free vacuum remain unchanged. Composition
    increases by the requested lateral multiplicity; no atom is added or deleted.
    """
    context = approved_context(approval)
    if context.get("shortlisted") is not True or not isinstance(reason, str) or not reason.strip():
        raise ValueError("A reconstruction probe requires a shortlist and explicit scientific reason")
    if tuple(supercell) not in ((2, 1), (1, 2), (2, 2)):
        raise ValueError("Only explicitly requested 2x1, 1x2 or 2x2 probes are supported")
    if type(seed) is not int:
        raise ValueError("The deterministic seed must be an integer")
    if isinstance(amplitude_A, bool) or not isinstance(amplitude_A, (int, float)) or not 0 < amplitude_A <= 0.1:
        raise ValueError("Choose an explicit small displacement amplitude in (0, 0.1] Angstrom")
    selected = list(surface_indices)
    if not selected or any(type(i) is not int or not 0 <= i < len(structure) for i in selected) or len(set(selected)) != len(selected):
        raise ValueError("Supply distinct valid site indices from the input slab")
    grouping = group_surface_layers(structure)
    if grouping["status"] == "AMBIGUOUS":
        raise ValueError("Ambiguous geometric surface groups require review before perturbation")
    groups = grouping["groups"]
    allowed = set(groups["upper_surface"]) | set(groups["lower_surface"])
    if not set(selected) <= allowed:
        raise ValueError("Only geometrically identified surface atoms may be perturbed")
    lattice = structure.lattice.matrix
    if not np.allclose(lattice[:2, 2], 0, atol=1e-10) or not np.allclose(lattice[2, :2], 0, atol=1e-10):
        raise ValueError("Reconstruction probes require the centered Stage 0 surface presentation")
    z = structure.cart_coords[:, 2]
    if not np.isclose((z.min() + z.max()) / 2, structure.lattice.c / 2, atol=1e-8):
        raise ValueError("Input slab must be centered")
    output = structure.copy()
    output.add_site_property("probe_source_index", list(range(len(structure))))
    factors = [int(supercell[0]), int(supercell[1]), 1]
    output.make_supercell(factors)
    mapping = output.site_properties["probe_source_index"]
    perturbed = [i for i, source in enumerate(mapping) if source in selected]
    rng = np.random.default_rng(seed)
    angles = rng.uniform(0, 2 * np.pi, len(perturbed))
    displacements = np.column_stack((amplitude_A * np.cos(angles), amplitude_A * np.sin(angles), np.zeros(len(angles))))
    for index, displacement in zip(perturbed, displacements):
        output.translate_sites(index, displacement, frac_coords=False, to_unit_cell=False)
    metadata = dict(context=context, reconstruction_supercell=factors,
                    source_site_indices=list(mapping), selected_source_indices=selected,
                    perturbed_site_indices=perturbed, displacements_A=displacements.tolist(),
                    amplitude_A=float(amplitude_A), seed=seed, reason=reason,
                    composition_multiplier=factors[0] * factors[1],
                    interpretation="Starting probe only; unchanged replication does not exclude reconstruction.")
    return output, metadata


def spin_probe_templates(approval, numerical, *, atom_count, initializations, reason, **options):
    """Return only the explicitly supplied collinear MAGMOM probes in site order."""
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("A spin probe requires a scientific question")
    if not isinstance(initializations, dict) or not initializations:
        raise ValueError("Supply one or more named MAGMOM initializations")
    probes = []
    for name, values in initializations.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Every initialization needs a reproducible name")
        moments = _magmom(values, atom_count)
        result = relaxation_template(approval, numerical, tier="shortlist_validation",
                                     validation_check="spin", reason=reason,
                                     magmom=moments, atom_count=atom_count, **options)
        result["spin_initialization_id"] = name
        probes.append(result)
    return probes
