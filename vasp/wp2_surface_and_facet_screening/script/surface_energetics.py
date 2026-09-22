"""Paired surface excess energies from explicit, matched ZnIn2S4 static records.

Inputs contain observed energies and calculation identities, never launcher or
queue state. PAW identity values must identify the actual datasets (for example,
their verified SHA256 values), not merely element names or a shared directory.
"""

import hashlib
import math
import re
from pathlib import Path

import numpy as np
from pymatgen.io.vasp import Kpoints, Poscar, Vasprun

FORMULA = {"Zn": 1, "In": 2, "S": 4}
MATCH_FIELDS = ("parent_phase", "xc", "dispersion", "paw_identities", "encut_eV",
                "energy_convention", "vasp_version", "formula_unit")
ENERGY_CONVENTIONS = ("energy(sigma->0)", "free_energy_TOTEN")


def formula_units(composition):
    """Return the integer ZnIn2S4 count, or None for a nonstoichiometric model."""
    if not isinstance(composition, dict) or set(composition) != set(FORMULA):
        return None
    counts = list(composition.values())
    if any(isinstance(n, bool) or not isinstance(n, (int, float))
           or not math.isfinite(n) or n <= 0 or int(n) != n for n in counts):
        return None
    units = composition["Zn"]
    if any(composition[element] != amount * units for element, amount in FORMULA.items()):
        return None
    return int(units)


def validate_paw_identities(identities):
    """Require explicit identities for each species; never read licensed data."""
    if not isinstance(identities, dict) or set(identities) != set(FORMULA):
        raise ValueError("PAW identities must explicitly cover Zn, In and S")
    if any(not isinstance(value, str) or not value.strip() for value in identities.values()):
        raise ValueError("Every PAW identity must be a nonempty verified dataset identifier")
    return dict(identities)


def validate_matched_reference(slab, bulk):
    """Reject unknown/mismatched method, parent, normalization or integration.

    Both records use total-cell energies. The common formula unit is ZnIn2S4.
    Distinct bulk/slab meshes are expected, so a reviewed common
    ``integration_compatibility_id`` records their compatibility rather than
    incorrectly requiring equal three-dimensional meshes.
    """
    for record in (slab, bulk):
        missing = [name for name in MATCH_FIELDS if record.get(name) in (None, "", {})]
        if missing:
            raise ValueError("Undefined matched-reference fields: " + ", ".join(missing))
        validate_paw_identities(record["paw_identities"])
        if record["parent_phase"] not in ("beta", "IIb", "spinel"):
            raise ValueError("Unknown parent phase")
        if record["vasp_version"] != "6.6.1":
            raise ValueError("Final WP2 references require VASP 6.6.1")
        if record["xc"] != "PBE" or record["dispersion"] != "D3(BJ)":
            raise ValueError("This WP2 convention requires PBE+D3(BJ)")
        if record["formula_unit"] != "ZnIn2S4":
            raise ValueError("Formula-unit normalization is undefined or incompatible")
        if record["energy_convention"] not in ENERGY_CONVENTIONS:
            raise ValueError("Energy convention must explicitly identify a supported VASP energy")
        if isinstance(record["encut_eV"], bool) or not isinstance(record["encut_eV"], (int, float)) or record["encut_eV"] <= 0:
            raise ValueError("ENCUT must be a positive value in eV")
        if not math.isfinite(record["encut_eV"]):
            raise ValueError("ENCUT must be finite")
        units = formula_units(record.get("composition"))
        if units is None:
            raise ValueError("Formula-unit normalization requires a stoichiometric composition")
        if record.get("n_formula_units", units) != units:
            raise ValueError("Declared formula-unit normalization disagrees with composition")
        if not isinstance(record.get("integration_compatibility_id"), str) or not record["integration_compatibility_id"].strip():
            raise ValueError("A reviewed bulk/slab integration compatibility identity is required")
        if not isinstance(record.get("integration_review_reference"), str) or not record["integration_review_reference"].strip():
            raise ValueError("Bulk/slab integration compatibility requires explicit review evidence")
    mismatches = [name for name in MATCH_FIELDS if slab[name] != bulk[name]]
    if slab["integration_compatibility_id"] != bulk["integration_compatibility_id"]:
        mismatches.append("integration_compatibility_id")
    if mismatches:
        raise ValueError("Unmatched bulk reference: " + ", ".join(mismatches))


def surface_energy(slab, bulk, *, faces_equivalent=None, face_equivalence_evidence=None):
    """Compute (E_slab - n E_bulk/n_bulk)/A, where A is ONE-face area in Å².

    ``slab`` supplies E_slab_eV, composition and surface_area_A2; ``bulk``
    supplies E_bulk_eV, composition and bulk_reference_id. Both also supply
    the identities checked by validate_matched_reference. Missing/unknown face
    equivalence leaves single-face energies None. Positive equivalence requires
    explicit evidence; two face IDs remain intact even when equivalent.
    """
    surface_id, faces = slab.get("surface_id"), slab.get("face_ids")
    if (not isinstance(surface_id, str) or not surface_id.strip()
            or not isinstance(faces, dict) or set(faces) != {"upper", "lower"}
            or any(not isinstance(value, str) or not value.strip() for value in faces.values())
            or faces["upper"] == faces["lower"]):
        raise ValueError("Supply the surface identity and both distinct Stage 0 face identities")
    units = formula_units(slab.get("composition"))
    empty = dict(surface_id=surface_id, face_context=dict(surface_id=surface_id, face_ids=dict(faces)),
                 parent_phase=slab.get("parent_phase"), n_formula_units=units,
                 surface_area_A2=slab.get("surface_area_A2"), E_slab_eV=slab.get("E_slab_eV"),
                 bulk_reference_id=bulk.get("bulk_reference_id"), e_bulk_eV_fu=None,
                 Gamma_pair_eV_A2=None, Gamma_pair_meV_A2=None,
                 gamma_single_eV_A2=None, gamma_single_meV_A2=None,
                 faces_equivalent=faces_equivalent, matched_reference_validated=False)
    if units is None:
        return dict(empty, status="UNSUPPORTED_NONSTOICHIOMETRIC",
                    notes="No chemical-potential correction or surface energy was assigned.")
    validate_matched_reference(slab, bulk)
    if not bulk.get("bulk_reference_id"):
        raise ValueError("An identified matched bulk static is required")
    area = slab.get("surface_area_A2")
    energy = slab.get("E_slab_eV")
    bulk_energy = bulk.get("E_bulk_eV")
    if any(isinstance(value, bool) or not isinstance(value, (int, float))
           or not math.isfinite(value) for value in (area, energy, bulk_energy)) or area <= 0:
        raise ValueError("Finite total energies and a positive one-face area are required")
    if faces_equivalent is not None and type(faces_equivalent) is not bool:
        raise ValueError("Face equivalence must be True, False or None")
    if faces_equivalent is True and (not isinstance(face_equivalence_evidence, str)
                                     or not face_equivalence_evidence.strip()):
        raise ValueError("Equivalent faces require explicit symmetry/structural evidence")
    bulk_per_unit = bulk_energy / formula_units(bulk["composition"])
    paired = (energy - units * bulk_per_unit) / area
    single = paired / 2 if faces_equivalent is True else None
    return dict(empty, status="COMPUTED_REQUIRES_HUMAN_REVIEW", matched_reference_validated=True,
                e_bulk_eV_fu=bulk_per_unit,
                Gamma_pair_eV_A2=paired, Gamma_pair_meV_A2=paired * 1000,
                gamma_single_eV_A2=single,
                gamma_single_meV_A2=None if single is None else single * 1000,
                notes="Gamma_pair is the combined two-face excess per one-face area; no surface is promoted.")


def surface_property_rows(result, *, properties=None):
    """Map a model's paired energetics to two canonical reaction-face records.

    Gamma_pair remains the same model-level quantity in both face rows; it must
    not be summed across them. Optional measured properties use the canonical
    surface_properties schema; unavailable fields remain None.
    """
    from result_schemas import schema_record

    shared = dict(properties or {})
    keys = ("surface_id", "parent_phase", "n_formula_units", "surface_area_A2", "E_slab_eV",
            "bulk_reference_id", "e_bulk_eV_fu", "Gamma_pair_meV_A2", "gamma_single_meV_A2",
            "faces_equivalent")
    shared.update({key: result[key] for key in keys})
    notes = "; ".join(value for value in (shared.get("notes"), result["status"], result["notes"])
                      if value)
    shared["notes"] = notes + " Same paired excess repeated for face context; do not sum face rows."
    return [schema_record("surface_properties", dict(shared, face_id=face_id, reaction_face_side=side))
            for side, face_id in result["face_context"]["face_ids"].items()]


def _paw_identities_from_file(path):
    """Hash complete local PAW blocks; return identifiers only, never dataset bytes."""
    data = path.read_bytes()
    blocks = re.findall(rb".*?End of Dataset[^\n]*\n", data, flags=re.S)
    if not blocks or b"".join(blocks).strip() != data.strip():
        raise ValueError("Cannot identify complete local PAW dataset blocks")
    identities = {}
    for block in blocks:
        title = re.search(rb"TITEL\s*=\s*PAW_PBE\s+(\S+)\s", block)
        if title is None:
            raise ValueError("Expected a verified PAW_PBE dataset")
        variant = title.group(1).decode()
        element = variant.split("_")[0]
        identity = "sha256:" + hashlib.sha256(block).hexdigest()
        if element in identities and identities[element] != identity:
            raise ValueError("Repeated species blocks contain inconsistent PAW datasets")
        identities[element] = identity
    return validate_paw_identities(identities)


def read_static_record(directory, template, *, kind, result_id, expected_structure=None):
    """Read one explicitly identified completed static; reject incomplete evidence.

    The caller supplies a template and the accepted input structure (bulk
    templates already contain it). Exact local POTCAR block hashes must equal
    the template's ``sha256:...`` identities; POTCAR bytes never leave this
    function. No directory is discovered automatically and no file is written.
    This reader accepts fixed-geometry statics only, never relaxations.
    """
    if kind not in ("surface", "bulk") or not isinstance(result_id, str) or not result_id.strip():
        raise ValueError("Explicit static kind and result identity are required")
    if template["incar"].get("NSW") != 0 or template["incar"].get("ICHARG") == 11:
        raise ValueError("Only the self-consistent fixed-geometry static protocol is supported")
    expected_structure = expected_structure if expected_structure is not None else template.get("structure")
    if expected_structure is None:
        raise ValueError("The accepted static input structure must be supplied explicitly")
    directory = Path(directory)
    with (directory / "OUTCAR").open("rb") as stream:
        stream.seek(0, 2)
        stream.seek(max(0, stream.tell() - 2_000_000))
        outcar = stream.read().decode(errors="replace")
    normal = "General timing and accounting informations for this job:" in outcar
    ediff = "aborting loop because EDIFF is reached" in outcar
    if not normal or not ediff:
        raise ValueError("Normal termination and final electronic-convergence evidence are required")
    run = Vasprun(directory / "vasprun.xml", parse_dos=False, parse_eigen=False,
                  parse_projected_eigen=False, parse_potcar_file=False, exception_on_bad_xml=True)
    if run.vasp_version != "6.6.1" or not run.converged_electronic or len(run.ionic_steps) != 1:
        raise ValueError("A completed VASP 6.6.1 electronic static is required")
    iterations = len(run.ionic_steps[-1]["electronic_steps"])
    if not 2 <= iterations < int(template["incar"]["NELM"]):
        raise ValueError("Static electronic iteration evidence is insufficient or exhausted")
    observed = dict(run.incar)
    observed.update(run.parameters)
    for key in ("NELECT", "NUPDOWN"):
        if key in run.incar and key not in template["incar"]:
            raise ValueError(f"Unreviewed charge or spin-population override: {key}")
    if run.incar.get("EFIELD", 0) != 0 and "EFIELD" not in template["incar"]:
        raise ValueError("An external electric field is outside the reviewed static protocol")
    for key, expected in template["incar"].items():
        actual = observed.get(key)
        if actual is None:
            raise ValueError(f"Missing observed static setting: {key}")
        equal = (str(actual).upper() == expected.upper() if isinstance(expected, str)
                 else np.allclose(actual, expected, rtol=1e-10, atol=1e-12))
        if key == "PREC":
            equal = str(actual).upper()[:6] == expected.upper()[:6]
        if not equal:
            raise ValueError(f"Observed static setting differs from the template: {key}")
    if (observed.get("LHFCALC", False) or observed.get("LDAU", False)
            or observed.get("LSORBIT", False) or observed.get("LNONCOLLINEAR", False)
            or str(observed.get("METAGGA", "none")).lower() not in ("none", "--", "")):
        raise ValueError("Unexpected hybrid, meta-GGA, +U or relativistic selector")
    for mesh in (Kpoints.from_file(directory / "KPOINTS"), run.kpoints):
        if (mesh.style.name != "Gamma" or len(mesh.kpts) != 1
                or list(mesh.kpts[0]) != list(template["k_mesh"])
                or not np.allclose(mesh.kpts_shift, 0, atol=0, rtol=0)):
            raise ValueError("Observed Gamma mesh differs from the explicit template")
    identities = _paw_identities_from_file(directory / "POTCAR")
    if identities != template.get("paw_identities"):
        raise ValueError("Actual PAW identities do not match the explicit template")
    poscar_bytes = (directory / "POSCAR").read_bytes()
    initial = Poscar.from_str(poscar_bytes.decode()).structure
    for structure in (initial, run.initial_structure, run.final_structure):
        if (structure.species != expected_structure.species
                or not np.allclose(structure.lattice.matrix, expected_structure.lattice.matrix, atol=1e-8, rtol=0)):
            raise ValueError("Static input/final cell or species differ from the accepted geometry")
        delta = structure.frac_coords - expected_structure.frac_coords
        delta -= np.rint(delta)
        if not np.allclose(delta, 0, atol=1e-8, rtol=0):
            raise ValueError("The static geometry differs from the accepted input")
    convention = template.get("energy_convention")
    if convention not in ENERGY_CONVENTIONS:
        raise ValueError("A supported final energy convention must be explicit")
    # Vasprun.final_energy corrects the known ionic e_0_energy XML inconsistency.
    energy = float(run.final_energy if convention == "energy(sigma->0)"
                   else run.ionic_steps[-1]["e_fr_energy"])
    if not math.isfinite(energy):
        raise ValueError("The observed static energy is not finite")
    composition = {element: int(value) for element, value in run.final_structure.composition.get_el_amt_dict().items()}
    context = template.get("context", {})
    record = {key: template.get(key) for key in MATCH_FIELDS if key != "parent_phase"}
    record.update(parent_phase=context.get("parent_phase", template.get("parent_phase")),
                  composition=composition, n_formula_units=formula_units(composition),
                  integration_compatibility_id=template.get("integration_compatibility_id"),
                  integration_review_reference=template.get("integration_review_reference"),
                  case_id=result_id, static_result_id=result_id, energy_eV=energy,
                  normal_termination=True, electronic_converged=True,
                  input_settings_verified=True,
                  electronic_iterations=iterations, input_structure_sha256=hashlib.sha256(poscar_bytes).hexdigest(),
                  source_directory=str(directory.resolve()), status="COMPLETED_STATIC_REQUIRES_HUMAN_REVIEW")
    if kind == "surface":
        record.update(surface_id=context.get("surface_id"), face_ids=context.get("face_ids"),
                      E_slab_eV=energy, surface_area_A2=float(np.linalg.norm(np.cross(initial.lattice.matrix[0], initial.lattice.matrix[1]))))
    else:
        record.update(bulk_reference_id=result_id, E_bulk_eV=energy)
    return record
