"""Canonical later-result columns, not scientific results or gate decisions.

None serializes as an empty CSV cell. Nested values use JSON; units occur in field
names. Unknown fields fail rather than silently diverging from these contracts.
"""

import csv
import io
import json


COLUMNS = {
    'convergence': '''case_id surface_id face_context dimension baseline_case_id target_value actual_value
        encut_eV k_mesh k_spacing_max_Ainv vacuum_A slab_thickness_A n_layers relaxed_region_mode
        energy_eV bulk_reference_id Gamma_pair_meV_A2 work_function_upper_eV work_function_lower_eV
        vacuum_upper_eV vacuum_lower_eV delta_Gamma_pair_meV_A2 delta_vacuum_upper_eV delta_vacuum_lower_eV
        delta_vacuum_reference_upper_eV delta_vacuum_reference_lower_eV geometry_change_flag
        reconstruction_flag spin_state_change_flag localisation_change_flag termination_change_flag
        numerical_status human_review_status notes'''.split(),
    'surface_properties': '''surface_id face_id parent_phase h k l termination_id reaction_face_side
        atom_count n_formula_units surface_area_A2 slab_thickness_A vacuum_A faces_equivalent spin_state
        E_slab_eV bulk_reference_id e_bulk_eV_fu Gamma_pair_meV_A2 gamma_single_meV_A2
        vacuum_upper_eV vacuum_lower_eV work_function_upper_eV work_function_lower_eV
        bulk_like_VBM_upper_eV bulk_like_CBM_upper_eV bulk_like_VBM_lower_eV bulk_like_CBM_lower_eV
        surface_state_flag electronic_localisation_notes relaxation_status numerical_convergence_status notes'''.split(),
    'relaxation': '''surface_id input_structure_id relaxation_tier spin_mode constraint_mode
        reconstruction_supercell initial_max_force_eVA final_max_force_eVA ionic_steps scf_iterations_total
        walltime_s normal_termination force_converged structure_change_rms_A max_displacement_A
        connectivity_changed termination_changed reconstruction_flag scientific_status notes'''.split(),
    'electronic': '''surface_id face_id static_result_id spin_mode fermi_eV sampled_integer_occupation_gap_eV gap_classification
        vacuum_upper_eV vacuum_lower_eV work_function_upper_eV work_function_lower_eV bulk_like_interior_available
        VBM_character CBM_character surface_state_present VBM_localisation CBM_localisation internal_field_flag
        alignment_status notes'''.split(),
}

# Shared meanings are defined once across tables. Unavailable evidence is None,
# including qualitative flags: absence of a reported change is not False.
FIELD_DOCS = {
    'case_id': 'Unique explicitly prepared numerical case identity.',
    'surface_id': 'Stable Stage 0 slab model identity, retained across numerical variants.',
    'face_id': 'Existing Stage 0 face identity for this reaction-facing record.',
    'face_context': 'JSON map containing surface_id and upper/lower face_ids.',
    'dimension': 'One varied dimension, or a single combined_verification.',
    'baseline_case_id': 'Explicit comparison case identity approved by human review.',
    'target_value': 'Requested value in the named dimension; structured thickness requests use JSON.',
    'actual_value': 'Realized value; complete-repeat rounding is reported, not hidden.',
    'encut_eV': 'Selected plane-wave kinetic-energy cutoff in eV.',
    'k_mesh': 'Gamma-centered integer mesh [N1,N2,1] for a slab.',
    'k_spacing_max_Ainv': 'Maximum in-plane |b_i|/N_i with b including 2*pi, inverse Angstrom.',
    'vacuum_A': 'Actual atom-free periodic separation along the slab normal, Angstrom.',
    'slab_thickness_A': 'Outer-atom normal-coordinate span, Angstrom.',
    'n_layers': 'Verified complete chemical layers where defined; not atomic-plane count.',
    'relaxed_region_mode': 'Full relaxation or explicitly reviewed fixed-interior treatment.',
    'energy_eV': 'Completed case total energy in its recorded matched energy convention, eV.',
    'bulk_reference_id': 'Matched completed parent bulk static result identity.',
    'Gamma_pair_meV_A2': 'Stoichiometric excess energy divided by ONE-face area, meV/Angstrom^2.',
    'gamma_single_meV_A2': 'Half Gamma_pair ONLY with established equivalent faces; otherwise missing.',
    'delta_Gamma_pair_meV_A2': 'Case minus baseline Gamma_pair in meV/Angstrom^2.',
    'geometry_change_flag': 'Explicit evidence of geometry change; alone does not classify reconstruction.',
    'reconstruction_flag': 'Explicit evidence of reconstruction, unknown if unassessed.',
    'spin_state_change_flag': 'Explicit evidence of a changed spin state, unknown if unassessed.',
    'localisation_change_flag': 'Explicit evidence of changed state localisation, unknown if unassessed.',
    'termination_change_flag': 'Explicit evidence of changed face termination, unknown if unassessed.',
    'numerical_status': 'Metric status only; never a scientific approval or job release.',
    'human_review_status': 'Human gate remains required regardless of scalar convergence.',
    'notes': 'Limitations, provenance, energy/alignment conventions and review rationale.',
    'parent_phase': 'Matched parent polymorph/polytype; not the lowest-energy unrelated phase.',
    'h': 'First Miller index in the accepted parent basis.',
    'k': 'Second Miller index in the accepted parent basis.',
    'l': 'Third Miller index in the accepted parent basis.',
    'termination_id': 'Stage 0 cut/termination family identity.',
    'reaction_face_side': 'Upper or lower side associated with face_id.',
    'atom_count': 'Number of atoms in the complete slab cell.',
    'n_formula_units': 'Verified integer ZnIn2S4 units in a stoichiometric slab.',
    'surface_area_A2': 'Magnitude of a cross b for ONE exposed face, Angstrom^2.',
    'faces_equivalent': 'Established face equivalence (True), established inequivalence (False), or unknown.',
    'spin_state': 'Evidence-based spin description; not guessed from initial MAGMOM.',
    'E_slab_eV': 'Final slab total energy in the declared energy convention, eV.',
    'e_bulk_eV_fu': 'Matched parent energy normalized to one ZnIn2S4 formula unit, eV.',
    'surface_state_flag': 'Projection-supported surface-state evidence; unknown if unassessed.',
    'electronic_localisation_notes': 'Projection evidence and limits of bulk-like/surface classification.',
    'relaxation_status': 'Reported relaxation completion evidence and limitations.',
    'numerical_convergence_status': 'Reported numerical convergence status without automatic approval.',
    'input_structure_id': 'Exact input model/geometry provenance identity.',
    'relaxation_tier': 'Candidate screening, question-driven shortlist validation, or final relaxation.',
    'spin_mode': 'Declared collinear/nonmagnetic spin convention.',
    'constraint_mode': 'Selective-dynamics convention and approved interior evidence.',
    'reconstruction_supercell': 'Explicit in-plane replication [N1,N2]; None if no probe.',
    'initial_max_force_eVA': 'Maximum initial atomic force norm in eV/Angstrom.',
    'final_max_force_eVA': 'Maximum final force norm over the relevant free degrees of freedom, eV/Angstrom.',
    'ionic_steps': 'Observed ionic step count.',
    'scf_iterations_total': 'Observed total electronic iterations across ionic steps.',
    'walltime_s': 'Observed elapsed wall time in seconds.',
    'normal_termination': 'Explicit normal termination evidence, not mere output existence.',
    'force_converged': 'Evidence against the applicable force criterion; unknown if unavailable.',
    'structure_change_rms_A': 'Mapped initial/final Cartesian displacement RMS, Angstrom.',
    'max_displacement_A': 'Largest mapped atom displacement, Angstrom.',
    'connectivity_changed': 'Comparison using a stated neighbor criterion; unknown if unassessed.',
    'termination_changed': 'Evidence-based change of termination, not inferred from RMS displacement.',
    'scientific_status': 'Evidence status; scientific acceptance remains a human decision.',
    'static_result_id': 'Self-consistent final-static source identity.',
    'fermi_eV': 'Fermi energy in the same SCF potential gauge as the vacuum levels, eV.',
    'sampled_integer_occupation_gap_eV': 'Separation from the projection reader ONLY when all sampled occupations are integer; distinct from the pilot frontier diagnostic, not a production gap.',
    'gap_classification': 'Sampled finite/near-zero/metallic-like/unresolved occupation diagnosis.',
    'bulk_like_interior_available': 'Defensible interior supported by geometry and independent evidence.',
    'VBM_character': 'Projection-supported occupied-edge character; not assigned by layer labels alone.',
    'CBM_character': 'Projection-supported empty-edge character; not assigned by layer labels alone.',
    'surface_state_present': 'Explicit state-localisation evidence; missing means unassessed.',
    'VBM_localisation': 'Occupied-edge layer/face/site projection evidence and classification.',
    'CBM_localisation': 'Empty-edge layer/face/site projection evidence and classification.',
    'internal_field_flag': 'Persistent interior electric-field evidence; missing means unassessed.',
    'alignment_status': 'Usable face reference or explicit limitation; never a forced alignment.',
}
for side in ('upper', 'lower'):
    FIELD_DOCS[f'vacuum_{side}_eV'] = f'Reviewed {side}-face LVHAR plateau in the source SCF gauge, eV.'
    FIELD_DOCS[f'work_function_{side}_eV'] = f'{side}-face vacuum minus same-source Fermi energy, eV.'
    FIELD_DOCS[f'delta_vacuum_{side}_eV'] = f'Case minus baseline {side} vacuum; only meaningful with a common potential gauge.'
    FIELD_DOCS[f'delta_vacuum_reference_{side}_eV'] = f'Case minus baseline gauge-invariant {side}-face quantity named in notes, eV.'
    for edge in ('VBM', 'CBM'):
        FIELD_DOCS[f'bulk_like_{edge}_{side}_eV'] = f'Candidate bulk-like {edge} energy minus {side}-face vacuum, eV; missing on any alignment limitation.'


def schema_descriptor(kind):
    """Return an ordered machine-readable schema, with every column documented."""
    return {'schema': kind, 'version': 2 if kind == 'electronic' else 1, 'missing_value': None,
            'fields': [{'name': key, 'description': FIELD_DOCS[key]} for key in COLUMNS[kind]]}


def schema_record(kind, values):
    """Normalize a partial evidenced record without inventing absent values."""
    if kind == 'electronic' and 'sampled_gap_eV' in values:
        raise ValueError('Legacy sampled_gap_eV is retired; recompute from identified source arrays')
    unknown = set(values) - set(COLUMNS[kind])
    if unknown:
        raise ValueError(f'Unknown {kind} columns: {sorted(unknown)}')
    return {key: values.get(key) for key in COLUMNS[kind]}


def csv_text(kind, rows):
    """Serialize supplied rows only; an empty iterable produces a schema header."""
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=COLUMNS[kind], lineterminator='\n')
    writer.writeheader()
    for row in rows:
        record = schema_record(kind, row)
        writer.writerow({key: json.dumps(value, sort_keys=True, allow_nan=False)
                         if isinstance(value, (dict, list, tuple)) else value
                         for key, value in record.items()})
    return stream.getvalue()
