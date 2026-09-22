"""Compare later completed cases using project working thresholds, never release jobs.

Energy-pair convergence uses 1 meV/A^2; applicable vacuum-referenced electronic
quantities use 0.05 eV. Raw LOCPOT offsets from separate calculations are not
comparable without an explicit common gauge. Work functions are gauge invariant.
"""

import math

from result_schemas import COLUMNS, schema_record

QUALITATIVE_FLAGS = ('reconstruction_flag', 'spin_state_change_flag',
                     'localisation_change_flag', 'termination_change_flag')


def compare_cases(baseline, case, *, electronic_relevant, electronic_quantity='work_function',
                  gamma_tolerance_meV_A2=1.0, electronic_tolerance_eV=0.05):
    """Return a canonical convergence record with no automatic scientific PASS.

    Inputs retain raw completion/energetics evidence in addition to schema fields:
    normal_termination, electronic_converged, and matched_reference_validated must
    each be True. These are supplied by completed-output analysis, never inferred
    from file presence. Qualitative flags are True/False/None with evidence notes.
    electronic_quantity may be work_function, bulk_like_VBM, or bulk_like_CBM.
    """
    if not isinstance(electronic_relevant, bool):
        raise ValueError('Explicit electronic_relevant decision required')
    if electronic_quantity not in ('work_function', 'bulk_like_VBM', 'bulk_like_CBM'):
        raise ValueError('Choose a gauge-invariant vacuum-referenced quantity')
    if any(not math.isfinite(x) or x <= 0 for x in (gamma_tolerance_meV_A2, electronic_tolerance_eV)):
        raise ValueError('Thresholds must be finite and positive')
    row = schema_record('convergence', {k: v for k, v in case.items() if k in COLUMNS['convergence']})
    row.update(baseline_case_id=baseline.get('case_id'), human_review_status='REQUIRED')
    notes = [case.get('notes') or '', 'Project working thresholds; no automatic scientific gate.']

    def finish(status, reason):
        row.update(numerical_status=status, notes=' '.join(x for x in [*notes, reason] if x))
        return row

    for key in QUALITATIVE_FLAGS:
        if case.get(key) is not None and not isinstance(case[key], bool):
            return finish('AMBIGUOUS', f'{key} requires True/False/None evidence.')
    if any(case.get(key) is True for key in QUALITATIVE_FLAGS):
        return finish('HUMAN_REVIEW_REQUIRED', 'Qualitative state/termination change overrides scalar comparison.')
    for key in ('case_id', 'surface_id', 'face_context'):
        if not baseline.get(key) or not case.get(key):
            return finish('INSUFFICIENT_DATA', f'Missing explicit {key}.')
    if (baseline['surface_id'] != case['surface_id'] or baseline['face_context'] != case['face_context']
            or baseline['case_id'] == case['case_id']):
        return finish('AMBIGUOUS', 'Comparison model/face identities disagree, or case compares with itself.')
    if case.get('baseline_case_id') not in (None, '', baseline['case_id']):
        return finish('AMBIGUOUS', 'Supplied comparison baseline does not match the prepared baseline identity.')
    dimension = case.get('dimension')
    if dimension not in ('kpoints', 'encut', 'vacuum', 'slab_thickness', 'relaxed_region', 'combined_verification'):
        return finish('AMBIGUOUS', 'An explicit prepared comparison dimension is required.')
    protocol = baseline.get('comparison_protocol')
    if not protocol or not case.get('comparison_protocol'):
        return finish('INSUFFICIENT_DATA', 'Observed/input-verified comparison protocol is required.')
    if protocol != case['comparison_protocol']:
        return finish('AMBIGUOUS', 'Spin/smearing/dipole/SCF/ionic/method protocols differ.')
    invariant = {'encut': 'encut_eV', 'kpoints': 'k_mesh', 'vacuum': 'vacuum_A',
                 'slab_thickness': 'slab_thickness_A', 'relaxed_region': 'relaxed_region_mode'}
    if dimension != 'combined_verification':
        for varied, key in invariant.items():
            if varied == dimension:
                continue
            left, right = baseline.get(key), case.get(key)
            if left is None or right is None:
                return finish('INSUFFICIENT_DATA', f'Missing invariant baseline/case {key}.')
            equal = (math.isclose(left, right, abs_tol=1e-7, rel_tol=0)
                     if isinstance(left, (int, float)) and isinstance(right, (int, float)) else left == right)
            if not equal:
                return finish('AMBIGUOUS', f'Unreviewed second dimension changed: {key}.')
    elif not case.get('combined_selection_review'):
        return finish('INSUFFICIENT_DATA', 'Combined verification requires the prior selection review.')
    for source in (baseline, case):
        if any(source.get(key) is not True for key in
               ('normal_termination', 'electronic_converged', 'matched_reference_validated')):
            return finish('INSUFFICIENT_DATA', 'Completed SCF and matched-parent energetics evidence required for both cases.')
        if protocol.get('incar', {}).get('NSW', 0) > 0 and source.get('force_converged') is not True:
            return finish('INSUFFICIENT_DATA', 'Relaxation comparison requires force convergence, not NSW exhaustion.')
        if electronic_relevant and source.get('alignment_status') != 'FACE_ALIGNED':
            return finish('AMBIGUOUS', 'Affirmative reviewed alignment of both face electronic references is required.')
    pairs = [('Gamma_pair_meV_A2', 'delta_Gamma_pair_meV_A2', gamma_tolerance_meV_A2)]
    if electronic_relevant:
        pairs += [(f'{electronic_quantity}_{side}_eV', f'delta_vacuum_reference_{side}_eV', electronic_tolerance_eV)
                  for side in ('upper', 'lower')]
        notes.append(f'Electronic convergence quantity: {electronic_quantity}, case minus baseline, both faces.')
    failed = False
    for key, output, tolerance in pairs:
        first, second = baseline.get(key), case.get(key)
        if first is None or second is None:
            return finish('INSUFFICIENT_DATA', f'Missing applicable {key}; unavailable values are not zero.')
        if (isinstance(first, bool) or isinstance(second, bool)
                or not isinstance(first, (int, float)) or not isinstance(second, (int, float))
                or not math.isfinite(first) or not math.isfinite(second)):
            return finish('AMBIGUOUS', f'Invalid/nonfinite {key}.')
        row[output] = second - first
        failed |= abs(row[output]) > tolerance + 1e-12
    # Separate SCF potential zeros can drift. Expose raw changes only when a
    # caller supplies a common alignment reference rather than hiding that issue.
    gauge = baseline.get('potential_gauge_reference')
    if gauge and gauge == case.get('potential_gauge_reference'):
        for side in ('upper', 'lower'):
            key = f'vacuum_{side}_eV'
            values = (baseline.get(key), case.get(key))
            if all(isinstance(v, (int, float)) and math.isfinite(v) for v in values):
                row[f'delta_{key}'] = values[1] - values[0]
    if any(case.get(key) is None for key in QUALITATIVE_FLAGS):
        notes.append('Unassessed qualitative flags remain unknown; review is still required.')
    return finish('NOT_CONVERGED' if failed else 'NUMERIC_PASS',
                  'Applicable scalar thresholds exceeded.' if failed else 'Applicable scalar thresholds satisfied.')


def convergence_record(prepared, completed, energetics, *, alignment=None, qualitative=None):
    """Join identified observed outputs to one prepared model, without writing rows.

    The completed reader must verify the input settings and identify its case_id.
    Qualitative evidence is supplied explicitly; no reconstruction/localisation
    classification is fabricated from one scalar or an unconverged geometry.
    The returned evidence-rich dictionary feeds compare_cases; schema_record or
    csv_text selects canonical CSV fields after comparison.
    """
    if completed.get('case_id') != prepared['case_id'] or completed.get('input_settings_verified') is not True:
        raise ValueError('Observed settings must be verified against this prepared case')
    if (energetics.get('surface_id') != prepared['surface_id']
            or energetics.get('face_context') != prepared['face_context']):
        raise ValueError('Energetics belongs to a different surface/face context')
    if qualitative and set(qualitative) - {*QUALITATIVE_FLAGS, 'geometry_change_flag', 'notes'}:
        raise ValueError('Unknown qualitative evidence fields')
    actual = prepared['actual_value']
    row = {key: prepared.get(key) for key in
           ('case_id', 'surface_id', 'face_context', 'dimension', 'baseline_case_id', 'target_value', 'actual_value')}
    row.update(actual)
    row.update(energy_eV=energetics.get('E_slab_eV'), bulk_reference_id=energetics.get('bulk_reference_id'),
               Gamma_pair_meV_A2=energetics.get('Gamma_pair_meV_A2'),
               matched_reference_validated=energetics.get('matched_reference_validated') is True,
               comparison_protocol=prepared['comparison_protocol'],
               normal_termination=completed.get('normal_termination'),
               electronic_converged=completed.get('electronic_converged'),
               force_converged=completed.get('force_converged'), human_review_status='REQUIRED')
    if prepared['dimension'] == 'combined_verification':
        row['combined_selection_review'] = prepared['prior_selections']
    if alignment:
        if (alignment.get('surface_id') != prepared['surface_id']
                or alignment.get('static_result_id') != completed.get('static_result_id')
                or alignment.get('face_ids') != prepared['face_context']['face_ids']):
            raise ValueError('Electronic alignment belongs to a different SCF model or face context')
        for key in ('alignment_status', 'vacuum_upper_eV', 'vacuum_lower_eV',
                    'work_function_upper_eV', 'work_function_lower_eV'):
            row[key] = alignment.get(key)
    row.update(qualitative or {})
    return row
