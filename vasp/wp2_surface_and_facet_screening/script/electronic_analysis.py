"""Pure Stage 05 projections, face alignment and explicit partial-charge settings.

All energies retain the source SCF zero until a reviewed LVHAR plateau is supplied.
PDOS sums retain spin and states/eV per cell; PAW projections need not sum to TDOS.
No file reader, calculation, scientific promotion or input materialization runs here.

LPARD semantics: https://vasp.at/wiki/LPARD and https://vasp.at/wiki/KPUSE
(checked 2026-09-12). Selected KPUSE changes weights and requires an ISYM=-1 source.
"""

import numpy as np


LAYER_LABELS = ('lower_surface', 'lower_subsurface', 'interior',
                'upper_subsurface', 'upper_surface')


def _finite_number(value):
    """Treat absent scalar evidence as unavailable instead of inventing a value."""
    return (isinstance(value, (int, float, np.integer, np.floating))
            and not isinstance(value, (bool, np.bool_)) and bool(np.isfinite(value)))


def _group_indices(grouping, atom_count):
    """Validate the geometry partition, retaining explicitly unassigned atoms."""
    groups = {label: list(grouping['groups'].get(label, [])) for label in LAYER_LABELS}
    groups['unassigned'] = list(grouping.get('unassigned_site_indices', []))
    indices = [index for values in groups.values() for index in values]
    if any(type(index) is not int for index in indices) or sorted(indices) != list(range(atom_count)):
        raise ValueError('Layer groups must partition all site indices exactly once')
    return groups


def _identity(surface_id, static_result_id, grouping):
    """Preserve the supplied model and Stage 0 face identities without renaming."""
    if not surface_id or not static_result_id or grouping.get('surface_id') != surface_id:
        raise ValueError('Explicit matching surface_id and static_result_id are required')
    faces = grouping.get('face_ids', {})
    if set(faces) != {'upper', 'lower'} or not all(isinstance(x, str) and x for x in faces.values()):
        raise ValueError('Supply both original Stage 0 face IDs')
    if faces['upper'] == faces['lower']:
        raise ValueError('The two sides retain distinct face IDs even when equivalent')
    return dict(faces)


def aggregate_pdos(energies_eV, total_dos, site_orbital_dos, grouping, *, surface_id,
                   static_result_id, undercoordinated_indices=()):
    """Sum arrays without interpolation, normalization, spin folding or gap inference.

    ``total_dos`` maps spin labels to length-N energy arrays. ``site_orbital_dos``
    is a site-ordered list of {orbital: {spin: length-N array}} mappings. Energies
    must increase strictly; densities are nonnegative states/eV. Empty groups
    are None (unavailable), whereas an observed zero-density group stays zero.
    Face projections include that face's surface and subsurface groups.
    """
    energies = np.asarray(energies_eV, dtype=float)
    if energies.ndim != 1 or len(energies) < 2 or not np.isfinite(energies).all() or np.any(np.diff(energies) <= 0):
        raise ValueError('DOS energy grid must be finite and strictly increasing')
    if not site_orbital_dos:
        raise ValueError('Site projections are required for DOS aggregation')
    faces = _identity(surface_id, static_result_id, grouping)
    groups = _group_indices(grouping, len(site_orbital_dos))
    spins = tuple(total_dos)
    if not spins or any(spin not in ('up', 'down') for spin in spins):
        raise ValueError('Use explicit up/down spin labels, including up for ISPIN=1')

    def checked(values):
        if set(values) != set(spins):
            raise ValueError('Every projection must preserve all source spin channels')
        result = {spin: np.asarray(values[spin], dtype=float) for spin in spins}
        if any(x.shape != energies.shape or not np.isfinite(x).all() or np.any(x < 0) for x in result.values()):
            raise ValueError('DOS densities must match the grid and be finite/nonnegative')
        return result

    total = checked(total_dos)
    orbitals = {}
    sites = []
    for site in site_orbital_dos:
        if not site:
            raise ValueError('Missing site projections must not be filled with zero')
        site_sum = {spin: np.zeros_like(energies) for spin in spins}
        for orbital, values in site.items():
            channel = checked(values)
            orbital_sum = orbitals.setdefault(str(orbital), {spin: np.zeros_like(energies) for spin in spins})
            for spin in spins:
                site_sum[spin] += channel[spin]
                orbital_sum[spin] += channel[spin]
        sites.append(site_sum)

    def summed(indices):
        return ({spin: np.sum([sites[index][spin] for index in indices], axis=0) for spin in spins}
                if indices else None)

    undercoord = list(undercoordinated_indices)
    if len(set(undercoord)) != len(undercoord) or any(type(x) is not int or x < 0 or x >= len(sites) for x in undercoord):
        raise ValueError('Under-coordinated atom indices must be unique source site indices')
    return dict(surface_id=surface_id, static_result_id=static_result_id, face_ids=faces,
                energies_eV=energies.copy(), total_dos=total, orbital_pdos=orbitals,
                site_pdos=sites, layer_pdos={label: summed(indices) for label, indices in groups.items()},
                face_pdos={faces[side]: summed(groups[side+'_surface'] + groups[side+'_subsurface'])
                           for side in ('lower', 'upper')},
                undercoordinated_site_pdos={index: sites[index] for index in undercoord},
                projection_note='PAW-projected sums exclude unprojected weight; grouping is geometric, not proof of bulk-like states.')


def complete_dos_arrays(complete_dos):
    """Adapt an already parsed pymatgen CompleteDos without reading any files."""
    def spin_arrays(values):
        return {('up' if int(spin) == 1 else 'down'): np.asarray(array, dtype=float)
                for spin, array in values.items()}

    sites = [{str(orbital): spin_arrays(values) for orbital, values in complete_dos.pdos[site].items()}
             for site in complete_dos.structure]
    return dict(energies_eV=np.asarray(complete_dos.energies, dtype=float),
                total_dos=spin_arrays(complete_dos.densities), site_orbital_dos=sites)


def classify_state(site_weights, grouping, *, evidence, surface_threshold=0.6, interior_threshold=0.6):
    """Describe a state using reviewed projections; geometry alone cannot make it bulk-like.

    Fractions normalize only the supplied nonnegative site-projected weight.
    The adjustable 0.6 thresholds are working descriptive cutoffs, not physical
    criteria. Adequate projection coverage and source identity require evidence.
    """
    weights = np.asarray(site_weights, dtype=float)
    if weights.ndim != 1 or not len(weights) or not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError('State projection weights must be finite, nonnegative and site ordered')
    if not 0.5 < surface_threshold <= 1 or not 0.5 < interior_threshold <= 1:
        raise ValueError('Character thresholds must lie in (0.5, 1]')
    groups = _group_indices(grouping, len(weights))
    result = dict(character='UNRESOLVED', localisation=None, layer_fractions=None)
    if (weights.sum() <= 0 or evidence.get('projection_coverage_verified') is not True
            or not evidence.get('projection_review_reference') or not evidence.get('static_result_id')
            or evidence.get('projection_source_id') != evidence['static_result_id']):
        return dict(result, reason='Missing reviewed projection coverage or matching SCF source')
    fractions = {label: float(weights[indices].sum()/weights.sum()) for label, indices in groups.items()}
    result['layer_fractions'] = fractions
    if fractions['unassigned'] > 0 or grouping.get('status') == 'AMBIGUOUS':
        return dict(result, reason='Projection weight falls in ambiguous geometry groups')
    lower = fractions['lower_surface'] + fractions['lower_subsurface']
    upper = fractions['upper_surface'] + fractions['upper_subsurface']
    if lower + upper >= surface_threshold:
        side = 'lower' if lower >= surface_threshold else 'upper' if upper >= surface_threshold else 'both'
        result.update(character='SURFACE_LOCALIZED', localisation=side)
    elif (grouping.get('interior_available') is True and fractions['interior'] >= interior_threshold
          and evidence.get('bulk_like_interior_verified') is True and evidence.get('bulk_like_review_reference')):
        result.update(character='CANDIDATE_BULK_LIKE_INTERIOR', localisation='interior')
    else:
        result.update(character='FINITE_SLAB_STATE', localisation='distributed_or_unverified_interior')
    return dict(result, reason='Projection-based character; no automatic scientific approval')


def sampled_band_edges(eigenvalues, projections, *, static_result_id, full_occupation=1.0,
                       occupation_tolerance=1e-3, degeneracy_tolerance_eV=1e-5):
    """Collect sampled edge states, preserving every degenerate spin/k/band state.

    Eigenvalues: {spin: (nk, nb, [energy_eV, occupation])}; projections: matching
    {spin: (nk, nb, nsites, norbitals)}. XML full occupancy is 1 by default.
    Partial occupations never produce a forced semiconducting VBM/CBM.
    Returned band/k indices are one-based for explicit later selection.
    """
    if (not static_result_id or not eigenvalues or set(eigenvalues) != set(projections)
            or full_occupation <= 0 or not 0 < occupation_tolerance < full_occupation/2
            or degeneracy_tolerance_eV < 0):
        raise ValueError('Matching eigenvalues/projections, source and valid tolerances are required')
    occupied, empty = [], []
    partial = False
    atom_count = None
    for spin, raw in eigenvalues.items():
        values = np.asarray(raw, dtype=float)
        weights = np.asarray(projections[spin], dtype=float)
        if (values.ndim != 3 or values.shape[-1] != 2 or weights.ndim != 4
                or weights.shape[:2] != values.shape[:2] or not np.isfinite(values).all()
                or not np.isfinite(weights).all() or np.any(weights < 0)):
            raise ValueError('Invalid eigenvalue/projection shapes or values')
        if atom_count is not None and atom_count != weights.shape[2]:
            raise ValueError('Spin projections must use the same site ordering and count')
        atom_count = weights.shape[2]
        for kpoint, band in np.ndindex(values.shape[:2]):
            energy, occupation = values[kpoint, band]
            if occupation < -occupation_tolerance or occupation > full_occupation + occupation_tolerance:
                raise ValueError('Occupation lies outside the supplied convention')
            state = dict(energy_eV=float(energy), spin=str(spin), kpoint_index=kpoint+1,
                         band_index=band+1, site_weights=weights[kpoint, band].sum(axis=1).tolist())
            if occupation >= full_occupation-occupation_tolerance:
                occupied.append(state)
            elif occupation <= occupation_tolerance:
                empty.append(state)
            else:
                partial = True
    result = dict(static_result_id=static_result_id, sampled_integer_occupation_gap_eV=None, gap_classification='UNRESOLVED', VBM=None, CBM=None)
    if partial:
        return dict(result, gap_classification='METALLIC_LIKE_OCCUPATIONS_OR_UNRESOLVED')
    if not occupied or not empty:
        return result
    vbm, cbm = max(x['energy_eV'] for x in occupied), min(x['energy_eV'] for x in empty)
    if cbm <= vbm:
        return dict(result, gap_classification='NO_POSITIVE_SAMPLED_GAP')
    result.update(sampled_integer_occupation_gap_eV=cbm-vbm, gap_classification='FINITE_SAMPLED_GAP')
    for label, edge, states in [('VBM', vbm, occupied), ('CBM', cbm, empty)]:
        result[label] = dict(energy_eV=edge, source_id=static_result_id,
                             states=[x for x in states if abs(x['energy_eV']-edge) <= degeneracy_tolerance_eV])
    return result


def vacuum_reference(vacuum, *, grouping, static_result_id, fermi_eV, source_evidence, plateau_review):
    """Calculate face-specific work functions Vvac - EF only for reviewed plateaus.

    A candidate plateau from Stage 01 is insufficient by itself. Fermi and LVHAR
    potential must share the same converged self-consistent result and geometry.
    A missing face remains None; an available opposite face retains its own ID.
    """
    surface_id = grouping.get('surface_id')
    faces = _identity(surface_id, static_result_id, grouping)
    result = dict(surface_id=surface_id, static_result_id=static_result_id, face_ids=faces,
                  fermi_eV=None, alignment_status='INSUFFICIENT_DATA', faces={})
    for side in ('upper', 'lower'):
        result['vacuum_'+side+'_eV'] = None
        result['work_function_'+side+'_eV'] = None
        result['faces'][faces[side]] = dict(face_id=faces[side], side=side, vacuum_eV=None, work_function_eV=None)
    if (not _finite_number(fermi_eV) or any(source_evidence.get(key) is not True for key in
            ('self_consistent', 'electronic_converged', 'normal_termination', 'geometry_identity_verified'))
            or any(source_evidence.get(key) != static_result_id for key in ('fermi_source_id', 'potential_source_id'))
            or source_evidence.get('potential_kind') != 'LVHAR'):
        return dict(result, notes='Matching converged SCF, geometry and LVHAR provenance are required')
    result['fermi_eV'] = float(fermi_eV)
    available = 0
    for side in ('upper', 'lower'):
        candidate = vacuum.get(side, {})
        selected = candidate.get('selected_window') or {}
        reviewed = plateau_review.get(side) is True and bool(plateau_review.get('review_reference'))
        if candidate.get('status') != 'CANDIDATE_PLATEAU' or selected.get('status') != 'CANDIDATE_PLATEAU' or not reviewed:
            continue
        level = selected.get('mean_potential_eV')
        if not _finite_number(level):
            continue
        level = float(level)
        work_function = level - float(fermi_eV)
        result['vacuum_'+side+'_eV'] = level
        result['work_function_'+side+'_eV'] = work_function
        result['faces'][faces[side]].update(vacuum_eV=level, work_function_eV=work_function)
        available += 1
    result['alignment_status'] = ('FACE_ALIGNED' if available == 2 else 'PARTIAL_FACE_ALIGNMENT'
                                  if available == 1 else 'NO_REVIEWED_PLATEAU')
    return dict(result, notes='Work functions use each face vacuum and the same SCF Fermi energy; semiconductor Fermi placement requires interpretation.')


def _require_current_edges(edges):
    """Reject legacy positive-gap records rather than guessing their semantics."""
    if 'sampled_gap_eV' in edges or 'sampled_integer_occupation_gap_eV' not in edges:
        raise ValueError('Legacy or insufficient edge data: recompute the diagnostic; sampled_gap_eV is retired')


def bulk_like_edges(edges, alignment, grouping, *, evidence):
    """Reference sampled candidate interior edges to each face vacuum, Eedge - Vvac.

    Finite-slab/surface states remain classified but do not become bulk-like edges.
    All degenerate edge states must have reviewed interior character. Persistent
    fields, reconstruction, missing evidence and unrelated energy zeros fail closed.
    """
    _require_current_edges(edges)
    result = dict(alignment_status='INSUFFICIENT_DATA', bulk_like_interior_available=None,
                  VBM_character=None, CBM_character=None, VBM_localisation=None, CBM_localisation=None)
    for side in ('upper', 'lower'):
        for label in ('VBM', 'CBM'):
            result['bulk_like_'+label+'_'+side+'_eV'] = None
    if (grouping.get('interior_available') is not True or grouping.get('status') == 'AMBIGUOUS'):
        return dict(result, alignment_status='NO_DEFENSIBLE_INTERIOR', bulk_like_interior_available=False)
    if evidence.get('bulk_like_interior_verified') is True and evidence.get('bulk_like_review_reference'):
        result['bulk_like_interior_available'] = True
    if evidence.get('internal_field_flag') is not False or evidence.get('electronic_reconstruction_flag') is not False:
        return dict(result, alignment_status='FIELD_OR_RECONSTRUCTION_UNRESOLVED')
    source = edges.get('static_result_id')
    if (not source or source != alignment.get('static_result_id') or source != evidence.get('static_result_id')
            or source != evidence.get('projection_source_id') or evidence.get('self_consistent') is not True):
        return dict(result, alignment_status='INCOMPATIBLE_SCF_SOURCE')
    if alignment.get('surface_id') != grouping.get('surface_id') or alignment.get('face_ids') != grouping.get('face_ids'):
        return dict(result, alignment_status='INCOMPATIBLE_FACE_CONTEXT')
    if alignment.get('alignment_status') != 'FACE_ALIGNED':
        return dict(result, alignment_status='NO_USABLE_FACE_PLATEAUS')
    if any(not _finite_number(alignment.get('vacuum_'+side+'_eV')) for side in ('upper', 'lower')):
        return dict(result, alignment_status='NO_USABLE_FACE_PLATEAUS')
    if edges.get('gap_classification') != 'FINITE_SAMPLED_GAP':
        return dict(result, alignment_status='NO_RESOLVED_SEMICONDUCTING_EDGES')
    for label in ('VBM', 'CBM'):
        edge = edges.get(label) or {}
        if edge.get('source_id') != source or not edge.get('states') or not _finite_number(edge.get('energy_eV')):
            return dict(result, alignment_status='INCOMPATIBLE_OR_MISSING_EDGE_PROJECTIONS')
        characters = [classify_state(state['site_weights'], grouping, evidence=evidence) for state in edge['states']]
        result[label+'_character'] = ','.join(sorted(set(x['character'] for x in characters)))
        result[label+'_localisation'] = ','.join(sorted(set(x['localisation'] or 'unresolved' for x in characters)))
    if any(result[label+'_character'] != 'CANDIDATE_BULK_LIKE_INTERIOR' for label in ('VBM', 'CBM')):
        return dict(result, alignment_status='NO_REVIEWED_BULK_LIKE_EDGES')
    result['bulk_like_interior_available'] = True
    for side in ('upper', 'lower'):
        for label in ('VBM', 'CBM'):
            result['bulk_like_'+label+'_'+side+'_eV'] = edges[label]['energy_eV'] - alignment['vacuum_'+side+'_eV']
    return dict(result, alignment_status='CANDIDATE_BULK_LIKE_ALIGNED', notes='Eedge - Vvac in eV; negative values lie below vacuum. Human review remains required.')


def electronic_summary(alignment, edges, character, *, spin_mode, internal_field_flag=None):
    """Return two canonical electronic-schema rows, one for each stable face ID."""
    from result_schemas import schema_record

    _require_current_edges(edges)
    values = dict(surface_id=alignment['surface_id'], static_result_id=alignment['static_result_id'],
                  spin_mode=spin_mode, internal_field_flag=internal_field_flag,
                  fermi_eV=alignment.get('fermi_eV'), sampled_integer_occupation_gap_eV=edges.get('sampled_integer_occupation_gap_eV'),
                  gap_classification=edges.get('gap_classification'), alignment_status=character.get('alignment_status'),
                  notes=character.get('notes', alignment.get('notes')))
    for key in ('vacuum_upper_eV', 'vacuum_lower_eV', 'work_function_upper_eV', 'work_function_lower_eV'):
        values[key] = alignment.get(key)
    for key in ('bulk_like_interior_available', 'VBM_character', 'CBM_character', 'VBM_localisation', 'CBM_localisation'):
        values[key] = character.get(key)
    # An absent positive finding is unknown, not evidence that surface states are absent.
    values['surface_state_present'] = (True if any('SURFACE_LOCALIZED' in (character.get(label+'_character') or '')
                                                  for label in ('VBM', 'CBM')) else None)
    return [schema_record('electronic', dict(values, face_id=alignment['face_ids'][side])) for side in ('lower', 'upper')]


def partial_charge_template(*, approval, source_wavecar, band_indices, kpoint_selection, rationale,
                            separate_bands=True, separate_kpoints=True):
    """Return LPARD settings only, requiring identified converged orbitals and selection.

    Source metadata: identity, static_result_id, surface_id, converged=True,
    isym, ispin, encut_eV, nbands, nkpoints, noncollinear=False. One-based band
    indices are always explicit. kpoint_selection is {'mode': 'all'} or
    {'mode': 'indices', 'indices': [...]}; the latter requires source ISYM=-1.
    Spin-polarized PARCHG contains total and magnetization, not a single spin file.
    """
    from stage_templates import approved_context

    context = approved_context(approval)
    required = ('identity', 'static_result_id', 'surface_id', 'converged', 'isym', 'ispin',
                'encut_eV', 'nbands', 'nkpoints', 'noncollinear')
    if any(key not in source_wavecar for key in required) or not rationale or not str(rationale).strip():
        raise ValueError('Identified WAVECAR metadata and explicit state-selection rationale are required')
    if (source_wavecar['surface_id'] != context['surface_id'] or source_wavecar['converged'] is not True
            or source_wavecar['noncollinear'] is not False or not source_wavecar['identity']
            or not source_wavecar['static_result_id'] or source_wavecar['ispin'] not in (1, 2)
            or not _finite_number(source_wavecar['encut_eV']) or source_wavecar['encut_eV'] <= 0):
        raise ValueError('WAVECAR must be a matching converged collinear SCF source')
    if (type(source_wavecar['isym']) is not int or source_wavecar['isym'] not in (-1, 0, 1, 2, 3)
            or type(source_wavecar['nkpoints']) is not int or source_wavecar['nkpoints'] <= 0
            or type(separate_bands) is not bool or type(separate_kpoints) is not bool):
        raise ValueError('Valid source symmetry, k-point count and explicit output booleans are required')

    def selected_indices(indices, count, label):
        values = list(indices)
        if (type(count) is not int or count <= 0 or not values or len(set(values)) != len(values)
                or any(type(index) is not int or not 1 <= index <= count for index in values)):
            raise ValueError(label+' must be explicit unique one-based source indices')
        return values

    bands = selected_indices(band_indices, source_wavecar['nbands'], 'Bands')
    incar = dict(LPARD=True, IBAND=bands, ISYM=source_wavecar['isym'], ISPIN=source_wavecar['ispin'],
                 ENCUT=source_wavecar['encut_eV'], LSEPB=bool(separate_bands), LSEPK=bool(separate_kpoints),
                 LCHARG=False, LWAVE=False)
    if kpoint_selection.get('mode') == 'indices' and set(kpoint_selection) == {'mode', 'indices'}:
        if source_wavecar['isym'] != -1:
            raise ValueError('Selected KPUSE requires ISYM=-1 in the source SCF and partial-charge run')
        incar['KPUSE'] = selected_indices(kpoint_selection['indices'], source_wavecar['nkpoints'], 'K points')
    elif kpoint_selection != {'mode': 'all'}:
        raise ValueError('Explicit all or indices k-point selection is required')
    return dict(context=context, incar=incar, source_wavecar=dict(source_wavecar),
                kpoint_selection=dict(kpoint_selection), state_selection_rationale=rationale,
                preparation_status='SPECIFICATION_ONLY',
                notes='Preserve source geometry, PAW and k-point sampling; no guessed bands. PARCHG spin channels need total/magnetization interpretation.')
