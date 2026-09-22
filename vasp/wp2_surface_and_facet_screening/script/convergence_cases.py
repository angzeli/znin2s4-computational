"""Prepare ONE explicitly reviewed convergence dimension; no queue or engine API.

API: prepare_dimension(structure, approval=..., baseline_case_id=...,
 baseline_settings=..., dimension='vacuum', values=[25], selections=...).
CLI: python convergence_cases.py --dimension vacuum --request reviewed.json
Preview is the default. --output NEW_DIRECTORY explicitly writes a bounded set of
scientific inputs only. A request names exactly one structure, approval, baseline,
values and previously reviewed selections. No registry scan or --all exists.
"""

import argparse
import json
import math
import re
from copy import deepcopy
from pathlib import Path

import numpy as np
from pymatgen.io.vasp import Incar, Poscar

from stage_templates import surface_template
from surface_builders import change_vacuum, extend_thickness, slab_geometry
from surface_layers import propose_relaxed_region

DIMENSIONS = ('kpoints', 'encut', 'vacuum', 'slab_thickness', 'relaxed_region', 'combined_verification')
CANDIDATE_VALUES = {'kpoints': (0.20, 0.16, 0.12), 'encut': (500, 600), 'vacuum': (20, 25, 30)}
PREREQUISITES = {
    'kpoints': (), 'encut': (), 'vacuum': ('kpoints', 'encut'),
    'slab_thickness': ('kpoints', 'encut', 'vacuum'),
    'relaxed_region': ('kpoints', 'encut', 'vacuum', 'slab_thickness'),
    'combined_verification': ('kpoints', 'encut', 'vacuum', 'slab_thickness', 'relaxed_region'),
}
CALCULATION_LAYOUT = (
    '02_numerical_convergence/01_kpoints_and_encut', '02_numerical_convergence/02_vacuum',
    '02_numerical_convergence/03_slab_thickness', '02_numerical_convergence/04_relaxed_region',
    '02_numerical_convergence/05_combined_verification', '03_surface_relaxation/01_candidate_screening',
    '03_surface_relaxation/02_shortlist_validation', '03_surface_relaxation/03_final_relaxation',
    '04_static_scf/01_bulk_reference', '04_static_scf/02_surface_static',
    '05_surface_electronic_structure/01_dos_pdos', '05_surface_electronic_structure/02_band_edge_charge',
)
RESULT_LAYOUT = ('02_numerical_convergence', '03_surface_relaxation', '04_surface_energetics',
                 '05_surface_electronic_structure', '06_final_selection')


def ensure_layout(wp2_root):
    """Explicitly create only empty skeleton directories; never scientific inputs."""
    for category, names in (('calculation', CALCULATION_LAYOUT), ('results', RESULT_LAYOUT)):
        for name in names:
            (Path(wp2_root) / category / name).mkdir(parents=True, exist_ok=True)


def mesh_record(structure, spacing_Ainv):
    """Match Stage 01's Gamma-centered ceil(|b|/spacing), b includes 2*pi."""
    if isinstance(spacing_Ainv, bool) or not isinstance(spacing_Ainv, (int, float)) or not math.isfinite(spacing_Ainv) or spacing_Ainv <= 0:
        raise ValueError('One finite positive k spacing is required')
    lengths = np.linalg.norm(structure.lattice.reciprocal_lattice.matrix, axis=1)
    mesh = [int(math.ceil(lengths[i] / spacing_Ainv)) for i in (0, 1)] + [1]
    return dict(target_spacing_Ainv=spacing_Ainv, reciprocal_lengths_Ainv=lengths.tolist(), k_mesh=mesh,
                actual_spacing_Ainv=[float(lengths[i] / mesh[i]) for i in (0, 1)],
                k_spacing_max_Ainv=max(float(lengths[i] / mesh[i]) for i in (0, 1)), centering='Gamma')


def _geometry(structure):
    measured = slab_geometry(structure)
    normal = np.cross(*structure.lattice.matrix[:2]); normal /= np.linalg.norm(normal)
    if not np.allclose(normal, [0, 0, 1], atol=1e-8):
        raise ValueError('The Stage 0 +z normal convention is required')
    return {key: measured[key] for key in ('slab_thickness_A', 'vacuum_A')}


def _selection_checks(approval, dimension, selections, settings, geometry):
    if set(selections) - set(DIMENSIONS[:-1]):
        raise ValueError('Unknown selected convergence dimension')
    for prior in set(PREREQUISITES[dimension]) | (set(selections) - {dimension}):
        record = selections.get(prior, {})
        if (record.get('surface_id') != approval['surface_id'] or not record.get('case_id')
                or not record.get('review_reference') or 'value' not in record):
            raise ValueError(f'Explicit reviewed {prior} selection for this surface is required')
        selected = record['value']
        if prior == 'kpoints' and list(selected) != list(settings['k_mesh']):
            raise ValueError('Baseline k mesh differs from the reviewed selection')
        if prior == 'encut' and selected != settings['encut_eV']:
            raise ValueError('Baseline ENCUT differs from the reviewed selection')
        if prior == 'vacuum' and not math.isclose(selected, geometry['vacuum_A'], abs_tol=1e-7):
            raise ValueError('Baseline actual vacuum differs from the reviewed selection')
        if prior == 'slab_thickness':
            # Record the REALIZED geometry of the selected complete-repeat model,
            # so an approximate target cannot masquerade as the selected thickness.
            actual = record.get('actual_thickness_A')
            if actual is None or not math.isclose(actual, geometry['slab_thickness_A'], abs_tol=1e-7):
                raise ValueError('Baseline thickness does not match the reviewed realized geometry')
        if prior == 'relaxed_region' and selected != settings['relaxed_region_mode']:
            raise ValueError('Baseline constraint convention differs from the reviewed selection')


def prepare_dimension(structure, *, approval, baseline_case_id, baseline_settings, dimension,
                      values=None, selections=None, parent=None, seed_metadata=None,
                      interior_evidence=None):
    """Return at most three in-memory cases along ONE axis, or ONE verification.

    Approval: approved=True, surface_id, face_ids={upper,lower}, review_reference.
    Baseline: explicit encut_eV, k_mesh, mode, relaxed_region_mode; n_layers when
    verified. Prior selections: dimension -> {surface_id,case_id,review_reference,
    value}; thickness also needs actual_thickness_A. kpoints/ENCUT can be reviewed
    separately; all later dimensions enforce the preceding selections.

    Thickness seeds must match Stage 0 provenance; relaxed distortions are never
    copied into new bulk repeats. Combined verification consumes an ALREADY built,
    reviewed thickness/constraint model, not a new multidimensional search.
    """
    if dimension not in DIMENSIONS or not re.fullmatch(r'[A-Za-z0-9_.-]+', baseline_case_id or ''):
        raise ValueError('One valid dimension and explicit safe baseline_case_id required')
    settings = deepcopy(baseline_settings)
    if settings.get('mode') not in ('static', 'relaxation') or settings.get('relaxed_region_mode') not in ('full', 'central_fixed'):
        raise ValueError('Explicit baseline static/relaxation and constraint conventions required')
    allowed = {'encut_eV', 'k_mesh', 'mode', 'relaxed_region_mode', 'n_layers', 'template_options',
               'selected_by_review', 'integration_compatibility_id', 'integration_review_reference'}
    if set(settings) - allowed:
        raise ValueError('Unknown baseline settings; independent variations cannot hide in the baseline')
    layers = settings.get('n_layers')
    if layers is not None and (type(layers) is not int or layers < 1):
        raise ValueError('Verified n_layers is one positive integer or unavailable, never a varied list')
    options = settings.get('template_options', {})
    if not isinstance(options, dict) or set(options) - {'dipole_correction', 'magmom', 'smearing', 'paw_identities'}:
        raise ValueError('Only explicit baseline dipole/spin/smearing/PAW options are supported')
    # The template validates scalar ENCUT, one mesh, and explicit surface approval.
    surface_template(approval, settings, mode=settings['mode'], atom_count=len(structure), **options)
    geometry = _geometry(structure)
    selections = deepcopy(selections or {})
    _selection_checks(approval, dimension, selections, settings, geometry)
    if dimension == 'combined_verification':
        if values not in (None, []):
            raise ValueError('Combined verification consumes selections; it accepts no varied values')
        requested = [None]
    else:
        if not isinstance(values, (list, tuple)) or not 1 <= len(values) <= 3:
            raise ValueError('Explicitly request 1-3 values of ONE dimension; no matrix mappings')
        requested = list(values)
        if len({json.dumps(value, sort_keys=True) for value in requested}) != len(requested):
            raise ValueError('Duplicate values are not separate adaptive cases')
    if dimension == 'relaxed_region' and settings['mode'] != 'relaxation':
        raise ValueError('Constraint sensitivity requires an explicit relaxation baseline, not a static')
    context = {key: deepcopy(approval[key]) for key in ('surface_id', 'face_ids')}
    cases = []
    for index, value in enumerate(requested, 1):
        slab = structure.copy(); numerical = deepcopy(settings); transform = {}
        if dimension == 'kpoints':
            transform = mesh_record(slab, value); numerical['k_mesh'] = transform['k_mesh']
        elif dimension == 'encut':
            if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
                raise ValueError('One finite positive ENCUT per value; no nested variation mapping')
            numerical['encut_eV'] = value
        elif dimension == 'vacuum':
            if isinstance(value, bool) or not isinstance(value, (float, int)):
                raise ValueError('One scalar vacuum per value')
            slab, transform = change_vacuum(slab, value, face_context=context)
        elif dimension == 'slab_thickness':
            allowed = {'n_layers', 'oriented_repeats', 'target_thickness_A'}
            if not isinstance(value, dict) or len(value) != 1 or not set(value) <= allowed:
                raise ValueError('Thickness values specify exactly one complete-repeat target')
            if parent is None or seed_metadata is None or seed_metadata.get('slab_id') != approval['surface_id']:
                raise ValueError('Explicit same-surface Stage 0 seed provenance and copied parent required')
            slab, transform = extend_thickness(parent, slab, seed_metadata,
                                               vacuum_A=geometry['vacuum_A'], **value)
            numerical['n_layers'] = transform.get('n_layers', transform.get('complete_layer_count'))
        elif dimension == 'relaxed_region':
            if value not in ('full', 'central_fixed'):
                raise ValueError('Explicit full or central_fixed treatment required')
            numerical['relaxed_region_mode'] = value
        # Refresh the geometry-derived mask after any thickness change; never
        # replicate a list of fixed atom indices into a different slab.
        if dimension == 'relaxed_region' or numerical['relaxed_region_mode'] == 'central_fixed':
            mask, region = propose_relaxed_region(slab, mode=numerical['relaxed_region_mode'],
                                                 explicitly_requested=True, interior_evidence=interior_evidence,
                                                 face_context=context)
            slab.add_site_property('selective_dynamics', mask)
            transform['relaxed_region'] = region
        elif 'selective_dynamics' in slab.site_properties:
            if not np.all(slab.site_properties['selective_dynamics']):
                raise ValueError('Baseline full-relaxation label conflicts with fixed input atoms')
        template = surface_template(approval, numerical, mode=numerical['mode'], atom_count=len(slab), **options)
        actual = _geometry(slab)
        lengths = np.linalg.norm(slab.lattice.reciprocal_lattice.matrix, axis=1)
        actual.update(encut_eV=numerical['encut_eV'], k_mesh=list(numerical['k_mesh']),
                      k_spacing_max_Ainv=max(float(lengths[i] / numerical['k_mesh'][i]) for i in (0, 1)),
                      n_layers=numerical.get('n_layers'), relaxed_region_mode=numerical['relaxed_region_mode'])
        cases.append(dict(case_id=f'{baseline_case_id}--{dimension}--{index:02d}',
                          surface_id=approval['surface_id'], face_context=context, dimension=dimension,
                          baseline_case_id=baseline_case_id, target_value=value, actual_value=actual,
                          structure=slab, template=template, transformation=transform,
                          comparison_protocol=comparison_protocol(template),
                          approval=deepcopy(approval), prior_selections=selections,
                          status='PREPARED_SPECIFICATION_ONLY', human_review_status='REQUIRED'))
    if dimension == 'kpoints' and len({tuple(c['actual_value']['k_mesh']) for c in cases}) != len(cases):
        raise ValueError('Requested k spacings map to duplicate meshes; choose a bounded distinct comparison')
    return cases


def comparison_protocol(template):
    """Scientific settings held fixed in a one-dimension comparison."""
    incar = template['incar']
    keys = ('GGA', 'IVDW', 'PREC', 'LASPH', 'LREAL', 'ISYM', 'ISPIN', 'ISMEAR', 'SIGMA',
            'LDIPOL', 'IDIPOL', 'DIPOL', 'EDIFF', 'EDIFFG', 'IBRION', 'ISIF', 'NSW', 'ALGO')
    return dict(incar={key: incar.get(key) for key in keys},
                energy_convention=template['energy_convention'], vasp_version=template['vasp_version'],
                paw_identities=template.get('paw_identities'))


def case_record(case):
    """JSON-compatible preparation metadata; no fabricated energy/result fields."""
    return {key: value.as_dict() if key == 'structure' else value for key, value in case.items()}


def materialize_cases(cases, output_directory):
    """Explicit future action: write public inputs into an entirely NEW directory.

    This is not a launcher, queue command, or licensed-input assembler. No POTCAR,
    job state, next-stage folder, or registry is generated. All validation and
    serialization complete before the first write; existing paths are refused.
    """
    if not cases or len(cases) > 3 or len({c['dimension'] for c in cases}) != 1 or len({c['surface_id'] for c in cases}) != 1:
        raise ValueError('One surface and one bounded convergence dimension required')
    if cases[0]['dimension'] == 'combined_verification' and len(cases) != 1:
        raise ValueError('Only one combined verification is permitted')
    files = {}
    for case in cases:
        identity = case['case_id']
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', identity) or identity in files:
            raise ValueError('Unique safe case identities required')
        mesh = case['template']['k_mesh']
        files[identity] = {
            'POSCAR': Poscar(case['structure'], comment=identity, sort_structure=False).get_str(significant_figures=16),
            'INCAR': str(Incar(case['template']['incar'])),
            'KPOINTS': f'{identity}; b includes 2*pi\n0\nGamma\n'+ ' '.join(map(str, mesh))+'\n0 0 0\n',
            'case.json': json.dumps(case_record(case), indent=2, sort_keys=True, allow_nan=False)+'\n',
        }
    root = Path(output_directory)
    root.mkdir(parents=True, exist_ok=False)
    for identity, inputs in files.items():
        directory = root / identity; directory.mkdir()
        for name, content in inputs.items():
            (directory / name).write_text(content)
    return [root / identity for identity in files]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension', required=True, choices=DIMENSIONS)
    parser.add_argument('--request', type=Path, required=True, help='Explicit one-surface reviewed JSON request')
    parser.add_argument('--output', type=Path, help='Optional NEW directory; default is an in-memory preview')
    args = parser.parse_args()
    request = json.loads(args.request.read_text())
    slab = Poscar.from_file(request.pop('structure_path'), check_for_potcar=False).structure
    if 'parent_path' in request:
        request['parent'] = Poscar.from_file(request.pop('parent_path'), check_for_potcar=False).structure
    if 'seed_metadata_path' in request:
        request['seed_metadata'] = json.loads(Path(request.pop('seed_metadata_path')).read_text())
    cases = prepare_dimension(slab, dimension=args.dimension, **request)
    if args.output:
        materialize_cases(cases, args.output)
    else:
        print(json.dumps([case_record(case) for case in cases], indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
