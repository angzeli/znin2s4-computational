"""Synthetic Stage 05 arrays and pure partial-charge specifications; no engine/files."""

import copy
import unittest
from unittest.mock import patch

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.electronic_structure.dos import CompleteDos, Dos

import electronic_analysis as electronic
from electrostatics import vacuum_windows


def groups(interior=True):
    labels = list(electronic.LAYER_LABELS) if interior else ['lower_surface', 'upper_surface']
    return dict(surface_id='fixture-surface', face_ids={'lower': 'original-lower', 'upper': 'original-upper'},
                groups={label: [i for i, item in enumerate(labels) if item == label] for label in electronic.LAYER_LABELS},
                unassigned_site_indices=[], site_labels=labels, interior_available=interior,
                status='GEOMETRY_GROUPED' if interior else 'NO_BULK_LIKE_INTERIOR')


def evidence():
    return dict(static_result_id='fixture-scf', projection_source_id='fixture-scf', self_consistent=True,
                projection_coverage_verified=True, projection_review_reference='synthetic projection review',
                bulk_like_interior_verified=True, bulk_like_review_reference='synthetic interior review',
                internal_field_flag=False, electronic_reconstruction_flag=False)


def alignment():
    z = np.linspace(0, 30, 601)
    vacuum = vacuum_windows(z, np.where(z < 15, 3., 4.), 10, 20, 30)
    source = dict(self_consistent=True, electronic_converged=True, normal_termination=True,
                  geometry_identity_verified=True, fermi_source_id='fixture-scf',
                  potential_source_id='fixture-scf', potential_kind='LVHAR')
    return electronic.vacuum_reference(vacuum, grouping=groups(), static_result_id='fixture-scf',
                                       fermi_eV=1., source_evidence=source,
                                       plateau_review={'lower': True, 'upper': True, 'review_reference': 'synthetic review'})


def edges():
    eigenvalues = {'up': np.array([[[-1., 1.], [1., 0.]], [[-.5, 1.], [1.5, 0.]]])}
    weights = np.zeros((2, 2, 5, 1))
    weights[:, :, 2, 0] = 1.
    return electronic.sampled_band_edges(eigenvalues, {'up': weights}, static_result_id='fixture-scf')


class ElectronicAnalysisTests(unittest.TestCase):
    def test_dos_sums_preserve_spin_faces_and_unprojected_total(self):
        energy = np.array([-1., 0., 1.])
        sites = [{'s': {'up': np.full(3, i+1.), 'down': np.full(3, (i+1.)/2)}} for i in range(5)]
        total = {'up': np.full(3, 20.), 'down': np.full(3, 10.)}
        result = electronic.aggregate_pdos(energy, total, sites, groups(), surface_id='fixture-surface',
                                           static_result_id='fixture-scf', undercoordinated_indices=[0, 4])
        np.testing.assert_equal(result['orbital_pdos']['s']['up'], [15., 15., 15.])
        np.testing.assert_equal(result['total_dos']['up'], [20., 20., 20.])
        np.testing.assert_equal(result['face_pdos']['original-lower']['down'], [1.5, 1.5, 1.5])
        self.assertEqual(set(result['undercoordinated_site_pdos']), {0, 4})
        thin = electronic.aggregate_pdos(energy, total, sites[:2], groups(False), surface_id='fixture-surface',
                                         static_result_id='fixture-scf')
        self.assertIsNone(thin['layer_pdos']['interior'])
        invalid = copy.deepcopy(sites)
        invalid[1]['s'].pop('down')
        with self.assertRaises(ValueError):
            electronic.aggregate_pdos(energy, total, invalid, groups(), surface_id='fixture-surface', static_result_id='fixture-scf')
        ambiguous = groups()
        ambiguous['groups']['interior'] = []
        ambiguous['unassigned_site_indices'] = [2]
        result = electronic.aggregate_pdos(energy, total, sites, ambiguous, surface_id='fixture-surface', static_result_id='fixture-scf')
        self.assertIsNotNone(result['layer_pdos']['unassigned'])

    def test_complete_dos_adapter_retains_order_orbitals_and_spin(self):
        structure = Structure(Lattice.cubic(10), ['S'], [[0, 0, 0]])
        density = {Spin.up: np.array([1., 2., 3.]), Spin.down: np.array([.5, 1., 1.5])}
        complete = CompleteDos(structure, Dos(0., [-1., 0., 1.], density), {structure[0]: {Orbital.s: density}})
        result = electronic.complete_dos_arrays(complete)
        np.testing.assert_equal(result['total_dos']['down'], density[Spin.down])
        self.assertEqual(list(result['site_orbital_dos'][0]), ['s'])

    def test_state_character_needs_projection_and_interior_evidence(self):
        surface = electronic.classify_state([.8, .1, .05, .03, .02], groups(), evidence=evidence())
        self.assertEqual(surface['character'], 'SURFACE_LOCALIZED')
        self.assertEqual(surface['localisation'], 'lower')
        interior = electronic.classify_state([.025, .025, .9, .025, .025], groups(), evidence=evidence())
        self.assertEqual(interior['character'], 'CANDIDATE_BULK_LIKE_INTERIOR')
        unverified = evidence(); unverified.pop('bulk_like_review_reference')
        self.assertEqual(electronic.classify_state([0, 0, 1, 0, 0], groups(), evidence=unverified)['character'], 'FINITE_SLAB_STATE')
        unverified.pop('projection_review_reference')
        self.assertEqual(electronic.classify_state([0, 0, 1, 0, 0], groups(), evidence=unverified)['character'], 'UNRESOLVED')

    def test_face_workfunctions_and_plateau_failures(self):
        result = alignment()
        self.assertEqual(result['work_function_lower_eV'], 2.)
        self.assertEqual(result['work_function_upper_eV'], 3.)
        self.assertEqual(result['faces']['original-upper']['vacuum_eV'], 4.)
        z = np.linspace(0, 30, 601)
        source = dict(self_consistent=True, electronic_converged=True, normal_termination=True,
                      geometry_identity_verified=True, fermi_source_id='fixture-scf',
                      potential_source_id='fixture-scf', potential_kind='LVHAR')
        args = dict(grouping=groups(), static_result_id='fixture-scf', fermi_eV=1., source_evidence=source,
                    plateau_review={'lower': True, 'upper': True, 'review_reference': 'synthetic'})
        for potential in [.08*z, np.where(z < 15, 3., 4.)]:
            vacuum = vacuum_windows(z, potential, 10, 20, 30)
            if np.allclose(potential, .08*z):
                self.assertEqual(electronic.vacuum_reference(vacuum, **args)['alignment_status'], 'NO_REVIEWED_PLATEAU')
            else:
                args['plateau_review'] = {}
                self.assertIsNone(electronic.vacuum_reference(vacuum, **args)['vacuum_upper_eV'])
        potential = np.where(z < 15, 3., 4.); potential[(z >= 4.5) & (z < 15)] += 1
        args['plateau_review'] = {'lower': True, 'upper': True, 'review_reference': 'synthetic'}
        ambiguous = electronic.vacuum_reference(vacuum_windows(z, potential, 10, 20, 30), **args)
        self.assertEqual(ambiguous['alignment_status'], 'PARTIAL_FACE_ALIGNMENT')
        self.assertIsNone(ambiguous['work_function_lower_eV'])
        source['fermi_source_id'] = 'unrelated-bulk'
        self.assertEqual(electronic.vacuum_reference(vacuum_windows(z, potential, 10, 20, 30), **args)['alignment_status'], 'INSUFFICIENT_DATA')
        args['fermi_eV'] = None
        self.assertIsNone(electronic.vacuum_reference(vacuum_windows(z, potential, 10, 20, 30), **args)['fermi_eV'])

    def test_sampled_edges_preserve_indices_and_partial_occupation_limit(self):
        result = edges()
        self.assertEqual(result['sampled_integer_occupation_gap_eV'], 1.5)
        self.assertEqual(result['VBM']['states'][0]['kpoint_index'], 2)
        self.assertEqual(result['VBM']['states'][0]['band_index'], 1)
        eigenvalues = {'up': np.array([[[-1., .9], [1., .1]]])}
        projections = {'up': np.ones((1, 2, 5, 1))}
        result = electronic.sampled_band_edges(eigenvalues, projections, static_result_id='fixture-scf')
        self.assertIsNone(result['VBM'])
        self.assertEqual(result['gap_classification'], 'METALLIC_LIKE_OCCUPATIONS_OR_UNRESOLVED')

    def test_bulk_like_edges_require_all_alignment_and_state_evidence(self):
        result = electronic.bulk_like_edges(edges(), alignment(), groups(), evidence=evidence())
        self.assertEqual(result['alignment_status'], 'CANDIDATE_BULK_LIKE_ALIGNED')
        self.assertEqual(result['bulk_like_VBM_upper_eV'], -4.5)
        self.assertEqual(result['bulk_like_CBM_lower_eV'], -2.)
        self.assertEqual(electronic.bulk_like_edges(edges(), alignment(), groups(False), evidence=evidence())['alignment_status'], 'NO_DEFENSIBLE_INTERIOR')
        for field in ('internal_field_flag', 'electronic_reconstruction_flag'):
            for value in (True, None):
                bad = evidence(); bad[field] = value
                self.assertEqual(electronic.bulk_like_edges(edges(), alignment(), groups(), evidence=bad)['alignment_status'], 'FIELD_OR_RECONSTRUCTION_UNRESOLVED')
        wrong = edges(); wrong['VBM']['source_id'] = 'unrelated-bulk'
        self.assertEqual(electronic.bulk_like_edges(wrong, alignment(), groups(), evidence=evidence())['alignment_status'], 'INCOMPATIBLE_OR_MISSING_EDGE_PROJECTIONS')
        surface = edges(); surface['VBM']['states'][0]['site_weights'] = [1, 0, 0, 0, 0]
        limited = electronic.bulk_like_edges(surface, alignment(), groups(), evidence=evidence())
        self.assertEqual(limited['VBM_character'], 'SURFACE_LOCALIZED')
        self.assertIsNone(limited['bulk_like_VBM_upper_eV'])
        self.assertTrue(limited['bulk_like_interior_available'])
        degenerate = edges()
        degenerate['VBM']['states'].append(dict(degenerate['VBM']['states'][0], site_weights=[1, 0, 0, 0, 0]))
        self.assertEqual(electronic.bulk_like_edges(degenerate, alignment(), groups(), evidence=evidence())['alignment_status'], 'NO_REVIEWED_BULK_LIKE_EDGES')

    def test_partial_charge_requires_explicit_source_bands_and_kpoints(self):
        approval = dict(approved=True, surface_id='fixture-surface', face_ids=groups()['face_ids'], review_reference='synthetic review')
        source = dict(identity='fixture-wavecar-identity', static_result_id='fixture-scf', surface_id='fixture-surface',
                      converged=True, isym=0, ispin=1, encut_eV=500., nbands=10, nkpoints=4, noncollinear=False)
        args = dict(approval=approval, source_wavecar=source, band_indices=[2, 3],
                    kpoint_selection={'mode': 'all'}, rationale='Inspect explicitly identified test states')
        with patch('pathlib.Path.write_text', side_effect=AssertionError('unexpected write')), \
             patch('pathlib.Path.write_bytes', side_effect=AssertionError('unexpected write')):
            result = electronic.partial_charge_template(**args)
        self.assertTrue(result['incar']['LPARD'])
        self.assertNotIn('KPUSE', result['incar'])
        self.assertEqual(result['incar']['IBAND'], [2, 3])
        args['kpoint_selection'] = {'mode': 'indices', 'indices': [1, 3]}
        with self.assertRaisesRegex(ValueError, 'ISYM=-1'):
            electronic.partial_charge_template(**args)
        source['isym'] = -1
        self.assertEqual(electronic.partial_charge_template(**args)['incar']['KPUSE'], [1, 3])
        for invalid in ([], [0], [11], [1, 1]):
            args['band_indices'] = invalid
            with self.assertRaises(ValueError):
                electronic.partial_charge_template(**args)
        args['band_indices'] = [2]
        source['noncollinear'] = True
        with self.assertRaises(ValueError):
            electronic.partial_charge_template(**args)

    def test_summary_uses_canonical_columns_and_never_fabricates_absence(self):
        result = electronic.bulk_like_edges(edges(), alignment(), groups(), evidence=evidence())
        rows = electronic.electronic_summary(alignment(), edges(), result, spin_mode='nonmagnetic', internal_field_flag=False)
        self.assertEqual([row['face_id'] for row in rows], ['original-lower', 'original-upper'])
        self.assertIsNone(rows[0]['surface_state_present'])
        self.assertTrue(rows[0]['bulk_like_interior_available'])


if __name__ == '__main__':
    unittest.main()
