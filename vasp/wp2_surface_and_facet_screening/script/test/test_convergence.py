"""One-axis preparation and evidence-based convergence; all writes use temp dirs."""

import copy
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from pymatgen.io.vasp import Poscar

import convergence_cases as prepare
from convergence_analysis import compare_cases, convergence_record
from result_schemas import COLUMNS, csv_text, schema_descriptor, schema_record
from surface_layers import group_surface_layers

ROOT = Path(__file__).resolve().parents[2]


class SequentialPreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        seed = ROOT/'structure/slab_candidates/beta/001/beta_001_t07'
        cls.slab = Poscar.from_file(seed/'POSCAR', check_for_potcar=False).structure
        cls.metadata = json.loads((seed/'metadata.json').read_text())
        cls.parent = Poscar.from_file(ROOT/'structure/bulk_parents/beta/POSCAR', check_for_potcar=False).structure
        cls.approval = dict(approved=True, surface_id='beta_001_t07',
                            face_ids={'upper':'beta_001_t07_upper', 'lower':'beta_001_t07_lower'},
                            review_reference='SYNTHETIC TEST ONLY; not an actual Stage 01 approval')
        cls.settings = dict(encut_eV=500, k_mesh=[10,10,1], mode='static', relaxed_region_mode='full', n_layers=1)

    def request(self, dimension, values=None, **options):
        args = dict(approval=self.approval, baseline_case_id='synthetic', baseline_settings=self.settings,
                    dimension=dimension, values=values)
        args.update(options)
        return prepare.prepare_dimension(self.slab, **args)

    def selections(self, **values):
        return {key: dict(surface_id=self.approval['surface_id'], case_id=f'test-{key}',
                          review_reference='SYNTHETIC reviewed selection', value=value) for key,value in values.items()}

    def test_stage01_kmesh_and_separate_encut_with_no_implicit_writes(self):
        with patch.object(Path, 'mkdir', side_effect=AssertionError('implicit write')), \
             patch.object(Path, 'write_text', side_effect=AssertionError('implicit write')):
            k = self.request('kpoints', [.20, .16, .12])
            e = self.request('encut', [500, 600])
        self.assertEqual([c['actual_value']['k_mesh'] for c in k], [[10,10,1], [12,12,1], [16,16,1]])
        self.assertTrue(all(c['actual_value']['encut_eV'] == 500 for c in k))
        self.assertTrue(all(c['actual_value']['k_mesh'] == [10,10,1] for c in e))
        self.assertAlmostEqual(k[0]['transformation']['reciprocal_lengths_Ainv'][0], 1.87072347, places=7)

    def test_no_matrix_no_unreviewed_next_stage_and_protocol_changes(self):
        bad = [('all', [1]), ('encut', {'encut':[500,600], 'vacuum':[20,25]}),
               ('kpoints', [{'spacing':.2, 'encut':600}]), ('vacuum', [25]),
               ('combined_verification', [20,25]), ('slab_thickness', [{'n_layers':2, 'vacuum_A':25}])]
        for dimension, values in bad:
            with self.subTest(dimension=dimension, values=values), self.assertRaises(ValueError):
                self.request(dimension, values)
        with self.assertRaises(ValueError):
            self.request('encut', [600], selections=self.selections(kpoints=[12,12,1]))
        with self.assertRaises(ValueError):
            self.request('kpoints', [.2], baseline_settings=dict(self.settings, vacuum=[20,25]))
        with self.assertRaises(ValueError):
            self.request('encut', [500], approval=dict(self.approval, approved=False))
        with self.assertRaises(ValueError):
            self.request('encut', [500], baseline_settings=dict(self.settings, n_layers=[2,4,6]))

    def test_vacuum_and_thickness_are_real_structure_changes(self):
        selections = self.selections(kpoints=[10,10,1], encut=500)
        vacuum = self.request('vacuum', [25,30], selections=selections)
        self.assertEqual(len(vacuum), 2)
        for case, target in zip(vacuum, [25,30]):
            self.assertAlmostEqual(case['actual_value']['vacuum_A'], target)
            np.testing.assert_allclose(case['structure'].cart_coords-self.slab.cart_coords,
                                       np.tile([0,0,(target-20)/2], (len(self.slab),1)), atol=1e-12)
        selections.update(self.selections(vacuum=20))
        thickness = self.request('slab_thickness', [{'n_layers':2}, {'n_layers':4}],
                                 selections=selections, parent=self.parent, seed_metadata=self.metadata)
        self.assertEqual([len(c['structure']) for c in thickness], [14,28])
        self.assertEqual([c['actual_value']['n_layers'] for c in thickness], [2,4])
        self.assertTrue(all(c['face_context']['face_ids'] == self.approval['face_ids'] for c in thickness))
        self.assertTrue(all(abs(c['actual_value']['vacuum_A']-20)<1e-10 for c in thickness))

    def test_single_combined_verification_and_fixed_region_gates(self):
        geometry = prepare._geometry(self.slab)
        selections = self.selections(kpoints=[10,10,1], encut=500, vacuum=20,
                                     slab_thickness={'n_layers':1}, relaxed_region='full')
        selections['slab_thickness']['actual_thickness_A'] = geometry['slab_thickness_A']
        cases = self.request('combined_verification', selections=selections)
        self.assertEqual(len(cases), 1)
        with self.assertRaises(ValueError):
            self.request('combined_verification', [dict(encut=600)], selections=selections)
        with self.assertRaises(ValueError):
            self.request('relaxed_region', ['central_fixed'], selections=selections,
                         baseline_settings=dict(self.settings, mode='relaxation'))
        grown = self.request('slab_thickness', [{'n_layers':4}], parent=self.parent, seed_metadata=self.metadata,
                             selections=self.selections(kpoints=[10,10,1], encut=500, vacuum=20))[0]
        slab = grown['structure']; grouping = group_surface_layers(slab, face_context=grown['face_context'])
        evidence = dict(bulk_like_established=True, reference='SYNTHETIC', criterion='SYNTHETIC matched coordination',
                        surface_id=self.approval['surface_id'], site_indices=grouping['groups']['interior'],
                        geometry_id=grouping['geometry_id'])
        selections['slab_thickness'].update(value={'n_layers':4}, actual_thickness_A=grown['actual_value']['slab_thickness_A'])
        selections.pop('relaxed_region')
        cases = prepare.prepare_dimension(slab, approval=self.approval, baseline_case_id='test-thick',
                                          baseline_settings=dict(self.settings, mode='relaxation', n_layers=4),
                                          dimension='relaxed_region', values=['full','central_fixed'],
                                          selections=selections, interior_evidence=evidence)
        self.assertTrue(np.all(cases[0]['structure'].site_properties['selective_dynamics']))
        self.assertTrue(np.any(np.logical_not(cases[1]['structure'].site_properties['selective_dynamics'])))

    def test_result_join_requires_case_and_face_source_identity(self):
        prepared = self.request('encut', [600])[0]
        completed = dict(case_id=prepared['case_id'], static_result_id='fixture-static',
                         input_settings_verified=True, normal_termination=True, electronic_converged=True)
        energy = dict(surface_id=prepared['surface_id'], face_context=prepared['face_context'],
                      E_slab_eV=-10, bulk_reference_id='fixture-bulk', Gamma_pair_meV_A2=1,
                      matched_reference_validated=True)
        record = convergence_record(prepared, completed, energy)
        self.assertTrue(record['matched_reference_validated'])
        self.assertEqual(record['comparison_protocol'], prepared['comparison_protocol'])
        self.assertEqual(record['encut_eV'], 600)
        with self.assertRaises(ValueError):
            convergence_record(prepared, dict(completed, case_id='unrelated'), energy)
        with self.assertRaises(ValueError):
            convergence_record(prepared, completed, dict(energy, face_context={}))

    def test_explicit_temp_materialization_and_no_overwrite(self):
        cases = self.request('encut', [600])
        with tempfile.TemporaryDirectory(dir='/private/tmp') as tmp:
            destination = Path(tmp)/'bounded'
            [path] = prepare.materialize_cases(cases, destination)
            self.assertEqual({p.name for p in path.iterdir()}, {'POSCAR','INCAR','KPOINTS','case.json'})
            self.assertEqual(json.loads((path/'case.json').read_text())['dimension'], 'encut')
            with self.assertRaises(FileExistsError):
                prepare.materialize_cases(cases, destination)


class ConvergenceEvidenceTests(unittest.TestCase):
    def records(self):
        baseline = dict(case_id='baseline', surface_id='test', face_context={'surface_id':'test','face_ids':{'upper':'test_upper','lower':'test_lower'}},
                        encut_eV=500, k_mesh=[4,4,1], vacuum_A=20, slab_thickness_A=30, relaxed_region_mode='full',
                        Gamma_pair_meV_A2=10., work_function_upper_eV=4., work_function_lower_eV=5.,
                        normal_termination=True, electronic_converged=True, matched_reference_validated=True,
                        alignment_status='FACE_ALIGNED', comparison_protocol={'incar':{'NSW':0},'energy_convention':'energy(sigma->0)'})
        case = dict(baseline, case_id='case', baseline_case_id='baseline', dimension='encut', encut_eV=600,
                    Gamma_pair_meV_A2=11., work_function_upper_eV=4.05)
        return baseline, case

    def test_numeric_pass_is_not_scientific_pass_and_missing_flags_stay_unknown(self):
        baseline, case = self.records()
        result = compare_cases(baseline, case, electronic_relevant=True)
        self.assertEqual(result['numerical_status'], 'NUMERIC_PASS')
        self.assertEqual(result['human_review_status'], 'REQUIRED')
        self.assertIsNone(result['reconstruction_flag'])
        self.assertIsNone(result['delta_vacuum_upper_eV'])
        self.assertAlmostEqual(result['delta_vacuum_reference_upper_eV'], .05)
        result = compare_cases(baseline, dict(case, Gamma_pair_meV_A2=11.01), electronic_relevant=False)
        self.assertEqual(result['numerical_status'], 'NOT_CONVERGED')

    def test_missing_ambiguous_qualitative_and_unmatched_protocol(self):
        baseline, case = self.records()
        for changes, expected in [({'reconstruction_flag':True},'HUMAN_REVIEW_REQUIRED'),
                                  ({'localisation_change_flag':True},'HUMAN_REVIEW_REQUIRED'),
                                  ({'normal_termination':False},'INSUFFICIENT_DATA'),
                                  ({'Gamma_pair_meV_A2':None},'INSUFFICIENT_DATA'),
                                  ({'Gamma_pair_meV_A2':float('nan')},'AMBIGUOUS'),
                                  ({'alignment_status':None},'AMBIGUOUS'),
                                  ({'vacuum_A':25},'AMBIGUOUS'),
                                  ({'comparison_protocol':{'incar':{'ISPIN':2}}},'AMBIGUOUS')]:
            with self.subTest(changes=changes):
                self.assertEqual(compare_cases(baseline, dict(case, **changes), electronic_relevant=True)['numerical_status'], expected)
        baseline['comparison_protocol']['incar']['NSW'] = 100
        case['comparison_protocol'] = copy.deepcopy(baseline['comparison_protocol'])
        self.assertEqual(compare_cases(baseline, case, electronic_relevant=False)['numerical_status'], 'INSUFFICIENT_DATA')

    def test_schemas_are_documented_ordered_and_missing_not_zero(self):
        for kind, columns in COLUMNS.items():
            descriptor = schema_descriptor(kind)
            self.assertEqual([f['name'] for f in descriptor['fields']], columns)
            self.assertTrue(all(f['description'] for f in descriptor['fields']))
            row = schema_record(kind, {'surface_id':'SYNTHETIC'})
            self.assertEqual(list(row), columns)
            self.assertIsNone(row['notes'])
            output = list(csv.DictReader(io.StringIO(csv_text(kind, [row]))))
            self.assertEqual(output[0]['notes'], '')
            self.assertEqual(len(csv_text(kind, []).splitlines()), 1)
            with self.assertRaises(ValueError):
                schema_record(kind, {'invented_field':1})


if __name__ == '__main__':
    unittest.main()
