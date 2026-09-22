"""Authentic five-pilot regression and explicit unsupported-data boundaries."""

import copy
import csv
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np

import analyse_runtime_pilots as pilot
import electronic_analysis as electronic
import pilot_frontier as frontier
from electrostatics import vacuum_window_sensitivity
from result_schemas import schema_record


FIXTURES = Path(__file__).parents[1]/'fixtures'


class PilotCorrectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURES/'stage01_pilot_evidence.npz', allow_pickle=False) as source:
            cls.data = {key: source[key] for key in source.files}
        with (FIXTURES/'stage01_pilot_expectations.csv').open() as stream:
            cls.expected = {r['pilot']: r for r in csv.DictReader(stream)}

    def arguments(self, name):
        def item(key):
            return self.data[name+'_'+key]
        return dict(nelect=float(item('nelect')), ispin=int(item('ispin')), nbands=int(item('nbands')),
                    kpoints=item('kpoints'), weights=item('weights'), noncollinear=bool(item('noncollinear')),
                    soc=bool(item('soc')), occupation_convention=frontier.OCCUPATION_CONVENTION)

    def states(self, name, **overrides):
        d = self.data
        args = dict(normal=bool(d[name+'_normal']), electronic_converged=True, timeout=False, fatal=False,
                    nsw=int(d[name+'_nsw']), ibrion=int(d[name+'_ibrion']), ediffg=float(d[name+'_ediffg']),
                    forces=d[name+'_forces'], ionic_marker=bool(d[name+'_ionic_marker']))
        args.update(overrides)
        return pilot.result_states(**args)

    def test_five_audit_frontiers_and_raw_occupations(self):
        for name, expected in self.expected.items():
            with self.subTest(pilot=name):
                values = self.data[name+'_eigenval']
                before = values.copy()
                result = frontier.electronic_diagnostic({1: values}, **self.arguments(name))
                self.assertEqual(result['validation_status'], 'VALIDATED')
                self.assertAlmostEqual(result['frontier_band_separation_eV'], float(expected['expected_frontier_separation_eV']), delta=1e-6)
                self.assertEqual(result['frontier_interpretation'], expected['expected_interpretation'])
                self.assertEqual(result['frontier_valence_band_index'], 31 if name.startswith('A') else 124)
                self.assertNotIn('sampled_gap_eV', result)
                self.assertEqual(result['partial_occupation_present'], name != 'B1')
                np.testing.assert_array_equal(values, before)
                if name in ('A0','A1','B0'):
                    self.assertLess(result['frontier_band_separation_eV'], 0)
        a2 = self.data['A2_eigenval'][0, 30:32]
        np.testing.assert_allclose(a2[:, 0], [-1.022448, -1.022438], atol=5e-7, rtol=0)
        np.testing.assert_allclose(a2[:, 1], [.498501, .498390], atol=5e-7, rtol=0)
        b1 = frontier.electronic_diagnostic({1:self.data['B1_eigenval']}, **self.arguments('B1'))
        self.assertAlmostEqual(b1['minimum_direct_frontier_separation_eV'], .774108, delta=1e-6)
        self.assertFalse(self.states('B1')['force_converged'])

    def test_xml_precision_order_and_electron_counts(self):
        for name in self.expected:
            d = self.data
            check = frontier.crosscheck_eigenval_xml(d[name+'_eigenval'], d[name+'_xml'],
                d[name+'_kpoints'], d[name+'_xml_kpoints'], d[name+'_weights'], d[name+'_xml_weights'],
                float(d[name+'_nelect']))
            self.assertEqual(check['eigenval_xml_status'], 'VALIDATED')
            self.assertLessEqual(check['eigenval_xml_max_difference'], 5.1e-5)
            result = frontier.electronic_diagnostic({1:d[name+'_eigenval']}, **self.arguments(name))
            self.assertLess(abs(result['nelect_check']-float(d[name+'_nelect'])), result['electron_count_tolerance'])
        args = [self.data['A2_'+key].copy() for key in ('eigenval','xml','kpoints','xml_kpoints','weights','xml_weights')]
        for position in (1, 3, 5):
            changed = copy.deepcopy(args)
            changed[position][0] += .01
            with self.assertRaises(ValueError):
                frontier.crosscheck_eigenval_xml(*changed, 62.)

    def test_unsupported_metadata_and_occupation_scales_fail_closed(self):
        base = self.arguments('A2')
        for change in [dict(ispin=2), dict(soc=True), dict(noncollinear=True), dict(noncollinear=None),
                       dict(nelect=63.), dict(nelect=62.000002), dict(nelect=float('nan')),
                       dict(nbands=31), dict(occupation_convention='FULL_2'), dict(weights=np.zeros(52))]:
            with self.subTest(change=change):
                result = frontier.electronic_diagnostic({1:self.data['A2_eigenval']}, **dict(base, **change))
                self.assertEqual(result['validation_status'], 'UNSUPPORTED_OR_UNRESOLVED')
                self.assertIsNone(result['frontier_band_separation_eV'])
                self.assertTrue(result['reason'])
        for scale in (.5, 2.):
            values = self.data['A2_eigenval'].copy(); values[:, :, 1] *= scale
            self.assertIsNone(frontier.electronic_diagnostic({1:values}, **base)['frontier_band_separation_eV'])
        for change in ('nonfinite', 'order', 'count'):
            values = self.data['A2_eigenval'].copy()
            if change == 'nonfinite': values[0,0,0] = np.nan
            if change == 'order': values[0,[0,1]] = values[0,[1,0]]
            if change == 'count': values[0,0,1] = 0
            self.assertIsNone(frontier.electronic_diagnostic({1:values}, **base)['frontier_band_separation_eV'])
        too_few = self.data['A2_eigenval'][:, :31]
        self.assertIsNone(frontier.electronic_diagnostic({1:too_few}, **dict(base,nbands=31))['frontier_band_separation_eV'])

    def test_reporting_resolution_does_not_certify_tiny_gap(self):
        args = dict(nelect=2,ispin=1,nbands=2,kpoints=[[0,0,0]],weights=[1],noncollinear=False,soc=False,
                    occupation_convention=frontier.OCCUPATION_CONVENTION)
        for separation, category in [(1e-5,'NEAR_TOUCHING_UNRESOLVED'),(1e-4,'NEAR_TOUCHING_UNRESOLVED'),
                                     (1.01e-4,'FINITE_SAMPLED_FRONTIER_SEPARATION')]:
            values = {1:np.array([[[0,.5],[separation,.5]]])}
            self.assertEqual(frontier.electronic_diagnostic(values,**args)['frontier_interpretation'],category)

    def test_legacy_positive_gap_is_rejected_by_all_consumers(self):
        old = dict(sampled_gap_eV=.276, gap_classification='FINITE_SAMPLED_GAP')
        with self.assertRaisesRegex(ValueError, 'Legacy'):
            schema_record('electronic', old)
        with self.assertRaisesRegex(ValueError, 'Legacy'):
            electronic.bulk_like_edges(old, {}, {}, evidence={})
        with self.assertRaisesRegex(ValueError, 'Legacy'):
            electronic.electronic_summary({}, old, {}, spin_mode='nonmagnetic')
        with self.assertRaisesRegex(ValueError, 'insufficient'):
            electronic.electronic_summary({}, {}, {}, spin_mode='nonmagnetic')

    def test_all_40_audited_windows_and_warnings_remain(self):
        total = 0
        for name in self.expected:
            z, potential, density = self.data[name+'_profile'].T
            records = vacuum_window_sensitivity(z,potential,*self.data[name+'_extrema'],density)
            self.assertEqual([(r['setback_A'],r['side']) for r in records],
                             [(s,side) for s in (2,4,5,6) for side in ('lower','upper')])
            self.assertEqual([r['status'] for r in records], self.data[name+'_expected_window_status'].tolist())
            for i,r in enumerate(records):
                values = [r[k] for k in self.data['window_metric_columns']]
                np.testing.assert_allclose(values,self.data[name+'_expected_window_metrics'][i],atol=1e-10,rtol=0)
                self.assertEqual(r['boundary_exclusion_A'],1.)
                self.assertEqual(r['automatic_status'], 'NO_CONFIRMED_PLATEAU')
                self.assertEqual(r['discontinuity_count'],0)
            self.assertEqual(records[0]['status'],'CHARGE_IN_VACUUM_WARNING')
            self.assertEqual(records[1]['status'],'CHARGE_IN_VACUUM_WARNING')
            self.assertEqual(all(r['status']=='CANDIDATE_PLATEAU' for r in records[-2:]),name!='A0')
            total += len(records)
        self.assertEqual(total,40)

    def test_invalid_narrow_missing_density_and_jumped_windows_do_not_pass(self):
        z=np.linspace(0,30,601); potential=np.where(z<15,3.,4.); density=np.zeros_like(z)
        for extrema in ((2.,28.,30.),(0.,30.,30.),(20.,10.,30.)):
            records=vacuum_window_sensitivity(z,potential,*extrema,density)
            self.assertTrue(all(r['status'] in ('INVALID_WINDOW','INSUFFICIENT_WINDOW') for r in records))
        records=vacuum_window_sensitivity(z,potential,10,20,30)
        self.assertTrue(all(r['status']=='INSUFFICIENT_DENSITY' for r in records))
        potential[(z>=2.5)&(z<15)] += 1
        records=vacuum_window_sensitivity(z,potential,10,20,30,density)
        for r in records:
            if r['side']=='lower':
                self.assertEqual(r['status'],'AMBIGUOUS_DISCONTINUITY')
                self.assertNotIn('slope_eV_per_A',r)
                self.assertTrue(r['reason'])

    def test_force_norms_executed_target_and_pilot_roles_are_separate(self):
        for name, expected in self.expected.items():
            states=self.states(name)
            forces=self.data[name+'_forces']; maxima=np.linalg.norm(forces,axis=2).max(axis=1)
            self.assertAlmostEqual(maxima[-1],float(expected['expected_final_force_eV_per_A']),delta=1e-8)
            self.assertEqual(len(forces),int(expected['expected_ionic_steps']))
            self.assertEqual(self.data[name+'_scf_iterations'].sum(),int(expected['expected_total_scf']))
            self.assertEqual(int(self.data[name+'_ediff_markers']),len(forces))
            if name in ('A0','A1','B0'):
                self.assertEqual(states['ionic_convergence_status'],'NOT_APPLICABLE')
                self.assertIsNone(states['force_converged'])
        self.assertTrue(self.states('A2')['force_converged'])
        self.assertFalse(self.states('A2',ediffg=-.01)['force_converged'])
        self.assertFalse(self.states('A2',ionic_marker=False)['force_converged'])
        b1=self.states('B1')
        self.assertEqual(b1['ionic_convergence_status'],'NSW_LIMIT_FORCE_UNCONVERGED')
        self.assertEqual(b1['minimum_force_step'],99)
        self.assertAlmostEqual(b1['minimum_max_force_eV_per_A'],.0358072,delta=5e-8)
        # Each component passes .02, but the vector norm fails it.
        self.assertFalse(self.states('A2',forces=[[[.015,.015,0]]])['force_converged'])
        z,potential,density=self.data['A0_profile'].T
        result=dict(self.states('A0'),normal_termination=True,electronic_converged=True,geometry_identity_verified=True,
                    vacuum_window_checks=vacuum_window_sensitivity(z,potential,*self.data['A0_extrema'],density))
        summary,_=pilot.pilot_summary('A0',pilot.WP2/self.expected['A0']['source_directory'],result)
        self.assertEqual(summary['accepted_role'],'ACCEPTED_UNCORRECTED_STATIC_CONTROL')
        self.assertEqual(summary['distant_window_status'],'NO_CONFIRMED_FIELD_FREE_REFERENCE')
        self.assertEqual(summary['production_surface_gap_status'],'NOT_ESTABLISHED')
        for flag in ('timeout_evidence','fatal_evidence'):
            rejected,_=pilot.pilot_summary('A0',pilot.WP2/self.expected['A0']['source_directory'],dict(result,**{flag:True}))
            self.assertEqual(rejected['accepted_role'],'UNRESOLVED')

    def test_import_and_help_cannot_write_or_launch_external_programs(self):
        guard = '''
import os, sys, runpy
def guard(event, args):
    if event == 'open':
        mode, flags = args[1:3]
        if (isinstance(mode,str) and any(c in mode for c in 'wax+')) or (flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC)):
            raise AssertionError('Unexpected write: '+str(args[0]))
    if event in ('os.mkdir','os.remove','os.rename','os.rmdir','subprocess.Popen','os.system','os.posix_spawn','socket.connect'):
        raise AssertionError('Unexpected side effect: '+event)
sys.addaudithook(guard)
'''
        for action in ("import analyse_runtime_pilots, pilot_frontier, electronic_analysis, electrostatics, result_schemas",
                       "sys.argv=['analyse_runtime_pilots','--help']; runpy.run_module('analyse_runtime_pilots',run_name='__main__')"):
            completed=subprocess.run([sys.executable,'-c',guard+'\n'+action],capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stderr)

    @unittest.skipUnless(os.environ.get('WP2_INTEGRITY_MANIFEST'),'Explicit local protected-evidence inventory required')
    def test_protected_native_inputs_and_historical_evidence_are_unchanged(self):
        with Path(os.environ['WP2_INTEGRITY_MANIFEST']).open() as stream:
            records=list(csv.DictReader(stream))
        self.assertGreater(len(records),75)
        for r in records:
            source=Path(r['path']); h=hashlib.sha256()
            with source.open('rb') as stream:
                for block in iter(lambda:stream.read(1024*1024),b''): h.update(block)
            self.assertEqual(source.stat().st_size,int(r['bytes']),str(source))
            self.assertEqual(h.hexdigest(),r['sha256'],str(source))


if __name__ == '__main__':
    unittest.main()
