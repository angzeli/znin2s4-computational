"""Focused input, human-gate and numerical-diagnostic tests; no engine is run."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Incar

import prepare_runtime_pilots as prepare
import analyse_runtime_pilots as analyse


class RuntimePilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.files,cls.specs=prepare.build_plan()

    def test_exactly_five_jobs_and_atomic_hold_plan(self):
        self.assertEqual([s['pilot_job'] for s in self.specs],['A0','A1','B0','A2','B1'])
        for spec in self.specs:
            is_relax=spec['pilot_job'] in ('A2','B1')
            self.assertEqual('--hold' in spec['cmw_add_arguments'],is_relax)
            self.assertEqual(spec['requires_explicit_human_release'],is_relax)
            self.assertEqual(spec['cmw_state'],str(prepare.STATE))
            self.assertEqual(spec['timeout_seconds'],172800 if is_relax else 21600)
            self.assertEqual(spec['launcher_argv'][spec['launcher_argv'].index('--timeout')+1],
                             str(spec['timeout_seconds']))
            self.assertNotIn('--allow-auxiliary',spec['cmw_add_arguments'])
            self.assertIn('--managed-foreground',spec['launcher_argv'])
            for flag,value in [('--ranks','8'),('--ncore','4'),('--kpar','1'),('--stop-before','0')]:
                self.assertEqual(spec['launcher_argv'][spec['launcher_argv'].index(flag)+1],value)
            for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','CMW_JOBS_OWN_SESSION'):
                self.assertEqual(spec['environment'][key],'1')
        self.assertFalse(any(p.name=='POTCAR' for p in self.files))

    def test_control_pair_and_sampling(self):
        a0,a1,b0,a2,b1=self.specs
        differences={k for k in set(a0['incar'])|set(a1['incar']) if a0['incar'].get(k)!=a1['incar'].get(k)}
        self.assertEqual(differences,{'SYSTEM','LDIPOL','IDIPOL','DIPOL'})
        self.assertEqual(a0['source_structure_hash'],a1['source_structure_hash'])
        self.assertEqual(a0['kpoints']['mesh'],[10,10,1])
        self.assertEqual(b0['kpoints']['mesh'],[9,3,1])
        self.assertEqual(a2['kpoints'],a0['kpoints'])
        self.assertEqual(b1['kpoints'],b0['kpoints'])
        for spec in self.specs:
            self.assertLessEqual(spec['kpoints']['maximum_inplane_spacing_inv_A'],.20)
        self.assertEqual(a0['potcar_dataset_order'],['S','In_d','Zn'])
        self.assertEqual(b0['potcar_dataset_order'],['S','In_d','Zn']*2)

    def test_fail_closed_on_input_changes(self):
        spec=self.specs[0]
        directory=Path(spec['input_directory']).relative_to(prepare.WP2.relative_to(prepare.REPO))
        inputs={name:self.files[directory/name] for name in ('POSCAR','INCAR','KPOINTS')}
        for name,content in [('POSCAR',inputs['POSCAR']+b'\n'),('INCAR',inputs['INCAR']+b'LHFCALC = .TRUE.\n'),
                             ('KPOINTS',inputs['KPOINTS'].replace(b'10 10 1',b'10 10 4'))]:
            changed=dict(inputs);changed[name]=content
            with self.assertRaises(ValueError):
                prepare.validate_inputs(spec,changed)
        changed=json.loads(json.dumps(self.specs[3]))
        changed['cmw_add_arguments'].remove('--hold')
        directory=Path(changed['input_directory']).relative_to(prepare.WP2.relative_to(prepare.REPO))
        with self.assertRaises(ValueError):
            prepare.validate_inputs(changed,{name:self.files[directory/name] for name in inputs})

    def test_determinism_and_no_implicit_writes(self):
        with patch.object(Path,'write_text',side_effect=AssertionError('implicit write')), \
             patch.object(Path,'write_bytes',side_effect=AssertionError('implicit write')), \
             patch.object(Path,'mkdir',side_effect=AssertionError('implicit directory')):
            files,specs=prepare.build_plan()
        self.assertEqual(files,self.files)
        self.assertEqual(specs,self.specs)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'existing'
            path.write_bytes(b'protected')
            with self.assertRaises(FileExistsError):
                prepare.write_unchanged_or_new(path,b'different')
            self.assertEqual(path.read_bytes(),b'protected')

    def test_flat_sloped_charged_and_ambiguous_vacuum(self):
        z=np.linspace(0,30,601)
        flat=np.where(z<15,3.,4.)
        density=np.zeros_like(z)
        result=analyse.vacuum_windows(z,flat,10,20,30,density)
        self.assertAlmostEqual(result['upper_minus_lower_potential_eV'],1.)
        self.assertAlmostEqual(result['lower']['selected_window']['slope_eV_per_A'],0.)
        sloped=analyse.vacuum_windows(z,.08*z,10,20,30)
        self.assertIsNone(sloped['upper_minus_lower_potential_eV'])
        self.assertAlmostEqual(sloped['lower']['geometric_window']['slope_eV_per_A'],.08)
        charged=analyse.vacuum_windows(z,flat,10,20,30,np.full_like(z,1e-3))
        self.assertEqual(charged['lower']['geometric_window']['status'],'CHARGE_IN_VACUUM_WARNING')
        jump=flat.copy();jump[(z>=4.5)&(z<15)]+=1
        ambiguous=analyse.vacuum_windows(z,jump,10,20,30)
        self.assertEqual(ambiguous['lower']['status'],'AMBIGUOUS_MULTIPLE_WINDOWS')
        self.assertIsNone(ambiguous['lower']['selected_window'])

    def test_sampled_occupation_diagnostic_is_bounded(self):
        # The retired positive gap discarded partial frontier states. The new
        # quantity keeps them and separates interpretation from occupations.
        finite={1:np.array([[[-1,1],[1,0]]])}
        args=dict(nelect=2,ispin=1,nbands=2,kpoints=[[0,0,0]],weights=[1],
                  noncollinear=False,soc=False,occupation_convention=analyse.OCCUPATION_CONVENTION)
        self.assertEqual(analyse.electronic_diagnostic(finite,**args)['frontier_interpretation'],'FINITE_SAMPLED_FRONTIER_SEPARATION')
        partial={1:np.array([[[-1,.9],[1,.1]]])}
        self.assertTrue(analyse.electronic_diagnostic(partial,**args)['partial_occupation_present'])
        self.assertEqual(analyse.electronic_diagnostic({})['validation_status'],'UNSUPPORTED_OR_UNRESOLVED')

    def test_completion_requires_electronic_and_termination_evidence(self):
        structure=Structure(Lattice.cubic(20),['S'],[[.5,.5,.5]])
        step=dict(electronic_steps=[{},{}],forces=[[0,0,0]],structure=structure)
        run=SimpleNamespace(ionic_steps=[step],converged_electronic=True,converged_ionic=True,
                            initial_structure=structure,final_structure=structure,final_energy=-1.,
                            parameters={'NELECT':6,'NELM':150,'NSW':0,'IBRION':-1},
                            incar={'NSW':0,'NELM':150},eigenvalues={},vasp_version='6.6.1')
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            (directory/'INCAR').write_text(str(Incar(dict(NSW=0,NELM=150))))
            from pymatgen.io.vasp import Poscar
            for name in ('POSCAR','CONTCAR'):
                (directory/name).write_text(str(Poscar(structure)))
            out='NELECT = 6\nIteration 1(1)\naborting loop because EDIFF is reached\nGeneral timing and accounting informations for this job:\nElapsed time (sec): 12.0\n'
            (directory/'OUTCAR').write_text(out)
            with patch.object(analyse,'Vasprun',return_value=run):
                self.assertEqual(analyse.output_evidence(directory)[0]['scientific_status'],'PASS')
                (directory/'OUTCAR').write_text(out.replace('aborting loop because EDIFF is reached',''))
                self.assertEqual(analyse.output_evidence(directory)[0]['scientific_status'],'INCONCLUSIVE')
                (directory/'OUTCAR').write_text(out)
                (directory/'RUN_METADATA.txt').write_text('reason: timeout\n')
                self.assertEqual(analyse.output_evidence(directory)[0]['scientific_status'],'INCONCLUSIVE')
        absent=Path('/private/tmp/wp2-nonexistent-test-runtime')
        with patch.object(Path,'exists',return_value=False):
            self.assertEqual(analyse.output_evidence(absent)[0]['scientific_status'],'PREPARED')

    def test_density_normalization_and_profile_geometry(self):
        structure=Structure(Lattice.from_parameters(4,4,30,90,90,90),['S','S'],[[0,0,1/3],[0,0,2/3]])
        charge_grid=np.full((2,2,120),12.)
        potential_grid=np.zeros((2,2,120))
        charge=analyse.Chgcar(structure,{'total':charge_grid})
        locpot=analyse.Locpot(structure,{'total':potential_grid})
        run=SimpleNamespace(final_structure=structure)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'INCAR').write_text('LVHAR = .TRUE.\n')
            with patch.object(analyse,'output_evidence',return_value=({'nelect':12.,'normal_termination':True,
                 'electronic_converged':True,'geometry_identity_verified':True},run)), \
                 patch.object(analyse.Locpot,'from_file',return_value=locpot), \
                 patch.object(analyse.Chgcar,'from_file',return_value=charge):
                result,profile=analyse.analyse_job(root)
            self.assertEqual(result['charge_density_status'],'VALIDATED_NELECT_INTEGRAL')
            self.assertAlmostEqual(profile[0]['mean_electron_density_e_per_A3'],12/structure.volume)
            self.assertEqual(result['vacuum']['lower']['geometric_window']['status'],'CHARGE_IN_VACUUM_WARNING')


if __name__=='__main__':
    unittest.main()
