"""Full frozen pre-refactor Stage 01 outputs and pure grid-adapter regression."""

import json
import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Chgcar, Locpot

import analyse_runtime_pilots as pilot
import electrostatics as shared


class ElectrostaticsRegressionTests(unittest.TestCase):
    def test_frozen_stage01_outputs_and_compatibility_exports(self):
        frozen = json.loads((Path(__file__).parents[1]/'fixtures/stage01_electrostatics.json').read_text())
        z = np.arange(0., 30.01, .5); potential = np.where(z < 15, 3., 4.)
        inputs = [('two_faces', potential, None), ('sloped', .08*z, None),
                  ('charged', potential, np.full_like(z, .001)),
                  ('ambiguous', potential+((z >= 4.5)&(z < 15)), None)]
        self.assertIs(pilot.vacuum_windows, shared.vacuum_windows)
        self.assertIs(pilot.fit_window, shared.fit_window)
        self.assertIs(pilot.same_geometry, shared.same_geometry)
        for name, values, density in inputs:
            self.assertEqual(shared.vacuum_windows(z, values, 10, 20, 30, density), frozen[name])

    def test_planar_profile_density_units_and_bad_grid(self):
        slab = Structure(Lattice.from_parameters(4, 4, 30, 90, 90, 90), ['S','S'], [[0,0,1/3], [0,0,2/3]])
        potential = Locpot(slab, {'total': np.zeros((2,2,120))})
        charge = Chgcar(slab, {'total': np.full((2,2,120), 12.)})
        result = shared.planar_profile(potential, charge=charge, nelect=12.)
        np.testing.assert_allclose(result['density_e_per_A3'], 12/slab.volume)
        self.assertEqual(result['charge_density_status'], 'VALIDATED_NELECT_INTEGRAL')
        charge = Chgcar(slab, {'total': np.full((2,2,60), 12.)})
        result = shared.planar_profile(potential, charge=charge, nelect=12.)
        self.assertIsNone(result['density_e_per_A3'])
        self.assertIn('no interpolation', result['charge_density_status'])
        # Equal z lengths do not establish identical three-dimensional grids.
        charge = Chgcar(slab, {'total': np.full((4,2,120), 12.)})
        result = shared.planar_profile(potential, charge=charge, nelect=12.)
        self.assertIsNone(result['density_e_per_A3'])
        self.assertIn('grids differ', result['charge_density_status'])


if __name__ == '__main__':
    unittest.main()
