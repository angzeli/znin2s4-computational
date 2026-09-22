"""Protect root resolution and export destinations after flattening script/."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import analyse_slab_interior as interior
import generate_surface_candidates as generate
import prepare_runtime_pilots as pilots


class ScriptLayoutTests(unittest.TestCase):
    def test_sources_resolve_within_wp2(self):
        wp2 = Path(__file__).resolve().parents[2]
        for module in (interior, generate, pilots):
            self.assertEqual(module.WP2, wp2)
        self.assertTrue(interior.CALC.is_dir())
        self.assertTrue(generate.AUDIT.is_file())
        self.assertTrue(pilots.POTCAR_DONOR.is_file())

    def test_exports_target_results_not_scripts(self):
        self.assertEqual(interior.RESULTS,
                         interior.WP2 / 'results/02_numerical_convergence')

    def test_csv_export_does_not_recreate_removed_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'profiles'
            with patch.dict('sys.modules', {'matplotlib': None}):
                interior.export({}, output_directory=output)
            self.assertEqual([p.name for p in output.iterdir()],
                             ['03_SLAB_INTERIOR_PROFILES.csv'])
            self.assertTrue((output / '03_SLAB_INTERIOR_PROFILES.csv')
                            .read_text().startswith('layers,z_A,'))


if __name__ == '__main__':
    unittest.main()
