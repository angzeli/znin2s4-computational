"""Synthetic energies validate units, face conventions and matched-reference failures."""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Kpoints, Poscar

from surface_energetics import read_static_record, surface_energy, surface_property_rows, validate_matched_reference
from stage_templates import surface_template


class SurfaceEnergeticsTests(unittest.TestCase):
    def setUp(self):
        method = dict(parent_phase="beta", xc="PBE", dispersion="D3(BJ)",
                      paw_identities=dict(Zn="fixture:Zn", In="fixture:In_d", S="fixture:S"),
                      encut_eV=500, energy_convention="energy(sigma->0)",
                      vasp_version="6.6.1", formula_unit="ZnIn2S4",
                      integration_compatibility_id="fixture reviewed integration family",
                      integration_review_reference="fixture review")
        self.slab = dict(method, surface_id="fixture_surface",
                         face_ids=dict(upper="fixture_upper", lower="fixture_lower"),
                         composition=dict(Zn=2, In=4, S=8), n_formula_units=2,
                         E_slab_eV=-39., surface_area_A2=20.)
        self.bulk = dict(method, bulk_reference_id="fixture_bulk_static",
                         composition=dict(Zn=2, In=4, S=8), n_formula_units=2, E_bulk_eV=-40.)

    def test_pair_formula_normalization_and_units(self):
        result = surface_energy(self.slab, self.bulk, faces_equivalent=False)
        self.assertAlmostEqual(result["e_bulk_eV_fu"], -20.)
        self.assertAlmostEqual(result["Gamma_pair_eV_A2"], .05)
        self.assertAlmostEqual(result["Gamma_pair_meV_A2"], 50.)
        self.assertIsNone(result["gamma_single_meV_A2"])
        self.assertEqual(result["face_context"]["face_ids"], self.slab["face_ids"])
        self.assertTrue(result["matched_reference_validated"])
        self.assertEqual(result["status"], "COMPUTED_REQUIRES_HUMAN_REVIEW")

    def test_equivalent_faces_need_evidence_and_inequivalent_remain_unassigned(self):
        with self.assertRaises(ValueError):
            surface_energy(self.slab, self.bulk, faces_equivalent=True)
        equal = surface_energy(self.slab, self.bulk, faces_equivalent=True,
                               face_equivalence_evidence="synthetic symmetry exchange")
        self.assertAlmostEqual(equal["gamma_single_meV_A2"], 25.)
        unknown = surface_energy(self.slab, self.bulk)
        self.assertIsNone(unknown["faces_equivalent"])
        self.assertIsNone(unknown["gamma_single_eV_A2"])

    def test_parent_method_dataset_and_energy_convention_mismatches_reject(self):
        differences = dict(parent_phase="spinel", xc="r2SCAN", dispersion="none", encut_eV=600,
                           energy_convention="free_energy_TOTEN", vasp_version="5.4.4",
                           formula_unit="unknown", integration_compatibility_id="other integration")
        differences["paw_identities"] = dict(Zn="other", In="fixture:In_d", S="fixture:S")
        for key, value in differences.items():
            with self.subTest(field=key), self.assertRaises(ValueError):
                surface_energy(self.slab, dict(self.bulk, **{key: value}))

    def test_undefined_normalization_and_missing_matching_evidence_reject(self):
        validate_matched_reference(self.slab, self.bulk)
        for change in (dict(n_formula_units=1), dict(energy_convention=None),
                       dict(integration_review_reference=""), dict(encut_eV=float("nan"))):
            with self.subTest(change=change), self.assertRaises(ValueError):
                surface_energy(self.slab, dict(self.bulk, **change))

    def test_nonstoichiometric_slab_gets_no_naive_surface_energy(self):
        result = surface_energy(dict(self.slab, composition=dict(Zn=2, In=4, S=7)), self.bulk)
        self.assertEqual(result["status"], "UNSUPPORTED_NONSTOICHIOMETRIC")
        self.assertIsNone(result["Gamma_pair_meV_A2"])
        self.assertIsNone(result["gamma_single_meV_A2"])
        self.assertIsNone(result["n_formula_units"])

    def test_total_energies_and_one_face_area_must_be_finite(self):
        for changes in (dict(surface_area_A2=0), dict(surface_area_A2=-1), dict(E_slab_eV=float("nan"))):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                surface_energy(dict(self.slab, **changes), self.bulk)

    def test_surface_properties_adapter_preserves_two_face_contexts(self):
        result = surface_energy(self.slab, self.bulk, faces_equivalent=False)
        rows = surface_property_rows(result, properties=dict(atom_count=14))
        self.assertEqual({row["face_id"] for row in rows}, {"fixture_upper", "fixture_lower"})
        self.assertTrue(all(row["Gamma_pair_meV_A2"] == 50 for row in rows))
        self.assertTrue(all(row["gamma_single_meV_A2"] is None for row in rows))
        self.assertTrue(all(row["vacuum_upper_eV"] is None for row in rows))

    def test_static_reader_requires_completed_observed_method_and_mesh(self):
        approval = dict(approved=True, surface_id="fixture_surface", parent_phase="beta",
                        face_ids=self.slab["face_ids"], review_reference="fixture approval")
        numerical = dict(encut_eV=500, k_mesh=[2, 2, 1],
                         integration_compatibility_id="fixture integration",
                         integration_review_reference="fixture review")
        template = surface_template(approval, numerical, paw_identities=self.slab["paw_identities"])
        structure = Structure(Lattice.orthorhombic(4, 4, 40), ["Zn", "In", "In", "S", "S", "S", "S"],
                              [[1, 1, z] for z in (10, 13, 15, 20, 25, 27, 30)], coords_are_cartesian=True)
        mesh = Kpoints.gamma_automatic((2, 2, 1))
        run = SimpleNamespace(vasp_version="6.6.1", converged_electronic=True,
                              ionic_steps=[dict(electronic_steps=[{}, {}], e_0_energy=-999, e_fr_energy=-19.51)],
                              final_energy=-19.5,
                              incar=template["incar"], parameters=template["incar"],
                              initial_structure=structure, final_structure=structure, kpoints=mesh)
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "OUTCAR").write_text("aborting loop because EDIFF is reached\nGeneral timing and accounting informations for this job:\n")
            (directory / "POSCAR").write_text(str(Poscar(structure)))
            (directory / "KPOINTS").write_text(str(mesh))
            with patch("surface_energetics.Vasprun", return_value=run), \
                 patch("surface_energetics._paw_identities_from_file", return_value=self.slab["paw_identities"]):
                record = read_static_record(directory, template, kind="surface", result_id="fixture-static", expected_structure=structure)
                self.assertEqual(record["E_slab_eV"], -19.5)
                self.assertEqual(record["case_id"], "fixture-static")
                self.assertTrue(record["input_settings_verified"])
                self.assertEqual(record["n_formula_units"], 1)
                run.kpoints = Kpoints.gamma_automatic((3, 3, 1))
                with self.assertRaises(ValueError):
                    read_static_record(directory, template, kind="surface", result_id="fixture-static", expected_structure=structure)
                run.kpoints = mesh
                run.parameters = dict(template["incar"], ENCUT=600)
                with self.assertRaises(ValueError):
                    read_static_record(directory, template, kind="surface", result_id="fixture-static", expected_structure=structure)
                run.parameters = template["incar"]
                for override in (dict(NELECT=100), dict(NUPDOWN=2), dict(EFIELD=.1)):
                    run.incar = dict(template["incar"], **override)
                    with self.subTest(override=override), self.assertRaises(ValueError):
                        read_static_record(directory, template, kind="surface", result_id="fixture-static", expected_structure=structure)
                run.incar = template["incar"]
                (directory / "OUTCAR").write_text("incomplete synthetic output")
                with self.assertRaises(ValueError):
                    read_static_record(directory, template, kind="surface", result_id="fixture-static", expected_structure=structure)


if __name__ == "__main__":
    unittest.main()
