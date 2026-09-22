"""Pure scientific-template checks; structures stay in memory and no job is made."""

import json
import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

from stage_templates import (approved_context, bulk_static_template, reconstruction_probe,
                             relaxation_template, spin_probe_templates,
                             surface_static_template, surface_template)


class StageTemplateTests(unittest.TestCase):
    def setUp(self):
        self.approval = dict(approved=True, surface_id="fixture_surface",
                             face_ids=dict(upper="fixture_upper", lower="fixture_lower"),
                             review_reference="synthetic fixture approval", parent_phase="beta",
                             shortlisted=True, final_relaxation_review="synthetic final review")
        self.numerical = dict(encut_eV=500, k_mesh=[8, 8, 1],
                              selected_by_review="synthetic selected settings",
                              integration_compatibility_id="synthetic integration family",
                              integration_review_reference="synthetic integration review")
        self.paws = dict(Zn="fixture:Zn-identity", In="fixture:In_d-identity", S="fixture:S-identity")

    def test_explicit_approval_and_slab_sampling(self):
        with self.assertRaises(ValueError):
            approved_context(dict(self.approval, approved=False))
        with self.assertRaises(ValueError):
            approved_context(dict(self.approval, face_ids={"upper": "new-face"}))
        with self.assertRaises(ValueError):
            surface_template(self.approval, dict(self.numerical, k_mesh=[8, 8, 2]))
        spec = surface_template(self.approval, self.numerical)
        self.assertEqual(spec["context"]["face_ids"], self.approval["face_ids"])
        self.assertEqual(spec["incar"]["ICHARG"], 2)
        self.assertEqual(spec["incar"]["IVDW"], 12)
        self.assertFalse(spec["incar"]["LWAVE"])
        self.assertEqual(spec["incar"]["ISMEAR"], 0)
        self.assertEqual(spec["incar"]["SIGMA"], .05)
        self.assertEqual(spec["incar"]["ALGO"], "Normal")
        self.assertEqual(spec["incar"]["AMIN"], .01)
        self.assertEqual(spec["incar"]["NELM"], 120)
        self.assertNotIn("AMIX", spec["incar"])
        self.assertNotIn("BMIX", spec["incar"])

    def test_relaxation_tiers_and_one_question_gate(self):
        screening = relaxation_template(self.approval, self.numerical)
        final = relaxation_template(self.approval, self.numerical, tier="final_relaxation")
        for spec in (screening, final):
            self.assertEqual(spec["incar"]["ISIF"], 2)
            self.assertEqual(spec["incar"]["ISYM"], 0)
            self.assertEqual(spec["incar"]["AMIN"], .01)
            self.assertEqual(spec["incar"]["NELM"], 120)
        self.assertEqual(screening["incar"]["EDIFF"], 1e-6)
        self.assertEqual(screening["incar"]["EDIFFG"], -.02)
        self.assertEqual(final["incar"]["EDIFF"], 1e-7)
        self.assertEqual(final["incar"]["EDIFFG"], -.01)
        with self.assertRaises(ValueError):
            relaxation_template(dict(self.approval, shortlisted=False), self.numerical, tier="final_relaxation")
        with self.assertRaises(ValueError):
            relaxation_template(self.approval, self.numerical, tier="shortlist_validation")
        with self.assertRaises(ValueError):
            relaxation_template(self.approval, dict(self.numerical, selected_by_review=None))
        valid = relaxation_template(self.approval, self.numerical, tier="shortlist_validation",
                                     validation_check="relaxed_region", reason="synthetic sensitivity question")
        self.assertEqual(valid["validation_check"], "relaxed_region")

    def test_explicit_spin_initializations_and_smearing(self):
        probes = spin_probe_templates(self.approval, self.numerical, atom_count=3,
                                      initializations={"supplied": [0, 1, -1]}, reason="synthetic spin question")
        self.assertEqual(len(probes), 1)
        self.assertEqual(probes[0]["incar"]["MAGMOM"], [0., 1., -1.])
        self.assertEqual(probes[0]["incar"]["ISPIN"], 2)
        with self.assertRaises(ValueError):
            spin_probe_templates(self.approval, self.numerical, atom_count=3,
                                 initializations={"bad": [1, -1]}, reason="synthetic")
        with self.assertRaises(ValueError):
            surface_template(self.approval, self.numerical, smearing=dict(ismear=1, sigma_eV=.1))
        changed = surface_template(self.approval, self.numerical,
                                   smearing=dict(ismear=1, sigma_eV=.1, reason="synthetic metallic-state sensitivity"))
        self.assertEqual(changed["incar"]["ISMEAR"], 1)

    def test_final_static_outputs_and_wavecar_gate(self):
        spec = surface_static_template(self.approval, self.numerical, paw_identities=self.paws)
        self.assertEqual(spec["incar"]["NSW"], 0)
        self.assertNotEqual(spec["incar"]["ICHARG"], 11)
        self.assertTrue(spec["incar"]["LCHARG"])
        self.assertTrue(spec["incar"]["LVHAR"])
        self.assertEqual(spec["incar"]["LORBIT"], 11)
        self.assertEqual(spec["incar"]["AMIN"], .01)
        self.assertEqual(spec["incar"]["NELM"], 120)
        self.assertFalse(spec["incar"]["LWAVE"])
        wave = surface_static_template(self.approval, self.numerical, paw_identities=self.paws,
                                       band_edge_charge_reason="synthetic explicitly selected states")
        self.assertTrue(wave["incar"]["LWAVE"])
        with self.assertRaises(ValueError):
            surface_static_template(dict(self.approval, shortlisted=False), self.numerical,
                                    paw_identities=self.paws, band_edge_charge_reason="synthetic")

    def test_accepted_bulk_geometry_and_matching(self):
        directory = Path(__file__).resolve().parents[2] / "structure/bulk_parents/beta"
        data = (directory / "POSCAR").read_bytes()
        metadata = json.loads((directory / "metadata.json").read_text())
        surface = surface_static_template(self.approval, self.numerical, paw_identities=self.paws)
        numerical = dict(self.numerical, k_mesh=[8, 8, 3])
        spec = bulk_static_template("beta", data, metadata, numerical=numerical,
                                    paw_identities=self.paws, matching_surface=surface,
                                    integration_review_reference="synthetic reciprocal-space review")
        self.assertEqual(len(spec["structure"]), 7)
        self.assertEqual(spec["incar"]["IBRION"], -1)
        self.assertEqual(spec["vasp_version"], "6.6.1")
        self.assertFalse(spec["incar"]["LDIPOL"])
        self.assertEqual(spec["incar"]["AMIN"], .01)
        self.assertEqual(spec["incar"]["NELM"], 120)
        candidate_numerical = dict(numerical, encut_eV=600, selected_by_review=None)
        candidate_surface = surface_template(self.approval, dict(candidate_numerical, k_mesh=[8, 8, 1]),
                                             paw_identities=self.paws)
        candidate = bulk_static_template("beta", data, metadata, numerical=candidate_numerical,
                                         paw_identities=self.paws, matching_surface=candidate_surface,
                                         integration_review_reference="synthetic compatibility review",
                                         purpose="convergence_reference", convergence_case_review="one 600 eV check")
        self.assertEqual(candidate["encut_eV"], 600)
        self.assertEqual(candidate["purpose"], "convergence_reference")
        for phase, payload, numerical_case in (("IIb", data, numerical),
                                                ("beta", data + b"\n", numerical),
                                                ("beta", data, dict(numerical, encut_eV=600))):
            with self.assertRaises(ValueError):
                bulk_static_template(phase, payload, metadata, numerical=numerical_case,
                                     paw_identities=self.paws, matching_surface=surface,
                                     integration_review_reference="synthetic review")

    def test_reconstruction_mapping_determinism_and_normal_geometry(self):
        slab = Structure(Lattice.orthorhombic(4, 4, 40), ["Zn", "In", "In", "S", "S", "S", "S"],
                         [[1, 1, z] for z in (10, 13, 15, 20, 25, 27, 30)], coords_are_cartesian=True)
        options = dict(supercell=(2, 1), surface_indices=[0, 6], amplitude_A=.03,
                       reason="synthetic reconstruction question", seed=17)
        first, record = reconstruction_probe(slab, self.approval, **options)
        second, repeated = reconstruction_probe(slab, self.approval, **options)
        self.assertTrue(np.array_equal(first.cart_coords, second.cart_coords))
        self.assertEqual(record, repeated)
        self.assertEqual(len(first), 14)
        self.assertEqual(first.composition, slab.composition * 2)
        self.assertTrue(np.array_equal(first.lattice.matrix[2], slab.lattice.matrix[2]))
        for i, source in enumerate(record["source_site_indices"]):
            self.assertAlmostEqual(first.cart_coords[i, 2], slab.cart_coords[source, 2])
        self.assertEqual(record["context"]["face_ids"], self.approval["face_ids"])
        self.assertEqual(set(record["source_site_indices"][i] for i in record["perturbed_site_indices"]), {0, 6})
        with self.assertRaises(ValueError):
            reconstruction_probe(slab, self.approval, **dict(options, surface_indices=[3]))
        with self.assertRaises(ValueError):
            reconstruction_probe(slab, self.approval, **dict(options, supercell=(3, 3)))


if __name__ == "__main__":
    unittest.main()
