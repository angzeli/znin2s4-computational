"""Synthetic result-identity joins; no native outputs or scientific engine."""
import unittest

from convergence_analysis import compare_cases, convergence_record


class ConvergenceHandoffTests(unittest.TestCase):
    def records(self, case_id="case", energy=-38.):
        faces = {"surface_id": "fixture", "face_ids": {"upper": "fixture_upper", "lower": "fixture_lower"}}
        prepared = dict(case_id=case_id, surface_id="fixture", face_context=faces,
                        dimension="encut", baseline_case_id="baseline", target_value=600,
                        actual_value=dict(encut_eV=600, k_mesh=[2, 2, 1], vacuum_A=20,
                                          slab_thickness_A=10, relaxed_region_mode="full"),
                        comparison_protocol={"incar": {"NSW": 0}, "energy_convention": "energy(sigma->0)"})
        completed = dict(case_id=case_id, static_result_id=case_id + "-static",
                         input_settings_verified=True, normal_termination=True,
                         electronic_converged=True, energy_eV=energy)
        energetics = dict(surface_id="fixture", face_context=faces,
                          slab_static_result_id=completed["static_result_id"], E_slab_eV=energy,
                          bulk_reference_id="fixture-bulk", Gamma_pair_meV_A2=(energy + 40.) / 20. * 1000.,
                          matched_reference_validated=True)
        return prepared, completed, energetics

    def test_matching_result_and_energy_preserve_numeric_comparison(self):
        baseline = convergence_record(*self.records("baseline", -39.))
        baseline["encut_eV"] = 500
        case = convergence_record(*self.records())
        self.assertEqual(case["energy_eV"], -38.)
        result = compare_cases(baseline, case, electronic_relevant=False)
        self.assertEqual(result["numerical_status"], "NOT_CONVERGED")
        self.assertEqual(result["human_review_status"], "REQUIRED")

    def test_foreign_static_cannot_join_on_shared_surface_and_faces(self):
        prepared, completed, _ = self.records()
        _, _, foreign = self.records("baseline", -39.)
        with self.assertRaisesRegex(ValueError, "static result"):
            convergence_record(prepared, completed, foreign)
        # Equal energy does not make a foreign source identity interchangeable.
        foreign["E_slab_eV"] = completed["energy_eV"]
        with self.assertRaisesRegex(ValueError, "static result"):
            convergence_record(prepared, completed, foreign)

    def test_missing_empty_and_nonstring_result_identities_refuse(self):
        for value in (None, "", " ", 1, True):
            for owner in ("completed", "energetics"):
                with self.subTest(value=value, owner=owner):
                    prepared, completed, energetics = self.records()
                    if owner == "completed":
                        completed["static_result_id"] = value
                    else:
                        energetics["slab_static_result_id"] = value
                    with self.assertRaisesRegex(ValueError, "static result"):
                        convergence_record(prepared, completed, energetics)
        for owner, key in ((1, "static_result_id"), (2, "slab_static_result_id")):
            records = self.records()
            del records[owner][key]
            with self.assertRaisesRegex(ValueError, "static result"):
                convergence_record(*records)

    def test_changed_energy_refuses_without_new_scientific_tolerance(self):
        prepared, completed, energetics = self.records()
        for energy in (-39., completed["energy_eV"] + 1e-10):
            with self.subTest(energy=energy):
                energetics["E_slab_eV"] = energy
                with self.assertRaisesRegex(ValueError, "energy"):
                    convergence_record(prepared, completed, energetics)

    def test_missing_nonfinite_and_non_numeric_energies_refuse(self):
        for value in (None, float("nan"), float("inf"), True, "-38"):
            for owner in ("completed", "energetics"):
                with self.subTest(value=value, owner=owner):
                    prepared, completed, energetics = self.records()
                    if owner == "completed":
                        completed["energy_eV"] = value
                    else:
                        energetics["E_slab_eV"] = value
                    with self.assertRaisesRegex(ValueError, "energy"):
                        convergence_record(prepared, completed, energetics)
        for owner, key in ((1, "energy_eV"), (2, "E_slab_eV")):
            records = self.records()
            del records[owner][key]
            with self.assertRaisesRegex(ValueError, "energy"):
                convergence_record(*records)

    def test_existing_completed_identity_settings_and_face_checks_remain(self):
        for change in ({"case_id": "wrong"}, {"input_settings_verified": False}):
            prepared, completed, energetics = self.records()
            completed.update(change)
            with self.assertRaises(ValueError):
                convergence_record(prepared, completed, energetics)
        prepared, completed, energetics = self.records()
        energetics["face_context"] = {"surface_id": "other"}
        with self.assertRaises(ValueError):
            convergence_record(prepared, completed, energetics)


if __name__ == "__main__":
    unittest.main()
