"""Temporary/in-memory geometry checks using frozen Stage 0 copies only."""

import json
import unittest
from pathlib import Path

import numpy as np
from pymatgen.io.vasp import Poscar

from surface_builders import change_vacuum, extend_thickness, slab_geometry
from surface_layers import group_surface_layers, propose_relaxed_region

WP2 = Path(__file__).resolve().parents[2]


def seed_fixture(phase="beta", facet="001", termination="t07"):
    """Read accepted copies under WP2, never current WP1 calculations."""
    root = WP2 / "structure/slab_candidates" / phase / facet / f"{phase}_{facet}_{termination}"
    slab = Poscar.from_file(root / "POSCAR", check_for_potcar=False).structure
    parent = Poscar.from_file(WP2 / "structure/bulk_parents" / phase / "POSCAR",
                             check_for_potcar=False).structure
    return parent, slab, json.loads((root / "metadata.json").read_text())


class SurfaceBuilderTests(unittest.TestCase):
    def test_vacuum_cartesian_invariance_and_roundtrip(self):
        _, seed, _ = seed_fixture()
        pairs = seed.cart_coords[:, None, :] - seed.cart_coords[None, :, :]
        for target in (20, 25, 30):
            result, metadata = change_vacuum(seed, target)
            np.testing.assert_array_equal(result.lattice.matrix[:2], seed.lattice.matrix[:2])
            np.testing.assert_allclose(result.cart_coords[:, None, :] - result.cart_coords[None, :, :],
                                       pairs, atol=1e-14, rtol=0)
            self.assertAlmostEqual(metadata["vacuum_A"], target)
            self.assertAlmostEqual(metadata["lower_A"], target / 2)
            returned, _ = change_vacuum(result, 20)
            np.testing.assert_allclose(returned.cart_coords, seed.cart_coords, atol=1e-14, rtol=0)
            self.assertEqual(result.species, seed.species)

    def test_vacuum_does_not_mutate_site_properties_or_seed(self):
        _, seed, _ = seed_fixture()
        seed.add_site_property("source_index", list(range(len(seed))))
        before = seed.cart_coords.copy()
        result, _ = change_vacuum(seed, 25)
        self.assertEqual(result.site_properties, seed.site_properties)
        np.testing.assert_array_equal(seed.cart_coords, before)
        with self.assertRaises(ValueError):
            change_vacuum(seed, 0)

    def test_periodically_wrapped_envelope_fails_closed(self):
        _, seed, _ = seed_fixture()
        wrapped = seed.copy()
        wrapped.translate_sites(range(len(seed)), [0, 0, 0.5], frac_coords=True, to_unit_cell=True)
        with self.assertRaisesRegex(ValueError, "unwrap or review"):
            change_vacuum(wrapped, 25)

    def test_beta_complete_layers_and_stable_face_identity(self):
        parent, seed, metadata = seed_fixture()
        for count in (2, 4, 6):
            result, record = extend_thickness(parent, seed, metadata, n_layers=count)
            self.assertEqual(len(result), 7 * count)
            self.assertEqual(record["n_layers"], count)
            self.assertEqual(record["oriented_repeat_count"], count)
            self.assertEqual(record["face_context"]["face_ids"]["upper"], metadata["upper_face_id"])
            self.assertEqual(record["face_context"]["face_ids"]["lower"], metadata["lower_face_id"])
            np.testing.assert_allclose(result.lattice.matrix[:2], seed.lattice.matrix[:2], atol=1e-12)
            self.assertAlmostEqual(record["vacuum_A"], 20)
            self.assertEqual(len(record["transformation"]["parent_site_indices"]), len(result))
            transform = record["transformation"]
            recovered = result.cart_coords + np.asarray(transform["inplane_wrap_integers"]) @ result.lattice.matrix[:2]
            recovered[:, 2] -= transform["translation_z_A"]
            recovered = recovered @ np.asarray(transform["row_cartesian_rotation"]).T
            recovered += metadata["raw_cut_shift"] * np.asarray(transform["oriented_unit_cell_lattice_A"])[2]
            mapped = parent.frac_coords[transform["parent_site_indices"]]
            difference = parent.lattice.get_fractional_coords(recovered) - mapped
            np.testing.assert_allclose(difference - np.rint(difference), 0, atol=1e-7)

    def test_iib_preserves_two_layer_stacking_repeat(self):
        parent, seed, metadata = seed_fixture("IIb")
        for repeats in (1, 2, 3):
            result, record = extend_thickness(parent, seed, metadata, oriented_repeats=repeats)
            self.assertEqual(record["n_layers"], 2 * repeats)
            self.assertEqual(len(result), 14 * repeats)
        with self.assertRaises(ValueError):
            extend_thickness(parent, seed, metadata, n_layers=3)

    def test_edges_round_up_without_trimming(self):
        parent, seed, metadata = seed_fixture("beta", "100", "t02")
        original_count = len(seed) // metadata["transformation"]["oriented_repeat_count"]
        for target in (12, 18, 24):
            result, record = extend_thickness(parent, seed, metadata, target_thickness_A=target)
            self.assertGreaterEqual(record["slab_thickness_A"] + 1e-8, target)
            self.assertEqual(len(result), original_count * record["oriented_repeat_count"])
            self.assertIsNone(record["n_layers"])
            repeated, second = extend_thickness(parent, seed, metadata, target_thickness_A=target)
            np.testing.assert_array_equal(result.cart_coords, repeated.cart_coords)
            self.assertEqual(record, second)

    def test_other_supported_edge_families_keep_their_terminations(self):
        for phase, facet, termination in (("beta", "110", "t01"), ("IIb", "100", "t01")):
            parent, seed, metadata = seed_fixture(phase, facet, termination)
            result, record = extend_thickness(parent, seed, metadata, target_thickness_A=24)
            self.assertGreaterEqual(record["slab_thickness_A"], 24 - 1e-8)
            self.assertEqual(record["termination_id"], metadata["termination_id"])
            self.assertEqual(result.composition.get_el_amt_dict()["S"], 4 * result.composition.get_el_amt_dict()["Zn"])

    def test_thickness_rejects_changed_geometry_and_multiple_targets(self):
        parent, seed, metadata = seed_fixture()
        changed = seed.copy()
        changed.translate_sites([0], [0.02, 0, 0], frac_coords=False)
        with self.assertRaisesRegex(ValueError, "relaxed geometry"):
            extend_thickness(parent, changed, metadata, n_layers=2)
        with self.assertRaisesRegex(ValueError, "exactly one"):
            extend_thickness(parent, seed, metadata, n_layers=2, oriented_repeats=2)
        vac, _ = change_vacuum(seed, 30)
        _, record = extend_thickness(parent, vac, metadata, n_layers=2)
        self.assertAlmostEqual(record["vacuum_A"], 30)


class SurfaceLayerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent, cls.thin, cls.metadata = seed_fixture()
        cls.thick, record = extend_thickness(cls.parent, cls.thin, cls.metadata, n_layers=4)
        cls.context = record["face_context"]

    def test_one_layer_does_not_acquire_bulk_like_interior(self):
        grouping = group_surface_layers(self.thin, face_context=self.context)
        self.assertFalse(grouping["interior_available"])
        self.assertEqual(grouping["groups"]["interior"], [])
        with self.assertRaisesRegex(ValueError, "central region"):
            propose_relaxed_region(self.thin, mode="central_fixed", explicitly_requested=True,
                                   face_context=self.context)

    def test_partition_and_normal_vacuum_translation_invariance(self):
        grouping = group_surface_layers(self.thick, face_context=self.context)
        self.assertEqual(grouping["status"], "GEOMETRY_GROUPED")
        self.assertTrue(grouping["interior_available"])
        self.assertIn("unverified", grouping["interior_interpretation"])
        assigned = [index for group in grouping["groups"].values() for index in group]
        self.assertEqual(sorted(assigned), list(range(len(self.thick))))
        moved, _ = change_vacuum(self.thick, 30)
        second = group_surface_layers(moved, face_context=self.context)
        self.assertEqual(grouping["groups"], second["groups"])

    def test_fixed_mask_requires_request_and_matching_bulk_like_evidence(self):
        grouping = group_surface_layers(self.thick, face_context=self.context)
        evidence = dict(bulk_like_established=True, reference="synthetic reviewed criterion",
                        criterion="synthetic fixture only", surface_id=self.context["surface_id"],
                        site_indices=grouping["groups"]["interior"], geometry_id=grouping["geometry_id"])
        with self.assertRaisesRegex(ValueError, "explicitly requested"):
            propose_relaxed_region(self.thick, mode="full")
        with self.assertRaisesRegex(ValueError, "bulk-like evidence"):
            propose_relaxed_region(self.thick, mode="central_fixed", explicitly_requested=True,
                                   grouping=grouping, face_context=self.context)
        mask, record = propose_relaxed_region(self.thick, mode="central_fixed", explicitly_requested=True,
                                             interior_evidence=evidence, grouping=grouping,
                                             face_context=self.context)
        self.assertEqual(np.asarray(mask).shape, (len(self.thick), 3))
        self.assertEqual(np.flatnonzero(~np.asarray(mask).all(axis=1)).tolist(),
                         sorted(grouping["groups"]["interior"]))
        self.assertEqual(record["scientific_status"], "PROPOSED_REQUIRES_HUMAN_REVIEW")
        evidence["surface_id"] = "unrelated"
        with self.assertRaises(ValueError):
            propose_relaxed_region(self.thick, mode="central_fixed", explicitly_requested=True,
                                   interior_evidence=evidence, grouping=grouping, face_context=self.context)

    def test_geometry_identity_binds_evidence_even_when_group_membership_is_unchanged(self):
        grouping = group_surface_layers(self.thick, face_context=self.context)
        evidence = dict(bulk_like_established=True, reference="synthetic reviewed criterion",
                        criterion="synthetic fixture only", surface_id=self.context["surface_id"],
                        site_indices=grouping["groups"]["interior"], geometry_id=grouping["geometry_id"])
        changed = self.thick.copy()
        changed.translate_sites([grouping["groups"]["interior"][0]], [0.01, 0, 0], frac_coords=False)
        updated = group_surface_layers(changed, face_context=self.context)
        self.assertEqual(updated["groups"], grouping["groups"])
        self.assertNotEqual(updated["geometry_id"], grouping["geometry_id"])
        with self.assertRaisesRegex(ValueError, "no longer matches"):
            propose_relaxed_region(changed, mode="central_fixed", explicitly_requested=True,
                                   interior_evidence=evidence, grouping=grouping, face_context=self.context)
        with self.assertRaisesRegex(ValueError, "bulk-like evidence"):
            propose_relaxed_region(changed, mode="central_fixed", explicitly_requested=True,
                                   interior_evidence=evidence, grouping=updated, face_context=self.context)
        annotated = self.thick.copy()
        annotated.add_site_property("selective_dynamics", [[True, True, True]] * len(annotated))
        self.assertEqual(group_surface_layers(annotated)["geometry_id"], grouping["geometry_id"])

    def test_full_mask_is_explicit_and_stale_grouping_is_rejected(self):
        mask, _ = propose_relaxed_region(self.thin, mode="full", explicitly_requested=True)
        self.assertTrue(np.asarray(mask).all())
        grouping = group_surface_layers(self.thick, face_context=self.context)
        grouping["groups"]["interior"] = []
        with self.assertRaisesRegex(ValueError, "no longer matches"):
            propose_relaxed_region(self.thick, mode="central_fixed", explicitly_requested=True,
                                   grouping=grouping, face_context=self.context)


if __name__ == "__main__":
    unittest.main()
