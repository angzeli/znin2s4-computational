"""Focused Stage 0 provenance, geometry, registry, and regeneration checks."""

import contextlib
import importlib
import io
import json
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar

import generate_surface_candidates as generate
from surface_model_utils import (
    PLANE_TOL, VACUUM, bond_graph, connectivity, diagnose_slab,
    faces_equivalent, present_slab, slabs_match,
)


class StageZeroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.files, cls.rows, cls.faces, cls.parents = generate.build_candidates()

    def test_parent_bytes_and_audited_handoff(self):
        for phase, metadata in self.parents.items():
            self.assertTrue(metadata["downstream_match"])
            data = self.files[Path("structure/bulk_parents") / phase / "POSCAR"]
            self.assertEqual(generate.sha256(data), generate.ACCEPTED[phase]["sha256"])
            self.assertEqual(data, (generate.REPO / metadata["source_wp1_path"]).read_bytes())

    def test_registries_and_duplicate_references(self):
        lookup = {r["slab_id"]: r for r in self.rows}
        self.assertEqual(len(lookup), len(self.rows))
        materialized = {r["slab_id"] for r in self.rows if r["materialized_path"]}
        self.assertEqual(Counter(f["slab_id"] for f in self.faces), Counter({s: 2 for s in materialized}))
        self.assertEqual(len({f["face_id"] for f in self.faces}), len(self.faces))
        for row in self.rows:
            if row["duplicate_of"]:
                self.assertIn(row["duplicate_of"], lookup)
                self.assertFalse(lookup[row["duplicate_of"]]["duplicate_of"])
                self.assertFalse(row["materialized_path"])
            if row["materialized_path"]:
                self.assertIn(row["status"], ("PASS", "HOLD"))
            self.assertEqual(row["broken_ZnS_estimate"], row["cut_diagnostics"]["crossing_ZnS"])
            self.assertEqual(row["broken_InS_estimate"], row["cut_diagnostics"]["crossing_InS"])
        for phase, facets in generate.FACETS.items():
            for facet in facets:
                group = [r for r in self.rows if r["phase"] == phase and r["facet_label"] == f"({facet})"]
                self.assertTrue(group)
                self.assertEqual([r["raw_cut_shift"] for r in group], sorted(r["raw_cut_shift"] for r in group))
                self.assertEqual([r["slab_id"] for r in group],
                                 [f"{phase}_{facet}_t{i:02d}" for i in range(1, len(group) + 1)])

    def test_retained_geometry_transform_and_vacuum(self):
        for row in self.rows:
            if not row["materialized_path"]:
                continue
            path = (generate.REPO / row["materialized_path"]).relative_to(generate.WP2)
            slab = Poscar.from_str(self.files[path].decode()).structure
            metadata = json.loads(self.files[path.with_name("metadata.json")])
            self.assertEqual(generate.sha256(self.files[path]), metadata["structure_sha256"])
            self.assertFalse(metadata["diagnostics"]["validation_errors"])
            self.assertNotIn("Selective", self.files[path].decode())
            self.assertAlmostEqual(slab.lattice.c - np.ptp(slab.cart_coords[:, 2]), VACUUM)
            self.assertGreater(slab.lattice.c, VACUUM)
            self.assertAlmostEqual(slab.cart_coords[:, 2].min(), VACUUM / 2)
            transform = row["transformation"]
            rotation = np.array(transform["row_cartesian_rotation"])
            np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-12)
            self.assertAlmostEqual(np.linalg.det(rotation), 1)
            # Undo display wrapping/translation/rotation and the original cut.
            coords = slab.cart_coords + np.array(transform["inplane_wrap_integers"]) @ slab.lattice.matrix[:2]
            coords[:, 2] -= transform["translation_z_A"]
            coords = coords @ rotation.T
            coords += row["raw_cut_shift"] * np.array(transform["oriented_unit_cell_lattice_A"])[2]
            parent, _, _ = generate.parent_record(row["phase"])
            parent_coords = parent.lattice.get_fractional_coords(coords)
            expected = parent.frac_coords[transform["parent_site_indices"]]
            difference = parent_coords - expected
            np.testing.assert_allclose(difference - np.rint(difference), 0, atol=1e-7)
            expected_plane = np.array(transform["oriented_unit_cell_lattice_A"])[:2] @ rotation
            np.testing.assert_allclose(slab.lattice.matrix[:2], expected_plane, atol=1e-10)

    def test_basal_completeness_and_bad_cut_rejection(self):
        for phase, layer_count in (("beta", 1), ("IIb", 2)):
            basal = [r for r in self.rows if r["phase"] == phase and r["facet_label"] == "(001)"]
            retained = [r for r in basal if r["materialized_path"]]
            self.assertEqual(len(retained), 1)
            self.assertEqual(retained[0]["complete_layer_count"], layer_count)
            self.assertTrue(retained[0]["cut_diagnostics"]["no_first_shell_bond_crosses"])
            self.assertGreater(retained[0]["cut_diagnostics"]["nearest_gap_width_A"], 2.7)
            bad = [r for r in basal if not r["cut_diagnostics"]["no_first_shell_bond_crosses"]]
            self.assertTrue(bad)
            self.assertTrue(all(r["status"] == "REJECT" for r in bad))

    def test_reversal_and_translation_do_not_merge_inequivalent_faces(self):
        parent, _, metadata = generate.parent_record("beta")
        generator = generate.SlabGenerator(parent, (0, 0, 1), **generate.construction_settings(True))
        shift = next(r["raw_cut_shift"] for r in self.rows if r["slab_id"] == "beta_001_t07")
        slab, _ = present_slab(generator.get_slab(shift), generator)
        coords = slab.frac_coords.copy()
        coords[:, 2] = 1 - coords[:, 2]
        reverse = Structure(slab.lattice, slab.species, coords, site_properties=slab.site_properties)
        criteria = metadata["coordination_criterion"]
        d = diagnose_slab(slab, parent, criteria, True)
        rd = diagnose_slab(reverse, parent, criteria, True)
        self.assertTrue(slabs_match(slab, reverse, d, rd))
        self.assertFalse(faces_equivalent(slab))
        self.assertNotEqual(d["faces"]["upper"]["surface_species_fingerprint"],
                            d["faces"]["lower"]["surface_species_fingerprint"])
        translated = slab.copy()
        translated.translate_sites(range(len(slab)), [0.237, 0.129, 0])
        self.assertTrue(slabs_match(slab, translated, d, d))
        # Regression: close spinel planes must not defeat normal-reversal matching.
        parent, _, metadata = generate.parent_record("spinel")
        generator = generate.SlabGenerator(parent, (0, 0, 1), **generate.construction_settings(False))
        shift = generator.gen_possible_terminations(ftol=PLANE_TOL)[0]
        slab, _ = present_slab(generator.get_slab(shift), generator)
        coords = slab.frac_coords.copy()
        coords[:, 2] = 1 - coords[:, 2]
        reverse = Structure(slab.lattice, slab.species, coords, site_properties=slab.site_properties)
        criteria = metadata["coordination_criterion"]
        self.assertTrue(slabs_match(slab, reverse, diagnose_slab(slab, parent, criteria, False),
                                    diagnose_slab(reverse, parent, criteria, False)))

    def test_finite_fragments_are_distinguished_from_periodic_components(self):
        finite = Structure(Lattice.cubic(20), ["Zn", "S"], [[0, 0, 0], [2.4, 0, 0]], coords_are_cartesian=True)
        criteria = {"Zn": {"cutoff_A": 2.9}, "In": {"cutoff_A": 2.9}}
        components = connectivity(finite, bond_graph(finite, criteria))
        self.assertEqual([(c["atom_count"], c["periodic_rank"]) for c in components], [(2, 0)])
        edges = [r for r in self.rows if r["phase"] == "beta" and r["facet_label"] == "(100)"
                 and r["materialized_path"]]
        self.assertTrue(all(c["periodic_rank"] == 1 for r in edges for c in r["connectivity_components"]))

    def test_deterministic_regeneration_and_existing_file_guard(self):
        regenerated, rows, faces, _ = generate.build_candidates()
        self.assertEqual(regenerated, self.files)
        self.assertEqual(rows, self.rows)
        self.assertEqual(faces, self.faces)
        with tempfile.TemporaryDirectory(prefix="wp2-stage0-test-") as tmp:
            root = Path(tmp)
            generate.write_outputs(self.files, root)
            for path, data in self.files.items():
                self.assertEqual((root / path).read_bytes(), data)
            self.assertFalse((root / "calculation").exists())
            protected = root / "structure/bulk_parents/beta/POSCAR"
            protected.write_text("different pre-existing parent\n")
            with self.assertRaises(FileExistsError):
                generate.write_outputs(self.files, root)
            self.assertEqual(protected.read_text(), "different pre-existing parent\n")

    def test_import_and_help_do_not_write(self):
        with patch.object(Path, "write_bytes", side_effect=AssertionError("import write")), \
             patch.object(Path, "write_text", side_effect=AssertionError("import write")), \
             patch.object(Path, "mkdir", side_effect=AssertionError("import mkdir")):
            importlib.reload(generate)
            with patch.object(sys, "argv", ["generate_surface_candidates.py", "--help"]), \
                 contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as result:
                generate.main()
            self.assertEqual(result.exception.code, 0)


if __name__ == "__main__":
    unittest.main()
