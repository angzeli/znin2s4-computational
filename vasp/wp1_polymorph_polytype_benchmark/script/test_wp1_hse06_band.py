"""Focused archival checks on the completed beta XML and temporary exports only."""
from pathlib import Path
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np

import export_wp1_hse06_band as hse
from plot_wp1_bands import LABEL_FONT, LABEL_SIZE, LABEL_Y, display_distances, path_ticks, plot_bands
from wp1_electronic_plot_data import load_phase


class HseArchiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, cls.metadata = hse.read_source()

    def test_real_optional_source_shape_and_provenance(self):
        self.assertEqual(self.data.raw_energy.shape, (180, 40))
        self.assertEqual(self.data.reference_ev, 4.0196)
        self.assertEqual(self.metadata["edges"]["sampled_cbm_raw_eV"], 5.1103)
        self.assertEqual(self.metadata["edges"]["fundamental_gap_eV"], 1.0907)
        self.assertEqual(self.metadata["edges"]["minimum_direct_gap_eV"], 1.0919)
        self.assertEqual(self.metadata["edges"]["direct_minus_fundamental_meV"], 1.2)
        self.assertEqual(self.metadata["source_files"]["vasprun.xml"]["sha256"],
                         "73d857dcc2a9f161071e87a4d2b459af350b3f5bc2db33b292f30f29cbd5dcde")
        serialized = json.dumps(self.metadata, allow_nan=False)
        self.assertNotIn("/Users/", serialized)
        self.assertNotIn("/Volumes/", serialized)

    def test_fail_closed_on_missing_partial_and_wrong_spin_xml(self):
        source = hse.ROOT / hse.SOURCE / "vasprun.xml"
        for mutation in ("missing", "partial", "spin"):
            tree = ET.parse(source)
            optional = tree.find("calculation/eigenvalues_kpoints_opt")
            if mutation == "missing":
                tree.find("calculation").remove(optional)
            elif mutation == "partial":
                spin = optional.find("eigenvalues/array/set/set")
                spin.remove(spin[-1])
            else:
                tree.find("parameters//i[@name='ISPIN']").text = "2"
            with self.subTest(mutation=mutation), patch.object(hse.ET, "parse", return_value=tree):
                with self.assertRaises(ValueError):
                    hse.read_source()

    def test_fail_closed_on_changed_gap_nonfinite_and_upper_window(self):
        raw = self.data.raw_energy.copy()
        raw[:, 31:] += .1
        with self.assertRaisesRegex(ValueError, "gaps differ"):
            hse.band_edges(raw)
        raw = self.data.raw_energy.copy()
        raw[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "finite"):
            hse.band_edges(raw)
        raw = self.data.raw_energy.copy()
        raw[0, 32:] = [6., 6.1, 6.2, 6.3, 6.4, 7.9, 10.1, 10.2]
        with self.assertRaisesRegex(ValueError, "problematic band"):
            hse.band_edges(raw)

    def test_csv_round_trip_all_rows_common_zero_and_repeated_points(self):
        with tempfile.TemporaryDirectory() as directory:
            csv_path = Path(directory) / hse.CSV_NAME
            hse.write_csv(csv_path, self.data)
            with csv_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 7200)
            self.assertEqual({int(r["path_point_index"]) for r in rows}, set(range(1, 181)))
            self.assertEqual({int(r["band_index"]) for r in rows}, set(range(1, 41)))
            # The portable path must not consult XML, geometry or JSON on reload.
            with patch.object(hse.ET, "parse", side_effect=AssertionError("Native XML read")), \
                 patch.object(hse.accepted, "parse_poscar", side_effect=AssertionError("Geometry read")):
                restored = hse.read_csv(csv_path)
                hse.save_pdf(restored, Path(directory) / hse.PDF_NAME)
            np.testing.assert_array_equal(restored.raw_energy, self.data.raw_energy)
            np.testing.assert_array_equal(restored.kpoints, self.data.kpoints)
            np.testing.assert_array_equal(restored.path.distances, self.data.path.distances)
            np.testing.assert_allclose(restored.band_energy, self.data.band_energy, atol=1e-12, rtol=0)
            self.assertEqual(hse.band_edges(restored.raw_energy), self.metadata["edges"])
            pairs, checks = hse.repeated_points(restored.kpoints, restored.raw_energy)
            self.assertEqual(len(pairs), 18)
            self.assertEqual(checks["maximum_difference_by_band_eV"],
                             {"31": 0., "32": 0., "38": .0005, "39": .0048, "40": .1421})
            self.assertEqual(checks["minimum_energy_minus_vbm_bands_38_to_40_eV"], 5.4992)
            points = rows[::40]
            self.assertEqual([int(points[i]["branch_index"]) for i in (0, 139, 140, 159, 160, 179)],
                             [1, 1, 2, 2, 3, 3])
            self.assertEqual(restored.path.disconnected_before, (7, 8))
            for segment in restored.path.segments:
                self.assertEqual(segment.stop - segment.start, 20)
            for index in (140, 160):
                self.assertEqual(restored.path.distances[index], restored.path.distances[index - 1])
                self.assertGreater(float(points[index]["display_x"]), float(points[index - 1]["display_x"]))

    def test_csv_rejects_row_loss_and_inconsistent_shift(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / hse.CSV_NAME
            hse.write_csv(path, self.data)
            text = path.read_text()
            path.write_text("\n".join(text.splitlines()[:-1]) + "\n")
            with self.assertRaisesRegex(ValueError, "7200"):
                hse.read_csv(path)
            rows = list(csv.reader(text.splitlines()))
            rows[1][-1] = str(float(rows[1][-1]) + .1)
            with path.open("w", newline="") as handle:
                csv.writer(handle).writerows(rows)
            with self.assertRaisesRegex(ValueError, "common VBM"):
                hse.read_csv(path)

    def test_figure_style_no_bridging_and_no_clipped_labels(self):
        fig, ax, info = plot_bands(self.data, window=hse.WINDOW)
        try:
            np.testing.assert_array_equal(fig.get_size_inches(), [7, 5])
            self.assertEqual(ax.get_ylim(), (-4, 4))
            self.assertEqual(ax.get_title(), "")
            self.assertEqual(ax.get_xlabel(), "")
            self.assertIn("VBM", ax.get_ylabel())
            self.assertFalse(fig.texts)
            self.assertEqual(info["tick_labels"], ["Γ", "M", "K", "Γ", "A", "L", "H", "A", "L", "M", "H", "K"])
            self.assertEqual(list(ax.texts), info["label_artists"])
            self.assertEqual(info["display_branch_gap_fraction"], .055)
            self.assertTrue(all(b < 37 for b in info["band_indices"]))
            shown = display_distances(self.data)
            for i, segment in enumerate(self.data.path.segments):
                for j, band in enumerate(info["band_indices"]):
                    line = ax.lines[i * len(info["band_indices"]) + j]
                    np.testing.assert_array_equal(line.get_xdata(), shown[segment])
                    np.testing.assert_array_equal(line.get_ydata(), self.data.band_energy[segment, band])
                    self.assertEqual(len(line.get_xdata()), 20)
                    self.assertEqual(line.get_linewidth(), 1.1)
            fig.canvas.draw()
            boxes = []
            for label in info["label_artists"]:
                self.assertEqual(label.get_fontfamily(), [LABEL_FONT])
                self.assertEqual(label.get_fontsize(), LABEL_SIZE)
                self.assertEqual(label.get_fontweight(), "bold")
                self.assertEqual(label.get_position()[1], LABEL_Y)
                self.assertEqual(label.get_rotation(), 0)
                box = label.get_window_extent()
                self.assertTrue(fig.bbox.contains(box.x0, box.y0) and fig.bbox.contains(box.x1, box.y1))
                boxes.append(box)
            self.assertFalse(any(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1:]))
        finally:
            plt.close(fig)

    def test_deterministic_csv_pdf_and_exact_pdf_page(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            for i in range(2):
                hse.write_csv(path / f"{i}.csv", self.data)
                hse.save_pdf(hse.read_csv(path / f"{i}.csv"), path / f"{i}.pdf")
            self.assertEqual((path / "0.csv").read_bytes(), (path / "1.csv").read_bytes())
            pdf = (path / "0.pdf").read_bytes()
            self.assertEqual(pdf, (path / "1.pdf").read_bytes())
            box = re.search(rb"/MediaBox\s*\[([^]]+)\]", pdf)
            self.assertEqual(list(map(float, box[1].split())), [0., 0., 504., 360.])
            _, again = hse.read_source()
            self.assertEqual(again, self.metadata)

    def test_frozen_beta_pbe_data_preparation_regression(self):
        pbe = load_phase(hse.ROOT, "beta", reference="vbm")
        self.assertEqual(pbe.band_energy.shape, (360, 40))
        self.assertEqual(pbe.endpoint_labels, self.data.endpoint_labels)
        self.assertEqual(pbe.path.disconnected_before, (7, 8))
        np.testing.assert_allclose(path_ticks(pbe)[0], path_ticks(self.data)[0], atol=1e-8, rtol=0)
        self.assertAlmostEqual(pbe.summary["path_gap_eV"], .303118, places=6)
        self.assertAlmostEqual(pbe.band_energy[:, 30].max(), 0., places=6)

    def test_import_and_help_have_no_analysis_or_output(self):
        script = Path(hse.__file__).resolve()
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
                           PYTHONPATH=os.pathsep.join((str(script.parent), os.environ.get("PYTHONPATH", ""))))
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable, "-B", "-c",
                "import export_wp1_hse06_band as h; "
                "h.read_source=lambda *a: (_ for _ in ()).throw(AssertionError('analysis')); "
                "h.main(['--help'])"], cwd=directory, env=environment,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            self.assertEqual(list(Path(directory).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
