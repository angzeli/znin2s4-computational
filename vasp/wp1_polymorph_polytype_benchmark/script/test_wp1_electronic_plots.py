"""Focused synthetic plotting contracts and read-only accepted-data integration.

Run sequentially with the plotting environment and numerical threads set to one.
Temporary figure exports are removed by TemporaryDirectory, never written to runs.
"""
from pathlib import Path
import os
import re
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from cmw.analysis.electronic import prepare_path
from plot_wp1_bands import BRANCH_GAPS, PATH_FAMILY, LABEL_Y, LABEL_SIZE, LABEL_FONT, display_distances, path_ticks, plot_bands
from plot_wp1_dos_pdos import (DEFAULT_WINDOW, FAMILIES, adaptive_y_limits,
                              common_y_limits, families_for_mode, plot_dos)
from wp1_electronic_plot_data import (EXPECTED_FU, GROUPS, PHASES, SELECTED, WP1,
    accepted, energy_window, load_all, load_phase, prepared_dos, project_root, source_mapping)
from wp1_electronic_plot_style import save_figure


class PlotContractTests(unittest.TestCase):
    def test_disconnected_single_axes_export_and_style(self):
        ends = np.array([[[0, 0, 0], [1.5, 0, 0]], [[0, 0, .5], [0, 0, 0]]])
        coordinates = np.concatenate([np.linspace(a, b, 4) for a, b in ends])
        path = prepare_path(np.eye(3), ends, coordinates % 1, 4)
        data = SimpleNamespace(phase='spinel', reference_mode='fermi', reference_ev=2.,
            band_energy=np.column_stack((np.linspace(-1, 1, 8), np.ones(8)*4)),
            path=path, endpoint_labels=(('GAMMA', 'U'), ('K', 'GAMMA')),
            summary={'primitive_cells': 4})
        before = dict(mpl.rcParams)
        fig, ax, metadata = plot_bands(data, window=(-2, 2))
        self.assertEqual(dict(mpl.rcParams), before)
        self.assertEqual(len(fig.axes), 1)
        self.assertEqual(ax.get_xlabel(), '')
        self.assertEqual(ax.get_title(), '')
        self.assertEqual(list(ax.texts), metadata['label_artists'])
        self.assertFalse(fig.texts)
        np.testing.assert_array_equal(fig.get_size_inches(), [7, 5])
        self.assertFalse(ax.get_xticklabels())
        for label in metadata['label_artists']:
            self.assertEqual(label.get_position()[1], LABEL_Y)
            self.assertEqual(label.get_rotation(), 0)
            self.assertEqual(label.get_horizontalalignment(), 'center')
            self.assertEqual(label.get_verticalalignment(), 'baseline')
            self.assertEqual(label.get_fontsize(), LABEL_SIZE)
            self.assertEqual(label.get_fontweight(), 'bold')
            self.assertEqual(label.get_fontfamily(), [LABEL_FONT])
        np.testing.assert_array_equal([label.get_position()[0] for label in metadata['label_artists']], ax.get_xticks())
        self.assertIn('folded-cell', metadata['caption'])
        self.assertEqual(metadata['tick_labels'], ['Γ', 'U', 'K', 'Γ'])
        np.testing.assert_array_equal(metadata['distances'], path.distances)
        shown = metadata['display_distances']
        self.assertFalse(np.shares_memory(shown, path.distances))
        self.assertAlmostEqual((shown[4] - shown[3]) / np.ptp(path.distances), .0225)
        self.assertEqual(metadata['band_indices'].tolist(), [0])
        band_lines = ax.lines[:2]
        for line, segment in zip(band_lines, path.segments):
            np.testing.assert_array_equal(line.get_xdata(), shown[segment])
            np.testing.assert_allclose(np.diff(shown[segment]), np.diff(path.distances[segment]), atol=1e-14)
            np.testing.assert_array_equal(line.get_ydata(), data.band_energy[segment, 0])
            self.assertEqual(len(line.get_xdata()), 4)
            self.assertEqual(line.get_linewidth(), 1.1)
        with tempfile.TemporaryDirectory() as directory:
            result = save_figure(fig, directory, 'spinel', 'band_structure')
            self.assertEqual([p.name for p in Path(directory).iterdir()], ['spinel_pbe_band_structure.pdf'])
            self.assertGreater(result.stat().st_size, 1000)
            self.assertEqual(len(re.findall(rb'/Type\s*/Page\b', result.read_bytes())), 1)
            with self.assertRaises(FileExistsError):
                save_figure(fig, directory, 'spinel', 'band_structure')
            with self.assertRaises(ValueError):
                save_figure(fig, directory, '../spinel', 'band_structure')
            from PIL import Image
            png = save_figure(fig, directory, 'spinel', 'band_structure', output_format='png')
            with Image.open(png) as raster:
                np.testing.assert_allclose(raster.info['dpi'], (600, 600), atol=.01)
                self.assertEqual(raster.convert('RGBA').getextrema()[3], (255, 255))
                self.assertEqual(raster.convert('RGB').getpixel((0, 0)), (255, 255, 255))
        plt.close(fig)

    def test_publication_inventory_and_displayed_only_y_scale(self):
        self.assertEqual(FAMILIES, ('selected_orbital_pdos',))
        self.assertEqual(families_for_mode('publication'), FAMILIES)
        self.assertEqual(families_for_mode('element-pdos'), ('element_pdos',))
        self.assertEqual(len(PHASES) * (len(FAMILIES) + 1), 10)
        self.assertEqual(DEFAULT_WINDOW, (-3, 3))
        data = SimpleNamespace(phase='beta', reference_mode='vbm', reference_ev=4.,
            formula_units=2, dos_energy=np.array([-6., -3., 0., 3., 4.]),
            curves={name: np.array([1000., 2., 0., 4., 500.])
                    for name in (*SELECTED, 'TDOS', 'Zn', 'In', 'S')})
        for kind in ('selected_orbital_pdos', 'element_pdos'):
            limits = adaptive_y_limits(data, kind, DEFAULT_WINDOW)
            self.assertGreaterEqual(limits[1], 1.05 * 2)
            self.assertLessEqual(limits[1], 1.10 * 2)
            fig, ax, metadata = plot_dos(data, kind)
            self.assertEqual(ax.get_ylim(), limits)
            self.assertFalse(ax.texts)
            self.assertFalse(fig.texts)
            self.assertEqual(ax.get_title(), '')
            self.assertNotIn('Energy,', ax.get_xlabel())
            self.assertIn('VBM', ax.get_xlabel())
            for name, values in metadata['curves'].items():
                np.testing.assert_array_equal(values, data.curves[name] / 2)
            plt.close(fig)

    def test_dos_exact_samples_single_channel_no_rc_leak(self):
        energy = np.array([-1., 0, 1.])
        density = np.array([2., 0, 4.])
        data = SimpleNamespace(phase='alpha1', reference_mode='fermi', reference_ev=3.,
            formula_units=3, dos_energy=energy, curves={'TDOS': density})
        before = dict(mpl.rcParams)
        fig, ax, metadata = plot_dos(data, 'tdos', window=(-1, 1), ylim=(0, 2))
        self.assertEqual(dict(mpl.rcParams), before)
        self.assertEqual(len(fig.axes), 1)
        self.assertEqual(tuple(metadata['curves']), ('TDOS',))
        np.testing.assert_array_equal(ax.lines[0].get_xdata(), energy)
        np.testing.assert_array_equal(ax.lines[0].get_ydata(), density / 3)
        np.testing.assert_array_equal(data.curves['TDOS'], [2, 0, 4])
        plt.close(fig)

    def test_import_and_root_resolution(self):
        script = str(Path(__file__).resolve().parent)
        # Preserve an explicitly supplied dependency checkout in child probes.
        python_path = os.pathsep.join(filter(None, (script, os.environ.get('PYTHONPATH', ''))))
        environment = dict(os.environ, PYTHONPATH=python_path, PYTHONDONTWRITEBYTECODE='1')
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable, '-B', '-c',
                "import pathlib,matplotlib as m; before=dict(m.rcParams); "
                "import plot_wp1_dos_pdos,plot_wp1_bands; "
                "assert before==dict(m.rcParams); assert not list(pathlib.Path('.').iterdir())"],
                cwd=directory, env=environment, check=True)
        root = project_root()
        for cwd in (root, root / WP1 / 'results/04_electronic_structure'):
            found = subprocess.check_output([sys.executable, '-B', '-c',
                'from wp1_electronic_plot_data import project_root; print(project_root())'],
                cwd=cwd, env=environment, text=True).strip()
            self.assertEqual(found, str(root))
        with self.assertRaises(ValueError):
            source_mapping(root, ('beta', 'beta'))

    def test_element_legend_and_subscript_rendering(self):
        data = SimpleNamespace(phase='alpha1', reference_mode='fermi', reference_ev=0.,
            formula_units=3, dos_energy=np.array([-1., 0, 1.]),
            curves={name: np.array([2., 0, 4.]) for name in ('TDOS', 'Zn', 'In', 'S')})
        fig, ax, _ = plot_dos(data, 'element_pdos', window=(-1, 1), ylim=(0, 2))
        with tempfile.TemporaryDirectory() as directory:
            save_figure(fig, directory, 'alpha1', 'element_pdos')
            legend = ax.get_legend().get_window_extent()
            self.assertGreaterEqual(legend.y0, ax.bbox.y0)
            self.assertLessEqual(legend.y1, ax.bbox.y1)
        plt.close(fig)


class RealDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = project_root()
        cls.data = load_all(cls.root)

    def test_mapping_gaps_references_and_normalization(self):
        self.assertEqual(tuple(d.phase for d in self.data), tuple(PHASES))
        for data in self.data:
            self.assertEqual(data.formula_units, EXPECTED_FU[data.phase])
            self.assertEqual(data.summary['validation'], 'PASS')
            raw_band = np.asarray(accepted.parse_eigenval(data.sources['band'] / 'EIGENVAL').energies)
            np.testing.assert_allclose(data.band_energy + data.reference_ev, raw_band, atol=1e-12)
            dos = accepted.parse_doscar(data.sources['dos'] / 'DOSCAR',
                                       accepted.parse_poscar(data.sources['dos'] / 'POSCAR').site_elements)
            self.assertEqual(data.reference_ev, dos.efermi)
            np.testing.assert_allclose(data.dos_energy + data.reference_ev, dos.energies, atol=.0011, rtol=0)
            _, normalized = prepared_dos(data)
            for name in data.curves:
                np.testing.assert_array_equal(normalized[name], data.curves[name] / data.formula_units)
                self.assertEqual(data.curves[name].shape, (3000,))
            self.assertTrue(set(SELECTED) <= set(data.curves))
            self.assertEqual(len(GROUPS['p']), 3)
            self.assertEqual(len(GROUPS['d']), 5)

    def test_paths_and_shared_plot_limits(self):
        dos_window, band_window = energy_window(self.data, 'dos'), energy_window(self.data, 'band')
        limits = common_y_limits(self.data, dos_window)
        self.assertEqual(set(limits), {'element_pdos', 'selected_orbital_pdos'})
        for data in self.data:
            self.assertEqual(len(data.path.distances), len(data.path.segments)*40)
            self.assertTrue(np.isfinite(data.band_energy).all())
            self.assertEqual(tuple(data.endpoint_labels), accepted.PHASE_META[data.phase]['segments'])
            _, labels = path_ticks(data)
            self.assertFalse(any(' | ' in label for label in labels))
            self.assertEqual(len(labels), len(data.path.segments) + 1 + len(data.path.disconnected_before))
            self.assertTrue(data.band_coverage[0] <= band_window[0] < band_window[1] <= data.band_coverage[1])

    def test_vbm_reference_preserves_arrays_sampling_and_interpretation(self):
        for fermi in self.data:
            vbm = load_phase(self.root, fermi.phase, reference='vbm')
            offset = fermi.reference_ev - vbm.reference_ev
            np.testing.assert_allclose(vbm.dos_energy, fermi.dos_energy + offset, atol=1e-12)
            np.testing.assert_allclose(vbm.band_energy, fermi.band_energy + offset, atol=1e-12)
            dense = accepted.parse_eigenval(vbm.sources['dos'] / 'EIGENVAL')
            raw = accepted.parse_eigenval(vbm.sources['band'] / 'EIGENVAL')
            dense_edges, path_edges = accepted.band_edges(dense), accepted.band_edges(raw)
            self.assertEqual(vbm.reference_ev, dense_edges['vbm'])
            expected = np.asarray(raw.energies) - dense_edges['vbm']
            np.testing.assert_array_equal(vbm.band_energy, expected)
            np.testing.assert_array_equal(vbm.path.distances, fermi.path.distances)
            np.testing.assert_array_equal(vbm.path.intended_fractional, fermi.path.intended_fractional)
            self.assertEqual(vbm.path.segments, fermi.path.segments)
            self.assertEqual(vbm.path.disconnected_before, fermi.path.disconnected_before)
            self.assertEqual(vbm.endpoint_labels, fermi.endpoint_labels)
            self.assertEqual(vbm.summary['primitive_cells'], 4 if vbm.phase == 'spinel' else 3 if vbm.phase == 'alpha1' else 1)
            curve_copies = {name: values.copy() for name, values in vbm.curves.items()}
            for kind, window in (('selected_orbital_pdos', DEFAULT_WINDOW), ('element_pdos', (-8, 4))):
                dos_fig, dos_ax, dos_metadata = plot_dos(vbm, kind, window=window)
                inside = (vbm.dos_energy >= window[0]) & (vbm.dos_energy <= window[1])
                maximum = max(values[inside].max() for values in dos_metadata['curves'].values())
                self.assertTrue(1.05 <= dos_ax.get_ylim()[1] / maximum <= 1.10)
                self.assertEqual(dos_metadata['reference_eV'], vbm.reference_ev)
                np.testing.assert_array_equal(dos_metadata['energy'], vbm.dos_energy)
                self.assertFalse(dos_ax.texts)
                for name, values in dos_metadata['curves'].items():
                    np.testing.assert_array_equal(values, fermi.curves[name] / vbm.formula_units)
                plt.close(dos_fig)
            for name in curve_copies:
                np.testing.assert_array_equal(vbm.curves[name], curve_copies[name])
            for field in ('uniform_gap_eV', 'path_gap_eV', 'classification'):
                self.assertEqual(vbm.summary[field], fermi.summary[field])
            self.assertAlmostEqual(path_edges['fundamental_gap'], vbm.summary['path_gap_eV'], places=6)
            self.assertEqual(energy_window([vbm], 'band', (-4, 4)), (-4, 4))
            if vbm.phase in ('spinel', 'IIa_prime'):
                occupied_index = int(raw.nelect / 2) - 1
                self.assertGreater(vbm.band_energy[:, occupied_index].max(), 0)
                self.assertAlmostEqual(vbm.band_energy[:, occupied_index].max(),
                                       path_edges['vbm'] - dense_edges['vbm'], places=12)
            if vbm.phase == 'beta':
                self.assertEqual(vbm.summary['classification'], 'AMBIGUOUS / NEAR-DEGENERATE')
            before = vbm.band_energy.copy()
            physical_before = vbm.path.distances.copy()
            fig, ax, metadata = plot_bands(vbm)
            np.testing.assert_array_equal(fig.get_size_inches(), [7, 5])
            self.assertEqual(ax.get_ylim(), (-4, 4))
            self.assertEqual((LABEL_SIZE, LABEL_Y, LABEL_FONT), (12.5, -.04, 'Arial Narrow'))
            fraction = {'spinel': .0225, 'alpha1': .0375, 'beta': .055, 'IIa_prime': .055, 'IIb': .055}[vbm.phase]
            self.assertEqual(BRANCH_GAPS[PATH_FAMILY[vbm.phase]], fraction)
            self.assertEqual(metadata['display_branch_gap_fraction'], fraction)
            fig.canvas.draw()
            label_artists = metadata['label_artists']
            self.assertEqual(list(ax.texts), label_artists)
            self.assertEqual([label.get_text() for label in label_artists], metadata['tick_labels'])
            self.assertEqual({label.get_position()[1] for label in label_artists}, {LABEL_Y})
            self.assertTrue(all(label.get_rotation() == 0 and label.get_horizontalalignment() == 'center'
                                and label.get_verticalalignment() == 'baseline'
                                and label.get_fontsize() == LABEL_SIZE and label.get_fontweight() == 'bold'
                                for label in label_artists))
            boxes = [label.get_window_extent() for label in label_artists]
            self.assertFalse(any(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i+1:]))
            self.assertTrue(all(box.x0 >= fig.bbox.x0 and box.x1 <= fig.bbox.x1
                                and box.y0 >= fig.bbox.y0 for box in boxes))
            np.testing.assert_array_equal(vbm.band_energy, before)
            np.testing.assert_array_equal(vbm.path.distances, physical_before)
            self.assertIn('VBM', ax.get_ylabel())
            offset = 0.
            gap_width = np.ptp(physical_before) * fraction
            for i, segment in enumerate(vbm.path.segments):
                if i in vbm.path.disconnected_before:
                    offset += gap_width
                    self.assertAlmostEqual(metadata['display_distances'][segment.start] -
                                           metadata['display_distances'][segment.start - 1], gap_width)
                np.testing.assert_allclose(metadata['display_distances'][segment] -
                                           physical_before[segment], offset, atol=1e-14)
                for j, band in enumerate(metadata['band_indices']):
                    line = ax.lines[i * len(metadata['band_indices']) + j]
                    self.assertEqual(line.get_linewidth(), 1.1)
                    np.testing.assert_array_equal(line.get_xdata(), metadata['display_distances'][segment])
                    np.testing.assert_allclose(np.diff(line.get_xdata()), np.diff(physical_before[segment]), atol=1e-14)
                    np.testing.assert_array_equal(line.get_ydata(), expected[segment, band])
            plt.close(fig)


if __name__ == '__main__':
    unittest.main()
