"""Ideal powder XRD from the five accepted PBE+D3(BJ) Stage 01 CONTCARs.

Provenance: results/01_geometry_optimisation/FINAL_GEOMETRY_OPT_AUDIT.md
and results/03_static_scf/STATIC_SCF_ANALYSIS.md. The latter confirms that
Stage 03 retained these geometries; its POSCARs are checked, never substituted.
No cell reduction, refinement or symmetrization is applied to diffraction.

CuKa is pymatgen's single effective wavelength convention, not an explicit
K-alpha doublet or a confirmed experimental source. Peaks include pymatgen's
atomic scattering and Lorentz-polarization factors, with no thermal factors.
One CSV row is a merged powder peak, not necessarily one unique hkl family:
semicolon-separated hkl and multiplicity entries correspond in order. Hexagonal
CSV indices are Miller-Bravais (h k i l); figure labels omit redundant i only.
Library peak merging and weak-peak thresholds are recorded in the summary.

Thin-film comparisons should emphasize positions/fingerprints, not quantitative
powder intensities: texture, strain, finite size, substrate, mixed phases and
instrument response can change observed patterns. Strong (00l) texture is not
powder phase abundance. DFT lattice error, temperature and strain can shift 2theta.
Sticks are unbroadened, individually normalized within the requested range, and
vertically offset for display only. No size/strain inference or experimental fit.

Only results/simulated_xrd is written (apart from temporary library caches).
Import and --help perform no analysis. Use --overwrite to regenerate owned files.
"""
from __future__ import annotations

import argparse
import csv
from importlib.metadata import version
from io import BytesIO, StringIO
import math
import os
from pathlib import Path
import tempfile

WP1 = Path(__file__).resolve().parents[1]
PHASES = ("spinel", "alpha1", "beta", "IIa_prime", "IIb")
LABELS = ("Spinel", "α₁", "β", "IIa′", "IIb")
COLOURS = ("#0072B2", "#D55E00", "#7A5195", "#009E73", "#C23B70")
# Audited full calculation cells, not conventional/primitive replacements.
AUDITED = {
    "spinel": ("spinel", 8, 10.617497, 10.617497, 90.),
    "alpha1": ("alpha1", 3, 3.883399, 36.494483, 120.),
    "beta": ("beta", 1, 3.878284, 12.148433, 120.),
    "IIa_prime": ("IIa_to_IIa_prime", 2, 3.889649, 24.560599, 120.),
    "IIb": ("IIb", 2, 3.883347, 24.351774, 120.),
}


def parser() -> argparse.ArgumentParser:
    """Build the dependency-free CLI; wavelength overrides the named radiation."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--radiation", default="CuKa", help="pymatgen radiation name (default: CuKa)")
    p.add_argument("--wavelength", type=float, help="explicit wavelength in angstrom; overrides --radiation")
    p.add_argument("--two-theta-min", type=float, default=5., help="lower 2theta limit in degrees")
    p.add_argument("--two-theta-max", type=float, default=80., help="upper 2theta limit in degrees")
    p.add_argument("--format", choices=("pdf", "png"), default="pdf", help="one comparison format; PDF is canonical")
    p.add_argument("--overwrite", action="store_true", help="replace only this analysis's named outputs")
    return p


def accepted_structure(phase: str):
    """Read the exact Stage 01 cell and verify its audited Stage 03 handoff."""
    import numpy as np
    from pymatgen.io.vasp import Poscar

    folder, units, a, c, gamma = AUDITED[phase]
    path = WP1 / "calculation/01_geometry_optimisation" / folder / "CONTCAR"
    structure = Poscar.from_file(path, check_for_potcar=False).structure
    reference = Poscar.from_file(WP1 / "calculation/03_static_scf" / phase / "POSCAR",
                                 check_for_potcar=False).structure
    expected = {"Zn": units, "In": 2 * units, "S": 4 * units}
    if not structure.is_ordered or structure.composition.get_el_amt_dict() != expected:
        raise ValueError(f"{phase}: composition/occupancy differs from accepted ZnIn2S4 cell")
    if not np.allclose(structure.lattice.parameters, (a, a, c, 90., 90., gamma), atol=2e-6, rtol=0):
        raise ValueError(f"{phase}: cell differs from rounded accepted geometry audit")
    if (structure.species != reference.species
            or not np.allclose(structure.lattice.matrix, reference.lattice.matrix, atol=1e-8, rtol=0)
            or not np.allclose(structure.frac_coords, reference.frac_coords, atol=1e-8, rtol=0)):
        raise ValueError(f"{phase}: Stage 01 CONTCAR no longer matches accepted Stage 03 POSCAR")
    return path, structure


def validate_pattern(pattern, limits, wavelength: float) -> None:
    """Enforce finite, ordered, normalized stick data and Bragg-law consistency."""
    import numpy as np

    x, y, d = map(np.asarray, (pattern.x, pattern.y, pattern.d_hkls))
    if (not len(x) or len(x) != len(pattern.hkls) or x.shape != y.shape or x.shape != d.shape
            or not np.isfinite([x, y, d]).all() or not (np.diff(x) > 0).all()
            or not ((x >= limits[0]) & (x <= limits[1])).all()
            or not (y >= 0).all() or not (d > 0).all()
            or not np.isclose(y.max(), 100., atol=1e-9, rtol=0)):
        raise ValueError("Invalid or empty normalized diffraction pattern")
    if not np.allclose(2 * d * np.sin(np.deg2rad(x / 2)), wavelength, atol=1e-8, rtol=0):
        raise ValueError("Reflection positions and d-spacings violate Bragg's law")


def family_text(families, compact=False) -> str:
    """Retain every contributing family; drop redundant hexagonal i only in labels."""
    labels = []
    for family in families:
        hkl = tuple(family["hkl"])
        if compact and len(hkl) == 4:
            hkl = (hkl[0], hkl[1], hkl[3])
        labels.append("(" + " ".join(map(str, hkl)) + ")")
    return "; ".join(labels)


def family_plot_label(families) -> str:
    """Format display-only three-index labels with overbars on negative indices."""
    labels = []
    for family in families:
        hkl = tuple(family["hkl"])
        if len(hkl) == 4:
            hkl = (hkl[0], hkl[1], hkl[3])
        if any(index < 0 for index in hkl):
            indices = [rf"\overline{{{abs(index)}}}" if index < 0 else str(index) for index in hkl]
            labels.append(r"$\mathbf{(" + r"\ ".join(indices) + r")}$")
        else:
            labels.append("(" + " ".join(map(str, hkl)) + ")")
    return "; ".join(labels)


def simulate(radiation, limits):
    """Calculate sticks on unchanged cells; symmetry is metadata-only."""
    import numpy as np
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    calculator = XRDCalculator(wavelength=radiation, symprec=0)
    results = []
    for phase, label in zip(PHASES, LABELS):
        path, structure = accepted_structure(phase)
        before = structure.as_dict()
        pattern = calculator.get_pattern(structure, scaled=True, two_theta_range=limits)
        validate_pattern(pattern, limits, calculator.wavelength)
        analyzers = [SpacegroupAnalyzer(structure, symprec=t) for t in (1e-4, 1e-3, 1e-2)]
        groups = [(analyzer.get_space_group_symbol(), analyzer.get_space_group_number())
                  for analyzer in analyzers]
        if structure.as_dict() != before:
            raise ValueError(f"{phase}: analysis unexpectedly changed source structure")
        robust = all(group == groups[0] for group in groups)
        strongest = int(np.argmax(pattern.y))
        summary = dict(phase=phase, display_name=label,
                       source_structure=path.relative_to(WP1).as_posix(), structure_file_type="VASP CONTCAR",
                       composition=structure.composition.formula, reduced_composition="ZnIn2S4",
                       atom_count=len(structure), **dict(zip(
                           ("a_angstrom", "b_angstrom", "c_angstrom", "alpha_deg", "beta_deg", "gamma_deg"),
                           structure.lattice.parameters)),
                       space_group=groups[0][0] if robust else "tolerance-dependent",
                       space_group_number=groups[0][1] if robust else "",
                       symmetry_tolerances_angstrom="0.0001;0.001;0.01", xrd_symprec=0,
                       radiation=radiation if isinstance(radiation, str) else "explicit wavelength",
                       wavelength_angstrom=calculator.wavelength,
                       two_theta_min_deg=limits[0], two_theta_max_deg=limits[1],
                       number_of_reflections=len(pattern.x), strongest_peak_2theta_deg=pattern.x[strongest],
                       strongest_peak_hkl=family_text(pattern.hkls[strongest]),
                       strongest_peak_d_angstrom=pattern.d_hkls[strongest],
                       hkl_convention="hkil" if structure.lattice.is_hexagonal() else "hkl",
                       two_theta_merge_tolerance_deg=calculator.TWO_THETA_TOL,
                       relative_intensity_cutoff=calculator.SCALED_INTENSITY_TOL,
                       debye_waller_factors="none", intensity_normalization="max=100 within requested range",
                       pymatgen_version=version("pymatgen"), numpy_version=version("numpy"),
                       matplotlib_version=version("matplotlib"), spglib_version=version("spglib"))
        rows = [dict({"2theta_deg": x}, relative_intensity=y, d_spacing_angstrom=d,
                     hkl=family_text(families),
                     multiplicity=";".join(str(f["multiplicity"]) for f in families),
                     total_multiplicity=sum(f["multiplicity"] for f in families),
                     wavelength_angstrom=calculator.wavelength)
                for x, y, d, families in zip(pattern.x, pattern.y, pattern.d_hkls, pattern.hkls)]
        results.append((summary, pattern, rows))
    return results


def comparison_figure(results, limits, output_format: str) -> bytes:
    """Render unchanged sticks with common vertical offsets and three major labels."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np

    rc = {
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "mathtext.fontset": "stixsans", "text.color": "black", "axes.edgecolor": "black",
        "axes.labelcolor": "black", "axes.labelsize": 22, "axes.labelweight": "bold",
        "axes.titlesize": 18, "axes.titleweight": "bold", "axes.linewidth": 1.8, "axes.grid": False,
        "xtick.labelsize": 14, "ytick.labelsize": 14, "xtick.color": "black", "ytick.color": "black",
        "xtick.direction": "in", "ytick.direction": "in", "xtick.major.width": 1.8, "ytick.major.width": 1.8,
        "xtick.major.size": 4, "ytick.major.size": 4, "xtick.bottom": True, "ytick.left": True,
        "xtick.top": False, "ytick.right": False, "xtick.minor.visible": False, "ytick.minor.visible": False,
        "lines.linewidth": 2., "legend.fontsize": 10, "legend.frameon": True, "legend.framealpha": 1.,
        "legend.facecolor": "white", "legend.edgecolor": "black", "savefig.dpi": 600,
        "savefig.bbox": "tight", "savefig.edgecolor": "white", "savefig.transparent": False, "pdf.fonttype": 42,
    }
    with mpl.rc_context(rc):
        fig, ax = plt.subplots(figsize=(8, 6))
        for index, ((summary, pattern, _), colour) in enumerate(zip(results, COLOURS)):
            base = (len(results) - 1 - index) * 150
            ax.hlines(base, *limits, color="0.8", linewidth=.6)
            ax.vlines(pattern.x, base, base + pattern.y, color=colour, linewidth=2.)
            label = r"α$_{\mathbf{1}}$" if summary["phase"] == "alpha1" else summary["display_name"]
            ax.text(.98, base + 85, label, transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=14, fontweight="bold", color=colour)
            # Labels never move reflection positions. Reserve lateral space by
            # selecting at most three strong peaks separated by 4.5 degrees.
            selected = []
            for i in np.argsort(-pattern.y, kind="stable"):
                if all(abs(pattern.x[i] - pattern.x[j]) >= 4.5 for j in selected):
                    selected.append(i)
                if len(selected) == 3:
                    break
            for i in selected:
                ax.text(pattern.x[i], base + pattern.y[i] + 6,
                        family_plot_label(pattern.hkls[i]), ha="center", va="bottom",
                        fontsize=10, fontweight="bold")
        ax.set(xlim=limits, ylim=(-15, 755), xlabel="2θ (°)", ylabel="Relative intensity", yticks=[])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.8)
        for tick in ax.get_xticklabels():
            tick.set_fontweight("bold")
        ax.grid(False, which="both")
        ax.minorticks_off()
        # ax.set_title("WP1 simulated powder XRD", fontsize=18, fontweight="bold")
        fig.tight_layout()
        buffer = BytesIO()
        metadata = {"CreationDate": None, "ModDate": None} if output_format == "pdf" else None
        fig.savefig(buffer, format=output_format, metadata=metadata, bbox_inches="tight",
                    facecolor="white", transparent=False, dpi=600)
        plt.close(fig)
    return buffer.getvalue()


def csv_bytes(rows) -> bytes:
    """Serialize UTF-8 CSV deterministically without an index or metadata preamble."""
    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def main(argv=None) -> None:
    """Validate and calculate in memory, then write only the seven owned outputs."""
    p = parser()
    args = p.parse_args(argv)
    limits = (args.two_theta_min, args.two_theta_max)
    if not all(map(math.isfinite, limits)) or not 0 < limits[0] < limits[1] < 180:
        p.error("Require 0 < two-theta-min < two-theta-max < 180 degrees")
    if args.wavelength is not None and (not math.isfinite(args.wavelength) or args.wavelength <= 0):
        p.error("Wavelength must be finite and positive")
    directory = WP1 / "results/simulated_xrd"
    names = [f"{phase}_xrd_reflections.csv" for phase in PHASES]
    names += ["wp1_simulated_xrd_summary.csv", f"wp1_simulated_xrd_comparison.{args.format}"]
    if directory.is_symlink() or any((directory / name).is_symlink() for name in names):
        p.error("Refusing symlinked output paths")
    if not args.overwrite and any((directory / name).exists() for name in names):
        p.error("Outputs exist; use --overwrite to regenerate only the named XRD outputs")
    # Keep Matplotlib/font caches away from scientific directories and active runtimes.
    with tempfile.TemporaryDirectory(prefix="wp1-xrd-") as cache:
        keys = {"MPLCONFIGDIR": cache, "XDG_CACHE_HOME": cache, "MPLBACKEND": "Agg"}
        previous = {key: os.environ.get(key) for key in keys}
        os.environ.update(keys)
        try:
            results = simulate(args.wavelength if args.wavelength is not None else args.radiation, limits)
            contents = [csv_bytes(rows) for _, _, rows in results]
            contents += [csv_bytes([summary for summary, _, _ in results]),
                         comparison_figure(results, limits, args.format)]
            directory.mkdir(parents=True, exist_ok=True)
            for name, content in zip(names, contents):
                (directory / name).write_bytes(content)
            for summary, _, _ in results:
                print(f"{summary['phase']}: {summary['number_of_reflections']} peaks; "
                      f"maximum at {summary['strongest_peak_2theta_deg']:.4f} degrees")
            print(f"Wrote {len(names)} files to {directory}")
        finally:
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


if __name__ == "__main__":
    main()
