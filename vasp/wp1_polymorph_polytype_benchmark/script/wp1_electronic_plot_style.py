"""Local single-axes style and explicit, narrowly named WP1 figure exports."""
from __future__ import annotations

import argparse
from importlib.metadata import version
from pathlib import Path
import platform
import subprocess
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt

from wp1_electronic_plot_data import LABELS, PHASES, WP1

PALETTE = ("#0072B2", "#D55E00", "#7A5195", "#009E73", "#C23B70", "#7A8F00")
NEUTRAL = "#4D4D4D"
RC = {
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.weight": "bold", "font.size": 14, "mathtext.fontset": "stixsans",
    # Explicit regular math avoids a Matplotlib 3.11 STIX bold-glyph fallback loop.
    "mathtext.default": "regular",
    "axes.labelsize": 22, "axes.labelweight": "bold", "axes.linewidth": 1.8,
    "axes.edgecolor": "black", "axes.facecolor": "white", "axes.grid": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "xtick.labelsize": 14, "ytick.labelsize": 14,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.major.size": 4, "ytick.major.size": 4,
    "xtick.major.width": 1.8, "ytick.major.width": 1.8,
    "xtick.top": False, "ytick.right": False,
    "xtick.minor.visible": False, "ytick.minor.visible": False,
    "legend.fontsize": 10, "legend.framealpha": 1, "legend.edgecolor": "black",
    "legend.facecolor": "white", "lines.linewidth": 2,
    "pdf.fonttype": 42, "ps.fonttype": 42,
}
KINDS = ("tdos", "element_pdos", "selected_orbital_pdos", "band_structure")


def axes_labels(ax, phase, *, note=""):
    """Finish one axes without modifying the caller's global plotting defaults."""
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.8)
    ax.minorticks_off()
    ax.grid(False)
    label = r"$\alpha_1$" if phase == "alpha1" else LABELS[phase]
    ax.text(.025, .97, label, transform=ax.transAxes, va="top", fontsize=14, fontweight="bold")
    if note:
        ax.text(.025, .915, note, transform=ax.transAxes, va="top", fontsize=10)
    # ax.set_title("Optional descriptive title", fontweight="bold")


def energy_label(reference):
    return r"$E-E_F$ (eV)" if reference == "fermi" else r"$E-E_{VBM}$ (eV)"


def save_figure(fig, directory, phase, kind, *, output_format="pdf", overwrite=False):
    """Save exactly one single-panel file; replacement requires explicit consent.

    Only workflow-owned basenames are accepted. PDF dates are omitted to make
    reruns byte-reproducible in the same environment. PNG is an opt-in 600 dpi.
    """
    if phase not in PHASES or kind not in KINDS or output_format not in ("pdf", "svg", "png"):
        raise ValueError("Unknown phase, plot kind or single output format")
    if len(fig.axes) != 1:
        raise ValueError("Each WP1 export must have exactly one Axes")
    directory = Path(directory)
    path = directory / f"{phase}_pbe_{kind}.{output_format}"
    if path.exists() and not overwrite:
        raise FileExistsError(f"Use explicit overwrite to replace {path.name}")
    directory.mkdir(parents=True, exist_ok=True)
    metadata = {"CreationDate": None, "ModDate": None} if output_format == "pdf" else None
    with mpl.rc_context(RC):
        # Resolve constrained layout and legend geometry before the export backend
        # computes its tight bounding box (important with mixed dashed handles).
        fig.canvas.draw()
        fig.savefig(path, format=output_format, dpi=600, facecolor="white", transparent=False,
                    bbox_inches="tight", metadata=metadata)
    return path


def environment_summary():
    """Expose the actual executable, dependency versions and imported CMW checkout."""
    import cmw.analysis.electronic as electronic
    module = Path(electronic.__file__).resolve()
    result = {"Python": platform.python_version(), "executable": sys.executable,
              **{name: version(name) for name in ("numpy", "matplotlib", "pymatgen", "spglib", "nbformat", "nbclient", "ipykernel")},
              "CMW module": str(module)}
    checkout = next((p for p in module.parents if (p / ".git").exists()), None)
    if checkout is None:
        raise RuntimeError("Notebook provenance requires the documented editable CMW checkout")
    result["CMW commit"] = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip()
    return result


def display_figure(fig):
    """Display an in-memory notebook preview, without a second external export."""
    from io import BytesIO
    from IPython.display import Image, display
    buffer = BytesIO()
    with mpl.rc_context(RC):
        fig.savefig(buffer, format="png", dpi=120, facecolor="white", bbox_inches="tight")
    display(Image(data=buffer.getvalue()))


def cli_parser(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--phases", nargs="+", choices=PHASES, default=list(PHASES))
    parser.add_argument("--reference", choices=("fermi", "vbm"), default="fermi")
    parser.add_argument("--energy-limits", type=float, nargs=2)
    parser.add_argument("--output-format", choices=("pdf", "svg", "png"), default="pdf")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--overwrite", action="store_true", help="Replace only this workflow's explicitly generated names")
    return parser


def default_output(root):
    return root / WP1 / "results/04_electronic_structure/figures"
