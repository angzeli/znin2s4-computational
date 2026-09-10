"""Standalone WP1 DOS/PDOS figures from validated, unbroadened sampled arrays.

Run this module explicitly to export; importing it creates no figure or file.
Use --overwrite on reruns to replace only its deterministic output names.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from wp1_electronic_plot_data import SELECTED, energy_window, load_all, prepared_dos, project_root
from wp1_electronic_plot_style import (NEUTRAL, PALETTE, RC, cli_parser,
                                      default_output, energy_label, save_figure)

COLORS = dict(zip(("Zn", "In", "S"), PALETTE)) | dict(zip(SELECTED, PALETTE)) | {"TDOS": NEUTRAL}
FAMILIES = ("selected_orbital_pdos",)
DEFAULT_WINDOW = (-3., 3.)
ELEMENT_WINDOW = (-8., 4.)


def families_for_mode(mode):
    """Publication is selected PDOS only; element PDOS requires explicit opt-in."""
    if mode == "publication":
        return FAMILIES
    if mode == "element-pdos":
        return ("element_pdos",)
    raise ValueError("Choose publication or element-pdos")


def channels_for(kind, selected=SELECTED):
    if kind == "tdos":
        return ("TDOS",)
    if kind == "element_pdos":
        return ("TDOS", "Zn", "In", "S")
    if kind != "selected_orbital_pdos" or not selected or len(set(selected)) != len(selected):
        raise ValueError("Choose a known DOS family and distinct selected channels")
    return tuple(selected)


def adaptive_y_limits(data, kind, window, *, normalization="formula_unit", selected=SELECTED):
    """Physical per-phase scale: 7–8% headroom over displayed sampled curves."""
    energy, curves = prepared_dos(data, normalization)
    inside = (energy >= window[0]) & (energy <= window[1])
    if not inside.any():
        raise ValueError("No sampled DOS values in the requested window")
    channels = channels_for(kind, selected)
    if any(channel not in curves for channel in channels):
        raise ValueError("Unavailable projection")
    maximum = max(float(np.max(curves[channel][inside])) for channel in channels)
    if maximum <= 0 or not np.isfinite(maximum):
        raise ValueError("Empty or nonfinite DOS figure family")
    quantum = 10. ** (np.floor(np.log10(maximum)) - 2)
    return (0., float(np.ceil(1.07 * maximum / quantum) * quantum))


def common_y_limits(data, window, *, normalization="formula_unit", selected=SELECTED,
                    families=("element_pdos", "selected_orbital_pdos")):
    """Explicit diagnostic common scale; never the standalone default."""
    return {kind: (0., max(adaptive_y_limits(phase, kind, window,
                normalization=normalization, selected=selected)[1] for phase in data))
            for kind in families}


def plot_dos(data, kind, *, window=DEFAULT_WINDOW, ylim=None, normalization="formula_unit", selected=SELECTED):
    """Return one axes and the exact plotted arrays; no smoothing or file writes."""
    energy, curves = prepared_dos(data, normalization)
    channels = channels_for(kind, selected)
    if ylim is None:
        ylim = adaptive_y_limits(data, kind, window, normalization=normalization, selected=selected)
    # Full channel identity determines colours, not the selected subset's order.
    extras = {name: PALETTE[i % len(PALETTE)] for i, name in enumerate(sorted(set(curves) - set(COLORS)))}
    colors = extras | COLORS
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(6.5, 5), layout="constrained")
        for channel in channels:
            if channel not in curves:
                plt.close(fig)
                raise ValueError(f"Unavailable projection: {channel}")
            ax.plot(energy, curves[channel], color=colors[channel], label=channel,
                    linestyle="--" if channel == "TDOS" and kind == "element_pdos" else "-", linewidth=2)
        ax.axvline(0, color="0.65", linestyle=":", linewidth=.9, zorder=0)
        ax.set(xlim=window, ylim=ylim, xlabel=energy_label(data.reference_mode))
        units = r"states eV$^{-1}$ f.u.$^{-1}$" if normalization == "formula_unit" else r"states eV$^{-1}$ cell$^{-1}$"
        ax.set_ylabel(f"PDOS ({units})")
        ax.minorticks_off()
        ax.grid(False)
        ax.legend(loc="best", prop={"weight": "bold", "size": 10})
    return fig, ax, {"energy": energy, "curves": {k: curves[k] for k in channels},
                     "reference_eV": data.reference_ev, "normalization": normalization}


def main():
    parser = cli_parser(__doc__)
    parser.set_defaults(reference="vbm")
    parser.add_argument("--mode", choices=("publication", "element-pdos"), default="publication")
    parser.add_argument("--common-y", action="store_true", help="Explicit across-phase diagnostic y scale")
    parser.add_argument("--normalization", choices=("formula_unit", "cell"), default="formula_unit")
    parser.add_argument("--selected", nargs="+", default=list(SELECTED))
    args = parser.parse_args()
    root = project_root(args.root)
    data = load_all(root, args.phases, args.reference)
    families = families_for_mode(args.mode)
    requested_window = args.energy_limits or (DEFAULT_WINDOW if args.mode == "publication" else ELEMENT_WINDOW)
    window = energy_window(data, "dos", requested_window)
    limits = (common_y_limits(data, window, normalization=args.normalization,
                             selected=args.selected, families=families) if args.common_y else {})
    for phase in data:
        for kind in families:
            fig, _, _ = plot_dos(phase, kind, window=window, ylim=limits.get(kind), normalization=args.normalization, selected=args.selected)
            try:
                print(save_figure(fig, args.output_dir or default_output(root), phase.phase, kind,
                                  output_format=args.output_format, overwrite=args.overwrite))
            finally:
                plt.close(fig)


if __name__ == "__main__":
    main()
