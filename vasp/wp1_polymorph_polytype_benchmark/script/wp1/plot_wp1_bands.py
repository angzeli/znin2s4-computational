"""Plot actual WP1 folded PBE eigenvalues on the ordered calculation-cell path."""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from wp1_electronic_plot_data import energy_window, load_all, project_root
from wp1_electronic_plot_style import RC, PALETTE, cli_parser, default_output, energy_label, save_figure

DEFAULT_WINDOW = (-4., 4.)
PATH_FAMILY = {"spinel": "cF2", "alpha1": "hR1", "beta": "hP2", "IIa_prime": "hP2", "IIb": "hP2"}
BRANCH_GAPS = {"cF2": .0225, "hR1": .0375, "hP2": .055}  # fractions of physical path span
LABEL_Y = -.04
LABEL_SIZE = 12.5
LABEL_FONT = "Arial Narrow"  # one consistent face preserves the short c-direction spans


def display_label(label):
    return {"GAMMA": "Γ", "H_2": r"H$_2$", "H_0": r"H$_0$", "S_0": r"S$_0$", "S_2": r"S$_2$"}.get(label, label)


def display_distances(data, branch_gap=None):
    """Copy physical distances and offset whole branches for display only.

    The offsets have no reciprocal-space meaning and must not be used for
    derivatives or effective masses. Scientific CMW path data remain untouched.
    """
    if branch_gap is None:
        branch_gap = BRANCH_GAPS[PATH_FAMILY[data.phase]]
    if not np.isfinite(branch_gap) or branch_gap <= 0:
        raise ValueError("The display-gap fraction must be finite and positive")
    gap_width = branch_gap * np.ptp(data.path.distances)
    distances = data.path.distances.copy()
    offset = 0.
    for index, segment in enumerate(data.path.segments):
        if index in data.path.disconnected_before:
            offset += gap_width
        distances[segment] += offset
    return distances


def path_ticks(data, branch_gap=None):
    """Keep separate ordered endpoint occurrences around display-only gaps."""
    distances = display_distances(data, branch_gap)
    positions, labels = [], []
    for index, (segment, endpoints) in enumerate(zip(data.path.segments, data.endpoint_labels)):
        start, end = map(display_label, endpoints)
        if index == 0 or index in data.path.disconnected_before:
            positions.append(distances[segment.start])
            labels.append(start)
        positions.append(distances[segment.stop - 1])
        labels.append(end)
    return positions, labels


def plot_bands(data, *, window=DEFAULT_WINDOW, x_title="", branch_gap=None):
    """Draw prepared energies unchanged; never re-zero the sampled path maximum.

    Reference selection belongs to the shared loader. The publication default
    uses its audited dense-mesh VBM, also available to DOS with reference='vbm'.
    Folded-cell interpretation is returned as caption metadata, not drawn.
    """
    if x_title not in ("", "Wave vector"):
        raise ValueError("The optional x-axis title must be 'Wave vector' or empty")
    bands = data.band_energy
    selected = np.flatnonzero((bands.min(axis=0) <= window[1]) & (bands.max(axis=0) >= window[0]))
    if not len(selected):
        raise ValueError("No calculated bands intersect this window")
    if branch_gap is None:
        branch_gap = BRANCH_GAPS[PATH_FAMILY[data.phase]]
    distances = display_distances(data, branch_gap)
    ticks, labels = path_ticks(data, branch_gap)
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
        for segment in data.path.segments:
            for index in selected:
                ax.plot(distances[segment], bands[segment, index], color=PALETTE[0], linewidth=1.1)
        ax.axhline(0, color="0.55", linestyle=":", linewidth=.6, zorder=2.1)
        for position in ticks[1:-1]:
            ax.axvline(position, color="0.78", linewidth=.8, zorder=2.2)
        ax.set(xlim=(ticks[0], ticks[-1]), ylim=window, ylabel=energy_label(data.reference_mode),
               xlabel=x_title, xticks=ticks)
        ax.tick_params(axis="x", labelbottom=False)
        label_artists = [ax.text(x, LABEL_Y, label, transform=ax.get_xaxis_transform(),
                                ha="center", va="baseline", rotation=0, clip_on=False,
                                fontsize=LABEL_SIZE, fontfamily=LABEL_FONT, fontweight="bold")
                         for x, label in zip(ticks, labels)]
        ax.minorticks_off()
        ax.grid(False)
        # ax.set_title("Optional descriptive title", fontsize=18, fontweight="bold")
        fold = data.summary['primitive_cells']
    caption = (f"Calculated folded-cell representation ({fold} primitive cells); "
               "primitive-BZ direct/indirect character cannot be inferred solely from this plot."
               if fold > 1 else "Calculated-cell band representation.")
    return fig, ax, {"band_indices": selected, "energy": bands, "distances": data.path.distances,
                     "display_distances": distances, "display_branch_gap_fraction": branch_gap,
                     "label_artists": label_artists,
                     "segments": data.path.segments, "tick_labels": labels, "reference_eV": data.reference_ev,
                     "reference_mode": data.reference_mode, "primitive_cells": fold, "caption": caption}


def main():
    parser = cli_parser(__doc__)
    parser.set_defaults(reference="vbm", energy_limits=DEFAULT_WINDOW)
    parser.add_argument("--x-title", choices=("", "Wave vector"), default="")
    parser.add_argument("--branch-gap", type=float, default=None,
                        help="Override the family default: display-only fraction of physical path span per break")
    args = parser.parse_args()
    root = project_root(args.root)
    data = load_all(root, args.phases, args.reference)
    window = energy_window(data, "band", args.energy_limits)
    for phase in data:
        fig, _, _ = plot_bands(phase, window=window, x_title=args.x_title, branch_gap=args.branch_gap)
        try:
            print(save_figure(fig, args.output_dir or default_output(root), phase.phase, "band_structure",
                              output_format=args.output_format, overwrite=args.overwrite))
        finally:
            plt.close(fig)


if __name__ == "__main__":
    main()
