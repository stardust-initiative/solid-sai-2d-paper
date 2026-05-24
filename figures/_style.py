"""Shared style helper for manuscript_figures/*/plot.py — ACP-ready defaults.

Single colour-source-of-truth: every diverging or sequential colormap used in
the manuscript is the master DIV_BWR (blue → white → yellow → orange → red,
modelled on ACP article 23/13665/2023 Fig 2 panel (a)) or a derived half of
it. Halves are derived programmatically so colour stops stay in lock-step
with the master.

Use:
  DIV_BWR        — diverging (data crosses a reference, e.g. ratio crossing 1)
  WARM_SEQ       — sequential warm with white at zero (e.g. mass mixing ratio)
  WARM_NOWHITE   — sequential warm without white anchor (e.g. ratio > 1)
  COOL_SEQ       — sequential cool with NAVY at vmin and white at vmax
                    (so "more cooling = darker", matching WARM_SEQ convention)
  COOL_NOWHITE   — sequential cool without white anchor

Usage at the top of each plot.py:

    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from _style import (apply, DIV_BWR, WARM_SEQ, WARM_NOWHITE, COOL_SEQ, COOL_NOWHITE,
                        centered_norm, panel_label, savefig_both)
    apply()
"""
import os
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm


# =============================================================================
# Master diverging colormap — single source of truth.
# Cool half (white at center → pale blue → light blue → mid blue → navy):
#   ColorBrewer Blues lightness ramp, mirrored against the warm half.
# Warm half (white at center → pale yellow → yellow → orange → red):
#   ColorBrewer YlOrRd palette, matching the warm half of NCL/ACP house style.
# Pair with TwoSlopeNorm(vmin=neg, vcenter=ref, vmax=pos) for diverging plots,
# or sample the appropriate half (see derived sub-maps below) for sequential.
# =============================================================================
DIV_BWR = LinearSegmentedColormap.from_list("div_bwr", [
    (0.00, "#08306B"),  # darkest navy
    (0.07, "#08519C"),  # very deep blue
    (0.14, "#2171B5"),
    (0.22, "#4292C6"),
    (0.29, "#6BAED6"),
    (0.37, "#9ECAE1"),
    (0.44, "#C6DBEF"),
    (0.49, "#F7FBFF"),  # almost-white pale blue
    (0.50, "#FFFFFF"),  # exact white at vcenter
    (0.51, "#FFF5BB"),  # very pale yellow — jumps off white quickly
    (0.56, "#FFEDA0"),
    (0.63, "#FED976"),
    (0.71, "#FEB24C"),
    (0.79, "#FD8D3C"),
    (0.86, "#FC4E2A"),
    (0.93, "#E31A1C"),
    (1.00, "#BD0026"),  # darkest red
], N=256)


def _derive_half(p_start, p_end, name):
    """Sample 256 colors along DIV_BWR[p_start:p_end] as a ListedColormap."""
    return ListedColormap(DIV_BWR(np.linspace(p_start, p_end, 256)), name=name)


# Sequential warm with WHITE at zero — for ranges where 0 is the meaningful
# floor (mass mixing ratio, lifetime >= 0). cmap(0)=white, cmap(1)=darkest red.
# Pair with Normalize(vmin=0, vmax=peak) or LogNorm; set_under('white') for
# below-vmin cells.
WARM_SEQ = _derive_half(0.50, 1.00, "warm_seq")
WARM_SEQ.set_under("white")

# Sequential warm WITHOUT a white anchor — for ranges where the lower bound is
# a meaningful nonzero value (monomer diameter, ratio just above 1.0). Lowest
# cell reads pale yellow (not white) so it stays distinct from NaN cut-out cells.
WARM_NOWHITE = _derive_half(0.51, 1.00, "warm_nowhite")

# Sequential cool with NAVY at the most-negative end (cmap(0)) and WHITE at the
# zero/no-anomaly end (cmap(1)). Convention: "more cooling = darker" matches
# WARM_SEQ's "more positive = darker red". Pair with
# Normalize(vmin=most_negative, vmax=0_or_slightly_negative). set_over('white')
# caps any cells above vmax (e.g., positive warming cells) to off-scale white.
COOL_SEQ = _derive_half(0.00, 0.50, "cool_seq")
COOL_SEQ.set_over("white")

# Sequential cool WITHOUT a white anchor at the high end — for ranges where the
# upper bound is a meaningful nonzero value (e.g., ratio just below 1.0).
COOL_NOWHITE = _derive_half(0.00, 0.49, "cool_nowhite")


def apply():
    """Install ACP-compliant matplotlib defaults. Call once at the top of each plot.py."""
    mpl.rcParams.update({
        "font.family":     "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size":       10,
        "axes.labelsize":  10,
        "axes.titlesize":  10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "pdf.fonttype":    42,
        "ps.fonttype":     42,
        "savefig.dpi":     300,
        "savefig.bbox":    "tight",
        "axes.unicode_minus": True,
    })


def centered_norm(vmin, vmax, vcenter=0.0):
    """TwoSlopeNorm for diverging cmaps centered on an arbitrary reference."""
    return TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)


def panel_label(ax, letter, x=0.02, y=0.97, **text_kw):
    """ACP convention: '(a)' lower-case, bold, top-left, on a white bbox."""
    defaults = dict(fontsize=12, fontweight="bold")
    defaults.update(text_kw)
    ax.text(x, y, f"({letter})", transform=ax.transAxes,
            ha="left", va="top",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=1.5),
            **defaults)


def savefig_both(fig, here, basename="fig"):
    """Write basename.pdf (vector, submission) AND basename.png (raster, preview)."""
    pdf_path = os.path.join(here, f"{basename}.pdf")
    png_path = os.path.join(here, f"{basename}.png")
    fig.savefig(pdf_path)
    fig.savefig(png_path)
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")
