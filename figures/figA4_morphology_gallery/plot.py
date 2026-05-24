#!/usr/bin/env python3
"""Figure B3 — solar-weighted volume-normalised backscatter vs effective
diameter for the silica_stardust monomer + 9 aggregate morphologies (R=250 nm).

Vertical layout:
  (a)  Scatter of beta_bar * sigma_sca / V vs d_eff, with a single dashed
       Mie reference line spanning d=300-900 nm.
  (b)  2 x 5 gallery of 3-D ball-and-stick thumbnails for the 10 morphologies,
       each stamped with the marker glyph that identifies it in panel (a).

beta_bar is defined by Eq. (XX) in the appendix body.  R_mono = 250 nm — stated
in the caption, no longer cluttering the scatter legend.

Run:
    PYTHONDONTWRITEBYTECODE=1 python3 plot.py
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import xarray as xr
import miepython as mp

HERE = Path(__file__).resolve().parent

STYLE_MORPH = {
    "monomer":                 ("monomer",          "o", "#009E73"),
    "dimer":                   ("dimer",            "o", "#0072B2"),
    "trimer_compact":          ("trimer compact",   "s", "#009E73"),
    "trimer_linear":           ("trimer linear",    "^", "#CC79A7"),
    "tetramer_compact":        ("tetramer cmp.",    "D", "#E69F00"),
    "tetramer_square":         ("tetramer sq.",     "p", "#9467BD"),
    "pentamer_semi_compact":   ("pent. compact",    "*", "#3F3F3F"),
    "pentamer_compact_planar": ("pent. planar",     "h", "#E377C2"),
    "pentamer_cross":          ("pent. cross",      "X", "#BCBD22"),
    "pentamer_square_pyramid": ("pent. sq. pyr.",   "P", "black"),
}
MARKER_GLYPH = {
    "o": "●", "s": "■", "^": "▲", "D": "◆", "p": "⬟",
    "*": "★", "h": "⬢", "X": "✕", "P": "✚",
}
CLUSTER_CENTERS = {
    "monomer": [(0,0,0)],
    "dimer":   [(-1,0,0), (1,0,0)],
    "trimer_compact": [
        (-1, -np.sqrt(3)/3, 0),
        ( 1, -np.sqrt(3)/3, 0),
        ( 0,  2*np.sqrt(3)/3, 0)],
    "trimer_linear": [(-2,0,0), (0,0,0), (2,0,0)],
    "tetramer_compact": [
        ( 1/np.sqrt(2),  1/np.sqrt(2),  1/np.sqrt(2)),
        ( 1/np.sqrt(2), -1/np.sqrt(2), -1/np.sqrt(2)),
        (-1/np.sqrt(2),  1/np.sqrt(2), -1/np.sqrt(2)),
        (-1/np.sqrt(2), -1/np.sqrt(2),  1/np.sqrt(2))],
    "tetramer_square": [(-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0)],
    "pentamer_semi_compact": [
        (-1, -np.sqrt(3)/3, 0),
        ( 1, -np.sqrt(3)/3, 0),
        ( 0,  2*np.sqrt(3)/3, 0),
        ( 0, 0,  2*np.sqrt(2/3)),
        ( 0, 0, -2*np.sqrt(2/3))],
    "pentamer_compact_planar": [
        (-2, 0, 0), (0, 0, 0), (2, 0, 0),
        (-1, np.sqrt(3), 0), (1, np.sqrt(3), 0)],
    "pentamer_cross": [
        (0,0,0), (2,0,0), (-2,0,0), (0,2,0), (0,-2,0)],
    "pentamer_square_pyramid": [
        (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0),
        (0, 0, np.sqrt(2))],
}
GALLERY_ORDER = [
    "monomer", "dimer",
    "trimer_compact", "trimer_linear",
    "tetramer_compact", "tetramer_square",
    "pentamer_semi_compact", "pentamer_compact_planar",
    "pentamer_cross", "pentamer_square_pyramid",
]

COLOR_MIE   = '#D55E00'
SPHERE_COLOR = '#3A7FB2'         # slightly cooler than silica blue to read as "geometry, not material"
BOND_COLOR  = '0.20'
SW_YLABEL = r'$\bar{\beta}\,\sigma_{\mathrm{sca}}/V$  [$\mu$m$^{-1}$]'

PENT_KEYS = ("pentamer_semi_compact", "pentamer_compact_planar",
             "pentamer_cross", "pentamer_square_pyramid")
TRI_KEYS  = ("trimer_compact", "trimer_linear")
TET_KEYS  = ("tetramer_compact", "tetramer_square")
def _xjitter(key, d):
    if key in PENT_KEYS:
        return d + 18.0 * (PENT_KEYS.index(key) - 1.5) / 1.5
    if key in TRI_KEYS:
        return d + 8.0 * (TRI_KEYS.index(key) - 0.5) / 0.5
    if key in TET_KEYS:
        return d + 8.0 * (TET_KEYS.index(key) - 1.0) / 1.0
    return d

plt.rcParams.update({
    'font.size': 10.0, 'axes.labelsize': 11.0, 'xtick.labelsize': 9.5,
    'ytick.labelsize': 9.5, 'legend.fontsize': 9.0, 'legend.frameon': True,
    'legend.framealpha': 0.92, 'axes.linewidth': 0.9, 'lines.linewidth': 1.8,
    'figure.dpi': 120,
})

# ---- nk loader + Mie reference -----------------------------------------
def load_yml_nk(path):
    lam, n, k = [], [], []
    with open(path) as f:
        in_data = False
        for line in f:
            if line.strip().startswith("data: |"):
                in_data = True; continue
            if not in_data: continue
            parts = line.strip().split()
            if len(parts) != 3: continue
            try:
                lam.append(float(parts[0])); n.append(float(parts[1])); k.append(float(parts[2]))
            except ValueError:
                pass
    return np.asarray(lam), np.asarray(n), np.asarray(k)

def make_nk_interp(yml_path):
    L, N, K = load_yml_nk(yml_path)
    def nk(wl_um):
        wl_um = np.atleast_1d(wl_um)
        return (np.interp(wl_um, L, N, left=N[0], right=N[-1]),
                np.interp(wl_um, L, K, left=K[0], right=0.0))
    return nk

def mie_sigb_over_V_solar(d_nm, wl_nm, I_sun_grid, nk_func):
    wl_um = np.asarray(wl_nm) / 1000.0
    n, k = nk_func(wl_um)
    theta = np.linspace(0.0, np.pi, 181); mu = np.cos(theta)
    R_um = (d_nm / 1000.0) / 2.0
    y = np.zeros_like(wl_um, dtype=float)
    for i, (wl_i, ni, ki) in enumerate(zip(wl_um, n, k)):
        m = complex(ni, -ki)
        x = 2 * np.pi * R_um / wl_i
        P = mp.i_unpolarized(m, x, mu, "4pi")
        beta = (1.0 / (2.0 * np.pi)) * np.trapezoid(P * np.sin(theta) * theta, theta)
        _, qsca, _, _ = mp.efficiencies(m, d_nm, wl_i * 1000.0)
        y[i] = (3.0 / 2.0) * (beta * qsca) / (d_nm / 1000.0)
    num = np.trapezoid(y * I_sun_grid, np.asarray(wl_nm))
    den = np.trapezoid(I_sun_grid, np.asarray(wl_nm))
    return float(num / den)

# ---- Scatter panel -----------------------------------------------------
def make_scatter_panel(ax, master_path, ref_yml_path):
    ds = xr.open_dataset(master_path)
    R_nm = float(ds.attrs["R_nm"])
    d_mono_nm = 2.0 * R_nm
    cfgs = [str(c) for c in ds.config.values]
    d_eff = ds.diameter_eff.values
    bq_cluster_solar = ds.bq_cluster_solar.values
    y_cluster = 1.5 * bq_cluster_solar / d_mono_nm * 1.0e3

    nk = make_nk_interp(ref_yml_path)
    sol_path = Path("/sessions/laughing-dreamy-fermi/mnt/y.lederer/cpdda/Solar_and_IR_spectrum.nc")
    wl_fine_nm = np.arange(300.0, 3001.0, 10.0)
    if sol_path.exists():
        sol = xr.open_dataset(sol_path)
        I_fine = sol.p_solar.interp(wavelength=wl_fine_nm/1000.0, method='linear',
                                    kwargs={'fill_value': 'extrapolate'}).values
        I_fine = np.maximum(I_fine, 0.0)
    else:
        # crude 6000 K Planck fallback (shouldn't be hit on the local machine)
        I_fine = np.ones_like(wl_fine_nm)

    d_sweep = np.arange(300.0, 901.0, 10.0)
    y_mie = np.array([mie_sigb_over_V_solar(d, wl_fine_nm, I_fine, nk) for d in d_sweep])
    ax.plot(d_sweep, y_mie, '--', color=COLOR_MIE, lw=2.0, zorder=2,
            label='Mie effective sphere')

    for i, cfg in enumerate(cfgs):
        if cfg not in STYLE_MORPH:
            continue
        label, marker, fc = STYLE_MORPH[cfg]
        ax.scatter(_xjitter(cfg, d_eff[i]), y_cluster[i],
                   marker=marker, s=130, facecolor=fc, edgecolor='black',
                   linewidth=0.9, zorder=3, label=label)

    ax.set_xlim(280.0, 900.0)
    y_all = np.concatenate([y_cluster, y_mie])
    pad = 0.04 * (y_all.max() - y_all.min())
    ax.set_ylim(y_all.min() - pad, y_all.max() + pad)
    ax.set_xlabel(r'$d_{\mathrm{eff}}$ [nm]')
    ax.set_ylabel(SW_YLABEL)
    ax.grid(alpha=0.3, zorder=0)

    # Compact legend — Mie line only.  R_mono noted in caption.
    from matplotlib.lines import Line2D
    mie_handle = Line2D([0], [0], ls='--', color=COLOR_MIE, lw=2.0,
                        label='Mie effective sphere')
    ax.legend(handles=[mie_handle], loc='upper right', fontsize=9.5,
              frameon=True, framealpha=0.92, handletextpad=0.5,
              borderpad=0.4)

# ---- Ball-and-stick thumbnail -----------------------------------------
UNIFORM_EXTENT = 3.0

def render_cluster(ax, centers, R=1.0, sphere_color=SPHERE_COLOR,
                   bond_color=BOND_COLOR, bond_lw=4.5, n_u=28, n_v=16):
    """Render a sphere cluster in a ball-and-stick style: thick gray bonds
    between touching monomers, then opaque shaded spheres on top.  Bonds make
    the cluster read as a 3-D structure instead of a flat coloured marker."""
    centers = np.asarray(centers, dtype=float)
    centroid = centers.mean(axis=0)

    # ----- bonds (drawn first, lower zorder so spheres overlay the bond ends)
    if len(centers) > 1:
        for i in range(len(centers)):
            for j in range(i + 1, len(centers)):
                d_ij = float(np.linalg.norm(centers[i] - centers[j]))
                if abs(d_ij - 2.0 * R) < 0.08 * R:        # touching pair
                    p, q = centers[i] - centroid, centers[j] - centroid
                    ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]],
                            color=bond_color, lw=bond_lw, alpha=0.9,
                            solid_capstyle='round', zorder=1)

    # ----- spheres
    u = np.linspace(0, 2*np.pi, n_u)
    v = np.linspace(0, np.pi, n_v)
    cu, su = np.cos(u), np.sin(u)
    sv, cv = np.sin(v), np.cos(v)
    for ci in centers:
        cxi, cyi, czi = ci - centroid
        x = R * np.outer(cu, sv) + cxi
        y = R * np.outer(su, sv) + cyi
        z = R * np.outer(np.ones_like(u), cv) + czi
        ax.plot_surface(x, y, z, color=sphere_color, alpha=1.0, linewidth=0,
                        antialiased=True, shade=True,
                        rcount=n_v, ccount=n_u, zorder=2)

    ax.set_xlim(-UNIFORM_EXTENT, UNIFORM_EXTENT)
    ax.set_ylim(-UNIFORM_EXTENT, UNIFORM_EXTENT)
    ax.set_zlim(-UNIFORM_EXTENT, UNIFORM_EXTENT)
    try:
        ax.set_box_aspect((1, 1, 1))
    except AttributeError:
        pass
    # Tilted view that breaks any "looks like a flat marker" reading: rotate
    # both elevation AND azimuth so even the planar pent-cross / square /
    # square-pyramid clusters show clear depth.
    ax.view_init(elev=22, azim=-55)
    ax.set_axis_off()
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor('white')

def make_gallery_panel(fig, gs_gallery, n_rows=2, n_cols=5):
    for i, key in enumerate(GALLERY_ORDER):
        r, c = divmod(i, n_cols)
        ax = fig.add_subplot(gs_gallery[r, c], projection='3d')
        render_cluster(ax, CLUSTER_CENTERS[key], R=1.0)
        label, marker, fc = STYLE_MORPH[key]
        glyph = MARKER_GLYPH.get(marker, marker)
        # Marker glyph in upper-left, coloured to match the scatter marker.
        # Drawn over a small white rounded box so the glyph reads clearly even
        # over a dark sphere.
        ax.text2D(0.04, 0.96, glyph, transform=ax.transAxes,
                  ha='left', va='top', fontsize=14, color=fc, fontweight='bold',
                  bbox=dict(facecolor='white', edgecolor='none',
                            boxstyle='round,pad=0.15', alpha=0.85))
        # Morphology label at the bottom
        ax.text2D(0.5, -0.04, label, transform=ax.transAxes,
                  ha='center', va='top', fontsize=9.5, color='black')


def main():
    """Absolute-positioned layout so panel (b)'s black frame aligns exactly
    with panel (a)'s DATA area (between the spines, not including the y-axis
    label/ticks).  Panel (b) is fixed at 3 cm tall."""
    PANEL_B_H_CM = 3.0
    PANEL_B_H_IN = PANEL_B_H_CM / 2.54     # = 1.181 in

    # Panel (a) DATA area (the rectangle inside the spines)
    A_DATA_W_IN = 5.7
    A_DATA_H_IN = 3.7

    # Margins (all in inches)
    L_M       = 0.85   # left  — fits y-axis label + ticks
    R_M       = 0.15
    TOP_M     = 0.15
    BOT_M     = 0.18
    XLABEL_M  = 0.50   # below panel (a) for the x-axis label + ticks
    SPACING   = 0.28   # gap between panel (a) bottom and panel (b) top

    FIG_W = L_M + A_DATA_W_IN + R_M
    FIG_H = TOP_M + A_DATA_H_IN + XLABEL_M + SPACING + PANEL_B_H_IN + BOT_M

    fig = plt.figure(figsize=(FIG_W, FIG_H))

    # ---- Panel (a) -----------------------------------------------------
    a_x = L_M / FIG_W
    a_y = (BOT_M + PANEL_B_H_IN + SPACING + XLABEL_M) / FIG_H
    a_w = A_DATA_W_IN / FIG_W
    a_h = A_DATA_H_IN / FIG_H
    ax_scatter = fig.add_axes([a_x, a_y, a_w, a_h])
    make_scatter_panel(
        ax_scatter,
        HERE / "data" / "R250_stardust" / "figures_appendix_B_master.nc",
        HERE / "materials" / "silica_stardust.yml",
    )

    # ---- Panel (b) black-framed box: width = panel (a) data width ------
    b_x = a_x
    b_w = a_w
    b_y_top = (BOT_M + PANEL_B_H_IN) / FIG_H
    b_h = PANEL_B_H_IN / FIG_H
    b_y = b_y_top - b_h

    rect = plt.Rectangle(
        (b_x, b_y), b_w, b_h,
        transform=fig.transFigure, fill=False,
        edgecolor="black", linewidth=1.0, zorder=0,
    )
    fig.patches.append(rect)

    # ---- 2x5 grid of thumbnails inside the box -------------------------
    n_rows, n_cols = 2, 5
    cell_w = b_w / n_cols
    cell_h = b_h / n_rows

    # Vertical split of each cell:  72 % sphere area, 28 % label area below.
    sphere_frac = 0.72
    label_frac  = 1.0 - sphere_frac

    # Small inset so 3-D axes don't touch the cell edges / box frame
    ix_pad = 0.015 * cell_w
    iy_pad = 0.015 * cell_h

    for i, key in enumerate(GALLERY_ORDER):
        r, c = divmod(i, n_cols)
        cx = b_x + c * cell_w
        cy_top = b_y_top - r * cell_h
        cy_bot = cy_top - cell_h

        ax_t = fig.add_axes(
            [cx + ix_pad,
             cy_bot + label_frac * cell_h,
             cell_w - 2 * ix_pad,
             sphere_frac * cell_h - iy_pad],
            projection="3d",
        )
        render_cluster(ax_t, CLUSTER_CENTERS[key])

        label, marker, fc = STYLE_MORPH[key]
        glyph = MARKER_GLYPH.get(marker, marker)
        ax_t.text2D(0.04, 0.96, glyph, transform=ax_t.transAxes,
                    ha="left", va="top", fontsize=9, color=fc,
                    fontweight="bold")

        # Morphology label inside the lower 28 % strip
        fig.text(cx + cell_w / 2,
                 cy_bot + label_frac * cell_h * 0.50,
                 label, ha="center", va="center",
                 fontsize=7.0, color="black")

    # ---- (a) and (b) panel letters, x-aligned --------------------------
    fig.text(a_x + 0.006, a_y + a_h - 0.005, "(a)",
             ha="left", va="top", fontsize=11, fontweight="bold")
    fig.text(b_x + 0.006, b_y_top - 0.005, "(b)",
             ha="left", va="top", fontsize=11, fontweight="bold")

    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"fig.{ext}", dpi=300, bbox_inches=None, pad_inches=0)
        print(f"saved fig.{ext}")
    plt.close(fig)


if __name__ == '__main__':
    main()
