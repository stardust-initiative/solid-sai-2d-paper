#!/usr/bin/env python3
"""Figure B1 — silica & calcite material comparison at monomer d = 500 nm.

Two-panel layout (1 row x 2 cols):
  (a)  beta_bar * sigma_sca / V over the solar range, linear y       [um^-1]
  (b)  sigma_abs / V over the thermal range, log y                   [um^-1]

The beta_bar definition is given by Eq. (XX) in the appendix body.

Usage:
    PYTHONDONTWRITEBYTECODE=1 python3 plot.py
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

COLORS = {'silica': '#0072B2', 'calcite': '#E69F00'}
MATS = ('silica', 'calcite')
PANEL_A_XLIM = (0.0, 3.0)
PANEL_B_XLIM = (0.0, 50.0)
PANEL_B_YMIN = 1e-3
WEIGHT_PEAK_FRAC = 0.90
WEIGHT_FILL = dict(color='0.6', alpha=0.14)
WEIGHT_LINE = dict(color='k', lw=0.9, alpha=0.5)
SW_YLABEL = r'$\bar{\beta}\,\sigma_{\mathrm{sca}}/V$  [$\mu$m$^{-1}$]'
LW_YLABEL = r'$\sigma_{\mathrm{abs}}/V$  [$\mu$m$^{-1}$]'

plt.rcParams.update({
    'font.size': 10.0, 'axes.labelsize': 11.0, 'xtick.labelsize': 9.5,
    'ytick.labelsize': 9.5, 'legend.fontsize': 10.0, 'legend.frameon': True,
    'legend.framealpha': 0.92, 'axes.linewidth': 0.9, 'lines.linewidth': 1.8,
    'figure.dpi': 120,
})

_H, _C, _KB = 6.62607015e-34, 2.99792458e8, 1.380649e-23
def planck_lambda(lam_um, T_K):
    lam_m = np.asarray(lam_um, float) * 1e-6
    u = _H * _C / (lam_m * _KB * T_K)
    out = np.zeros_like(lam_m); ok = (u > 1e-8) & (u < 700.0)
    out[ok] = 1.0 / lam_m[ok] ** 5 / np.expm1(u[ok])
    return out

def _trapz(y, x):
    return float(np.trapezoid(y, x)) if hasattr(np, 'trapezoid') else float(np.trapz(y, x))

def sw_weight(lam_um, T_sw, cutoff_um):
    w = planck_lambda(lam_um, T_sw)
    return np.where(np.asarray(lam_um) >= cutoff_um, w, 0.0)

def lw_windows(bounds_cm, T_list, lam_min, lam_max):
    nu_edges = np.concatenate([[0.0], np.asarray(bounds_cm, float), [np.inf]])
    bands = []
    for k in range(len(T_list)):
        nu_lo, nu_hi = nu_edges[k], nu_edges[k + 1]
        lam_hi = (1e4 / nu_lo) if nu_lo > 0 else lam_max
        lam_lo = (1e4 / nu_hi) if np.isfinite(nu_hi) else lam_min
        lam_hi, lam_lo = min(lam_hi, lam_max), max(lam_lo, lam_min)
        bands.append((lam_lo, lam_hi, float(T_list[k])) if lam_hi > lam_lo else None)
    flux = []
    for b in bands:
        if b is None:
            flux.append(0.0); continue
        xx = np.linspace(b[0], b[1], 2000); flux.append(_trapz(planck_lambda(xx, b[2]), xx))
    flux = np.asarray(flux); wk = flux / flux.sum()
    return [(b[0], b[1], b[2], wk[k]) for k, b in enumerate(bands) if b is not None]

def lw_weight_curve(lam_um, windows):
    lam = np.asarray(lam_um, float); w = np.zeros_like(lam)
    for lo, hi, T_k, wk in windows:
        m = (lam >= lo) & (lam <= hi)
        if m.sum() < 2:
            continue
        Bk = planck_lambda(lam[m], T_k); nrm = _trapz(Bk, lam[m])
        if nrm > 0:
            w[m] = wk * Bk / nrm
    return w

def overlay_weight_linear(ax, lam, w, xlim, peak_frac=WEIGHT_PEAK_FRAC):
    sel = (lam >= xlim[0]) & (lam <= xlim[1])
    y0, y1 = ax.get_ylim(); wmax = float(np.nanmax(w[sel])) or 1.0
    ws = y0 + (w / wmax) * peak_frac * (y1 - y0)
    ax.fill_between(lam[sel], y0, ws[sel], zorder=0, **WEIGHT_FILL)
    ax.plot(lam[sel], ws[sel], zorder=0.5, **WEIGHT_LINE)

def overlay_weight_logpanel(ax, lam, w, xlim, peak_frac=WEIGHT_PEAK_FRAC):
    sel = (lam >= xlim[0]) & (lam <= xlim[1])
    wmax = float(np.nanmax(w[sel])) or 1.0
    ymin, ymax = ax.get_ylim()
    w_norm = np.clip(w / wmax, 1e-6, 1.0)
    log_min = np.log10(ymin); log_top = np.log10(peak_frac * ymax)
    y_overlay = 10.0 ** (log_min + (log_top - log_min) * w_norm)
    ax.fill_between(lam[sel], ymin, y_overlay[sel], zorder=0, **WEIGHT_FILL)
    ax.plot(lam[sel], y_overlay[sel], zorder=0.5, **WEIGHT_LINE)


def make_fig(D, out_base, panel_b_logy):
    a_um = float(D['r_um'])
    lam_sw, lam_lw = D['lam_sw_um'], D['lam_lw_um']
    fac = 0.75 / a_um
    T_sw, cut = float(D['T_sw_K']), float(D['sw_cutoff_um'])
    windows = lw_windows(D['lw_bounds_cm'], D['T_lw_K'], float(lam_lw[0]), float(lam_lw[-1]))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.6, 3.6), constrained_layout=True)

    ymax_a = 0.0
    for m in MATS:
        y = fac * D[f'{m}_betabar_sw'] * D[f'{m}_qsca_sw']
        sel = np.isfinite(y)
        axA.plot(lam_sw[sel], y[sel], color=COLORS[m], zorder=3, label=m)
        in_x = sel & (lam_sw >= PANEL_A_XLIM[0]) & (lam_sw <= PANEL_A_XLIM[1])
        ymax_a = max(ymax_a, float(np.nanmax(y[in_x])))
    axA.set_xlim(*PANEL_A_XLIM); axA.set_ylim(0.0, 1.06 * ymax_a)
    overlay_weight_linear(axA, lam_sw, sw_weight(lam_sw, T_sw, cut), PANEL_A_XLIM)
    axA.axvline(cut, color='0.45', ls='--', lw=0.7, alpha=0.6, zorder=0.6)
    axA.set_xlabel(r'$\lambda$ [$\mu$m]')
    axA.set_ylabel(SW_YLABEL)
    axA.legend(loc='upper right')
    axA.grid(alpha=0.3, zorder=0)

    ymax_b = 0.0
    for m in MATS:
        y = fac * D[f'{m}_qabs_lw']
        sel = np.isfinite(y) & (y > 0)
        axB.plot(lam_lw[sel], y[sel], color=COLORS[m], zorder=3, label=m)
        in_x = sel & (lam_lw >= PANEL_B_XLIM[0]) & (lam_lw <= PANEL_B_XLIM[1])
        ymax_b = max(ymax_b, float(np.nanmax(y[in_x])))
    axB.set_xlim(*PANEL_B_XLIM)
    wlw = lw_weight_curve(lam_lw, windows)
    if panel_b_logy:
        axB.set_yscale('log'); axB.set_ylim(PANEL_B_YMIN, 1.6 * ymax_b)
        overlay_weight_logpanel(axB, lam_lw, wlw, PANEL_B_XLIM)
        axB.grid(alpha=0.3, which='major', zorder=0); axB.grid(alpha=0.12, which='minor', zorder=0)
    else:
        axB.set_ylim(0.0, 1.06 * ymax_b)
        overlay_weight_linear(axB, lam_lw, wlw, PANEL_B_XLIM)
        axB.grid(alpha=0.3, zorder=0)
    for lo, hi, T_k, wk in windows:
        if PANEL_B_XLIM[0] < lo < PANEL_B_XLIM[1]:
            axB.axvline(lo, color='0.55', lw=0.6, alpha=0.5, zorder=0.6)
    axB.set_xlabel(r'$\lambda$ [$\mu$m]')
    axB.set_ylabel(LW_YLABEL)
    axB.legend(loc='upper right')

    for ax, lab in ((axA, 'a'), (axB, 'b')):
        ax.text(0.025, 0.96, f'({lab})', transform=ax.transAxes, ha='left', va='top',
                fontsize=12, fontweight='bold', zorder=5)

    for ext in ('pdf', 'png'):
        fig.savefig(f'{out_base}.{ext}', dpi=300, bbox_inches='tight')
        print(f'saved {out_base}.{ext}')
    plt.close(fig)


def main():
    D = dict(np.load(HERE / 'mie_d500_data.npz', allow_pickle=True))
    make_fig(D, str(HERE / 'fig'),       panel_b_logy=True)
    make_fig(D, str(HERE / 'fig_linB'),  panel_b_logy=False)


if __name__ == '__main__':
    main()
