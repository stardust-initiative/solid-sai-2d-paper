#!/usr/bin/env python3
"""Fig 8 (fig:swlw_split) — §3.3 Radiative forcing and stratospheric heating.

2 row x 3 col layout:
  Row 1 (top):    SW SARF = global cosine-weighted ASR anomaly  (negative = cooling)
  Row 2 (bottom): LW SARF = global cosine-weighted -OLR anomaly (positive = warming)
  Cols:           sulfate d_0=0.5 um, silica d_0=0.5 um, calcite d_0=0.3 um

Both rows: full diverging vik (Crameri) with TwoSlopeNorm. Same envelope and
centre per row, so the three columns of each row are directly eye-comparable.
  SW row: envelope -15..0, centre -7.5  (data midpoint = -7.15; rounded)
  LW row: envelope   0..2.4, centre  1.2 (data midpoint =  1.16; rounded)

ERA5 10-year climatological tropopause (2008-2017) overlay: 2.2 pt white under 1.0 pt black.
Standalone (a)-(f) panel tags upper-left, white bbox.
Panel pairing: (a,b) sulfate, (c,d) silica, (e,f) calcite.

Self-contained: reads data.nc, writes fig.png, both in the same folder.

Usage:
    cd manuscript_figures/fig08_swlw_split
    python plot.py
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import xarray as xr
from matplotlib.ticker import FixedLocator
from matplotlib.colors import Normalize

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, COOL_SEQ, WARM_SEQ, panel_label, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))

ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
lat_inj  = ds.lat_inj.values
lev_inj  = ds.lev_inj.values
SW       = ds.SW.values    # (material, lev_inj, lat_inj)
LW       = ds.LW.values
mat_name = ds.material_name.values
mat_size = ds.material_size.values
tp_lat   = ds.tp_lat.values
tp_pres  = ds.tp_pres_clim.values
ds.close()

P_MAX = 120
ALT_TICKS = [15, 17, 19, 21]
SW_TAGS = ['(a)', '(c)']
LW_TAGS = ['(b)', '(d)']
TITLES = [
    r'Silica, $d_0 = 0.5\ \mu$m',
    r'Calcite, $d_0 = 0.3\ \mu$m',
]


def p_to_z(p): return 7.0 * np.log(1000.0 / p)
def z_to_p(z): return 1000.0 * np.exp(-z / 7.0)


def fmt_axes(ax, ylabel=False, altaxis=False, xlabel=True):
    ax.set_yscale('log'); ax.invert_yaxis()
    ax.set_ylim(P_MAX * 1.02, lev_inj.min() * 0.98)
    ax.set_xlim(lat_inj.min() - 3, lat_inj.max() + 3)
    ax.set_xticks(np.arange(-60, 61, 20))
    if xlabel:
        ax.set_xlabel('Injection latitude (°)')
    ax.yaxis.set_major_locator(FixedLocator([50, 70, 100, 120]))
    ax.yaxis.set_minor_locator(FixedLocator([60, 80, 90, 110]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{int(v)}'))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.plot(tp_lat, tp_pres, 'w-', lw=2.2)
    ax.plot(tp_lat, tp_pres, 'k-', lw=1.0)
    if ylabel:
        ax.set_ylabel('Injection pressure (hPa)')
    if altaxis:
        ax2 = ax.secondary_yaxis('right', functions=(p_to_z, z_to_p))
        ax2.set_ylabel('Altitude (km)', labelpad=12)
        ax2.set_yticks(ALT_TICKS)
        ax2.set_yticklabels([str(v) for v in ALT_TICKS])
        ax2.minorticks_off()


# Per-row sequential envelopes (SW <= 0, LW >= 0).
# SW vmax kept at -1 (not 0) so the densely-populated near-zero band reads
# as the lightest in-scale colour rather than collapsing to white. extend='max'
# caps the few cells in (-1, 0] to the off-scale white.
levels_sw = np.arange(-12.0, -1.0 + 0.01, 1.0)
contour_sw = np.array([-11, -9, -7, -5, -3])
norm_sw   = Normalize(vmin=-12.0, vmax=-1.0)

levels_lw = np.arange(0.0, 2.41, 0.2)
contour_lw = np.array([0.4, 0.8, 1.2, 1.6, 2.0])
norm_lw   = Normalize(vmin=0.0, vmax=2.4)

fig, axes = plt.subplots(2, 2, figsize=(12.5, 9), sharey=True)

for c_idx in range(SW.shape[0]):
    sw_grid = np.clip(SW[c_idx], levels_sw[0], levels_sw[-1])
    lw_grid = np.clip(LW[c_idx], levels_lw[0], levels_lw[-1])

    # SW (top row)
    ax = axes[0, c_idx]
    pm_top = ax.contourf(lat_inj, lev_inj, sw_grid,
                         levels=levels_sw, cmap=COOL_SEQ, norm=norm_sw, extend='max')
    cs = ax.contour(lat_inj, lev_inj, sw_grid, levels=contour_sw,
                    colors='k', linewidths=0.6, alpha=0.85)
    ax.clabel(cs, inline=True, fmt='%.0f', fontsize=8)
    fmt_axes(ax,
             ylabel=(c_idx == 0),
             altaxis=(c_idx == SW.shape[0] - 1),
             xlabel=False)
    ax.set_title(TITLES[c_idx], fontsize=11)
    panel_label(ax, SW_TAGS[c_idx].strip('()'))

    # LW (bottom row)
    ax = axes[1, c_idx]
    pm_bot = ax.contourf(lat_inj, lev_inj, lw_grid,
                         levels=levels_lw, cmap=WARM_SEQ, norm=norm_lw)
    cs = ax.contour(lat_inj, lev_inj, lw_grid, levels=contour_lw,
                    colors='k', linewidths=0.6, alpha=0.85)
    ax.clabel(cs, inline=True, fmt='%.1f', fontsize=8)
    fmt_axes(ax,
             ylabel=(c_idx == 0),
             altaxis=(c_idx == SW.shape[0] - 1),
             xlabel=True)
    panel_label(ax, LW_TAGS[c_idx].strip('()'))

# Two row colorbars (one per row) sitting to the right of the altitude axes.
fig.subplots_adjust(left=0.06, right=0.84, top=0.95, bottom=0.08,
                    hspace=0.30, wspace=0.22)
cax_top = fig.add_axes([0.92, 0.55, 0.013, 0.40])
cb_top  = fig.colorbar(pm_top, cax=cax_top,
                       ticks=[-12, -10, -8, -6, -4, -2, -1])
cb_top.set_label(r'SW SARF: ASR anomaly (W m$^{-2}$)')

cax_bot = fig.add_axes([0.92, 0.08, 0.013, 0.40])
cb_bot  = fig.colorbar(pm_bot, cax=cax_bot,
                       ticks=[0, 0.6, 1.2, 1.8, 2.4])
cb_bot.set_label(r'LW SARF: $-$OLR anomaly (W m$^{-2}$)')

savefig_both(fig, HERE)
