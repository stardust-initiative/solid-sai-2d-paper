#!/usr/bin/env python3
"""Fig 10 (fig:silica_d0_optimum) — §3.3.

SARF ratio silica d_0=0.3 / d_0=0.5 side by side:
  (a) no coagulation
  (b) with coagulation

Diverging DIV_BWR with explicit 12-colour sampling. Boundaries at 0.4..1.6
in steps of 0.1 with the boundary EXACTLY on 1.0 — the pale-blue band
[0.9, 1.0) meets the pale-yellow band [1.0, 1.1) sharply at the no-effect
threshold (no broad central white band). Shared envelope for direct
eye-comparability.

Self-contained: reads data.nc, writes fig.png, both in the same folder.
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import xarray as xr
from matplotlib.ticker import FixedLocator

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, DIV_BWR, panel_label, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))

ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
inj_lat   = ds.inj_lat.values
inj_pres  = ds.inj_pres.values
ratio_nc  = ds.ratio_nocoag.values
ratio_co  = ds.ratio_coag.values
tp_lat    = ds.tp_lat.values
tp_pres   = ds.tp_pres_clim.values
ds.close()

P_MAX = 120
ALT_TICKS = [15, 17, 19, 21]


def p_to_z(p): return 7.0 * np.log(1000.0 / p)
def z_to_p(z): return 1000.0 * np.exp(-z / 7.0)


def fmt_axes(ax, ylabel=False, altaxis=False):
    ax.set_yscale('log'); ax.invert_yaxis()
    ax.set_ylim(P_MAX * 1.02, inj_pres.min() * 0.98)
    ax.set_xlim(inj_lat.min() - 3, inj_lat.max() + 3)
    ax.set_xticks(np.arange(-60, 61, 20))
    ax.set_xlabel('Injection latitude (deg)')
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


BOUNDARIES = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
COLORS = [DIV_BWR(p) for p in np.linspace(0.0, 1.0, len(BOUNDARIES) - 1)]
ticks = list(BOUNDARIES)
contour_levs = np.array([0.6, 0.8, 1.0, 1.2, 1.4])

fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
panels = [
    (axes[0], ratio_nc, 'a', 'No coagulation'),
    (axes[1], ratio_co, 'b', 'With coagulation'),
]
pm = None
for c_idx, (ax, data, tag, title) in enumerate(panels):
    data_c = np.clip(data, BOUNDARIES[0], BOUNDARIES[-1])
    pm = ax.contourf(inj_lat, inj_pres, data_c, levels=BOUNDARIES,
                     colors=COLORS, extend='both')
    cs = ax.contour(inj_lat, inj_pres, data_c, levels=contour_levs,
                    colors='k', linewidths=0.6, alpha=0.85)
    ax.clabel(cs, inline=True, fmt='%.1f', fontsize=8)
    fmt_axes(ax, ylabel=(c_idx == 0), altaxis=(c_idx == len(panels) - 1))
    ax.set_title(title, fontsize=11)
    panel_label(ax, tag)

fig.subplots_adjust(right=0.88, wspace=0.12)
cax = fig.add_axes([0.933, 0.18, 0.018, 0.66])
cb  = fig.colorbar(pm, cax=cax, ticks=ticks)
cb.set_label(r'Net SARF ratio: $d_0=0.3 / d_0=0.5$')

savefig_both(fig, HERE)
