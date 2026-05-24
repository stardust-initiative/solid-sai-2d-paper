#!/usr/bin/env python3
"""Fig 8 (fig:ratio_03_07_lifetime_coag) — §3.2.3 Role of monomer size.

Sensitivity of silica with coagulation to monomer core diameter.
  (a) Lifetime ratio  tau(d_0=0.3) / tau(d_0=0.7)  — all >= 1
  (b) d_eff ratio     d_eff(d_0=0.3) / d_eff(d_0=0.7)  — crosses 1.0

Panel (a) uses WARM_NOWHITE (sequential warm, no white anchor); 8 discrete
bands of width 0.1 between 1.0 and 1.8.

Panel (b) uses DIV_BWR diverging, 12 discrete bands of width 0.1 from 0.4 to
1.6 with the boundary EXACTLY at 1.0 — the pale-blue band [0.9, 1.0) meets
the pale-yellow band [1.0, 1.1) sharply at the no-effect threshold (no broad
central white band).

Self-contained: reads data.nc from the same folder, writes fig.png to
the same folder.
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import xarray as xr
from matplotlib.ticker import FixedLocator
from matplotlib.colors import Normalize

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, WARM_NOWHITE, DIV_BWR, panel_label, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))

ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
inj_lat   = ds.inj_lat.values
inj_pres  = ds.inj_pres.values
ratio_tau = ds.ratio_lifetime.values
ratio_de  = ds.ratio_d_eff.values
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


# Panel (a) — lifetime ratio, all >= 1.
VMIN_A, VMAX_A = 1.0, 1.8
norm_a   = Normalize(vmin=VMIN_A, vmax=VMAX_A)
levels_a = np.linspace(VMIN_A, VMAX_A, 9)      # 8 bands of width 0.1
ticks_a  = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]
clev_a   = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])

# Panel (b) — d_eff ratio crossing 1.0. 12 bands width 0.1, boundary at 1.0.
BOUNDARIES_B = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                          1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
COLORS_B = [DIV_BWR(p) for p in np.linspace(0.0, 1.0, len(BOUNDARIES_B) - 1)]
ticks_b  = list(BOUNDARIES_B)
clev_b   = np.array([0.5, 0.7, 0.9, 1.0, 1.1])

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)

# Panel (a)
ax = axes[0]
data = np.clip(ratio_tau, VMIN_A, VMAX_A)
pm_a = ax.contourf(inj_lat, inj_pres, data, levels=levels_a,
                   cmap=WARM_NOWHITE, norm=norm_a, extend='max')
cs   = ax.contour(inj_lat, inj_pres, data, levels=clev_a,
                  colors='k', linewidths=0.6, alpha=0.85)
ax.clabel(cs, inline=True, fmt='%.1f', fontsize=8)
fmt_axes(ax, ylabel=True, altaxis=False)
ax.set_title(r'Lifetime ratio  $\tau(d_0{=}0.3\,\mu\mathrm{m})/\tau(d_0{=}0.7\,\mu\mathrm{m})$',
             fontsize=11)
panel_label(ax, 'a')

# Panel (b)
ax = axes[1]
data = np.clip(ratio_de, BOUNDARIES_B[0], BOUNDARIES_B[-1])
pm_b = ax.contourf(inj_lat, inj_pres, data, levels=BOUNDARIES_B,
                   colors=COLORS_B, extend='both')
cs   = ax.contour(inj_lat, inj_pres, data, levels=clev_b,
                  colors='k', linewidths=0.6, alpha=0.85)
ax.clabel(cs, inline=True, fmt='%.1f', fontsize=8)
fmt_axes(ax, ylabel=False, altaxis=True)
ax.set_title(r'Effective-diameter ratio  $d_{\rm eff}(d_0{=}0.3)/d_{\rm eff}(d_0{=}0.7)$',
             fontsize=11)
panel_label(ax, 'b')

fig.subplots_adjust(left=0.07, right=0.82, top=0.92, bottom=0.13, wspace=0.40)
cax_a = fig.add_axes([0.395, 0.13, 0.014, 0.79])
cax_b = fig.add_axes([0.92,  0.13, 0.014, 0.79])
cb_a  = fig.colorbar(pm_a, cax=cax_a, ticks=ticks_a)
cb_b  = fig.colorbar(pm_b, cax=cax_b, ticks=ticks_b)
cb_a.set_label(r'$\tau(0.3)/\tau(0.7)$')
cb_b.set_label(r'$d_{\rm eff}(0.3)/d_{\rm eff}(0.7)$')

savefig_both(fig, HERE)
