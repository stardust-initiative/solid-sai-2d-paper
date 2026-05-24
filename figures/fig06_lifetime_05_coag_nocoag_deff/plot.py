#!/usr/bin/env python3
"""Fig 5 (fig:lifetime_coag_nocoag_deff) — §3.2.1 Role of coagulation.

1 row x 3 cols:
  (a) Lifetime, no coagulation, silica d_0 = 0.5 um (years)
  (b) Lifetime, with coagulation D_f=1.6, silica d_0 = 0.5 um (years)
  (c) Effective diameter d_eff (um) with coagulation

All three panels: full diverging vik (Crameri) with TwoSlopeNorm centred at
1.0 (1 yr threshold for lifetimes; 1 um threshold for d_eff). Values below
the centre render in blue; values above render in red. White at the centre
gives a clear "1.0 isoline" reading on each panel.

ERA5 10-year climatological tropopause (2008-2017) overlay: 2.2 pt white under 1.0 pt black.

Self-contained: reads data.nc from the same folder, writes fig.png to
the same folder. cmap (full vik) defined inline.

Usage:
    cd manuscript_figures/fig05_lifetime_05_coag_nocoag_deff
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
from _style import apply, WARM_SEQ, WARM_NOWHITE, panel_label, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))

ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
inj_lat  = ds.inj_lat.values
inj_pres = ds.inj_pres.values
tp_lat   = ds.tp_lat.values
tp_pres  = ds.tp_pres_clim.values
tau_nc   = ds.lifetime_nocoag.values
tau_co   = ds.lifetime_coag.values
d_eff_co = ds.d_eff_coag.values
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


# Lifetime panels: sequential warm with white-at-zero; envelope [0, 2.5] yr.
# 10-step linspace gives 9 discrete colour bands between adjacent levels.
levels_tau = np.linspace(0.0, 2.5, 10)
labels_tau = [0.5, 1.0, 1.5, 2.0]
norm_tau   = Normalize(vmin=0.0, vmax=2.5)
ticks_tau  = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]

# d_eff panel: sequential warm without white anchor; envelope [0.4, 1.4] um.
levels_de = np.linspace(0.4, 1.4, 10)
labels_de = [0.6, 0.8, 1.0, 1.2]
norm_de   = Normalize(vmin=0.4, vmax=1.4)
ticks_de  = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]

fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), sharey=True)

# Panel (a) — lifetime no coag
ax = axes[0]
data = np.clip(tau_nc, levels_tau[0], levels_tau[-1])
pm_a = ax.contourf(inj_lat, inj_pres, data, levels=levels_tau,
                   cmap=WARM_SEQ, norm=norm_tau, extend='max')
cs = ax.contour(inj_lat, inj_pres, data, levels=labels_tau,
                colors='k', linewidths=0.6, alpha=0.85)
ax.clabel(cs, inline=True, fmt='%.2f', fontsize=8)
fmt_axes(ax, ylabel=True, altaxis=False)
ax.set_title('Lifetime, no coagulation', fontsize=11)
panel_label(ax, 'a')

# Panel (b) — lifetime with coag
ax = axes[1]
data = np.clip(tau_co, levels_tau[0], levels_tau[-1])
pm_b = ax.contourf(inj_lat, inj_pres, data, levels=levels_tau,
                   cmap=WARM_SEQ, norm=norm_tau, extend='max')
cs = ax.contour(inj_lat, inj_pres, data, levels=labels_tau,
                colors='k', linewidths=0.6, alpha=0.85)
ax.clabel(cs, inline=True, fmt='%.2f', fontsize=8)
fmt_axes(ax, ylabel=False, altaxis=False)
ax.set_title(r'Lifetime, with coagulation $(D_f{=}1.6,\ 20\ \mathrm{Tg\ yr}^{-1})$',
             fontsize=11)
panel_label(ax, 'b')

# Panel (c) — d_eff with coag
ax = axes[2]
data = np.clip(d_eff_co, levels_de[0], levels_de[-1])
pm_c = ax.contourf(inj_lat, inj_pres, data, levels=levels_de,
                   cmap=WARM_NOWHITE, norm=norm_de, extend='max')
cs = ax.contour(inj_lat, inj_pres, data, levels=labels_de,
                colors='k', linewidths=0.6, alpha=0.85)
ax.clabel(cs, inline=True, fmt='%.2f', fontsize=8)
fmt_axes(ax, ylabel=False, altaxis=True)
ax.set_title(r'Effective diameter $d_{\rm eff}\ (\mu\mathrm{m})$, with coagulation',
             fontsize=11)
panel_label(ax, 'c')

# Two colorbars: shared lifetime cbar between panels b and c, d_eff cbar far right
fig.subplots_adjust(left=0.05, right=0.88, top=0.93, bottom=0.13, wspace=0.30)
cax_tau = fig.add_axes([0.608, 0.13, 0.010, 0.79])
cax_de  = fig.add_axes([0.92,  0.13, 0.010, 0.79])
cb_tau  = fig.colorbar(pm_b, cax=cax_tau, ticks=ticks_tau)
cb_de   = fig.colorbar(pm_c, cax=cax_de, ticks=ticks_de)
cb_tau.ax.yaxis.set_label_position('left')
cb_tau.set_label('Lifetime (yr)')
cb_de.set_label(r'$d_{\rm eff}\ (\mu\mathrm{m})$')

savefig_both(fig, HERE)
