#!/usr/bin/env python3
"""Fig A2 (fig:seasonal_Psi_resTEM) — Appendix A.

Seasonal-mean TEM residual stream function from ERA5 (2008-2017
climatology), aggregated to four seasons (DJF/MAM/JJA/SON).

Diverging field around 0 (Brewer-Dobson NH cell positive, SH cell
negative). Full diverging vik (Crameri) with TwoSlopeNorm centred at 0,
symmetric envelope ±3000 m Pa s^-1. Rendered with contourf over 21
evenly-spaced levels for smooth filled contours; black isoline overlay
at +/-1500, +/-500.

Self-contained: reads data.nc, writes fig.png, both in the same folder.
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import xarray as xr
from matplotlib.ticker import FixedLocator

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, DIV_BWR, centered_norm, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))

ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
lat       = ds.latitude.values
lev       = ds.level.values
seasonal  = ds.Psi_res_TEM.values
sname     = ds.season_name.values
ds.close()

vmax = 3000.0
levels = np.linspace(-vmax, vmax, 21)
contour_levels = np.array([-1500, -500, 500, 1500])
norm = centered_norm(-vmax, vmax, vcenter=0.0)

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
# Panel order: (a) DJF top-left, (b) MAM top-right,
#              (c) JJA bottom-left, (d) SON bottom-right.
positions = [('a', 'DJF', axes[0, 0]),
             ('b', 'MAM', axes[0, 1]),
             ('c', 'JJA', axes[1, 0]),
             ('d', 'SON', axes[1, 1])]
pms = []
for letter, s, ax in positions:
    si = int(np.where(sname == s)[0][0])
    Z = np.clip(seasonal[si], -vmax, vmax)
    pm = ax.contourf(lat, lev, Z, levels=levels, cmap=DIV_BWR,
                     norm=norm, extend='both')
    cs = ax.contour(lat, lev, Z, levels=contour_levels,
                    colors='k', linewidths=0.6, alpha=0.85)
    ax.clabel(cs, inline=True, fmt=lambda v: f'{int(v):d}', fontsize=8)
    ax.set_yscale('log'); ax.invert_yaxis()
    ax.set_ylim(100, 10)
    ax.set_xlim(-90, 90)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.yaxis.set_major_locator(FixedLocator([10, 30, 100]))
    ax.yaxis.set_minor_locator(FixedLocator([20, 40, 50, 60, 70, 80]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{int(v)}'))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_title(f'({letter}) {s}', fontsize=11)
    pms.append(pm)

for ax in axes[1, :]:
    ax.set_xlabel('Latitude (°)')
for ax in axes[:, 0]:
    ax.set_ylabel('Pressure (hPa)')

fig.suptitle(r'Seasonal Mean $\Psi_{\rm res,\,TEM}$ — 10-year mean (2008–2017)',
             fontsize=12)
fig.subplots_adjust(left=0.08, right=0.88, top=0.92, bottom=0.08,
                    hspace=0.22, wspace=0.10)
cax = fig.add_axes([0.91, 0.10, 0.013, 0.78])
cb  = fig.colorbar(pms[0], cax=cax,
                   ticks=[-3000, -2000, -1000, 0, 1000, 2000, 3000])
cb.set_label(r'$\Psi_{\rm res,\,TEM}$ (m Pa s$^{-1}$)')

savefig_both(fig, HERE)
