#!/usr/bin/env python3
"""Fig 7 (fig:ratio_Df16_Df30_lifetime) — §3.2.2 Role of agglomerate shape.

Single-panel: lifetime ratio tau(D_f=1.6)/tau(D_f=3.0), silica d_0=0.5 um
with coagulation. Data is essentially all >= 1 (loose aggregates outlive
compact ones across nearly the whole grid), so we use a sequential warm
colormap rather than a diverging one. WARM_NOWHITE (pale yellow at the low
end, no white anchor) so the lowest band reads pale-yellow rather than
white — this avoids the confusion of the central white band colliding with
the white NaN cut-out cells where no injection runs were performed.

8 discrete bands of width 0.05 between 1.00 and 1.40.

Layout adjusted so the colorbar + tick labels clear the secondary altitude
axis "Altitude (km)" label.

Self-contained: reads data.nc, writes fig.png, both in the same folder.

Usage:
    cd manuscript_figures/fig06_ratio_Df16_Df30_lifetime
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
from _style import apply, WARM_NOWHITE, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))

ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
inj_lat  = ds.inj_lat.values
inj_pres = ds.inj_pres.values
ratio    = ds.ratio.values
tp_lat   = ds.tp_lat.values
tp_pres  = ds.tp_pres_clim.values
ds.close()

P_MAX = 120
ALT_TICKS = [15, 17, 19, 21]


def p_to_z(p): return 7.0 * np.log(1000.0 / p)
def z_to_p(z): return 1000.0 * np.exp(-z / 7.0)


# Sequential warm scale; vmin=1.0 (no-effect reference), vmax=1.40.
# 9 levels = 8 discrete bands of width 0.05. Lowest band [1.00, 1.05) reads
# pale yellow under WARM_NOWHITE — explicitly NOT white, to keep NaN cut-out
# cells (which render as white) visually distinct from "ratio just above 1".
VMIN, VMAX = 1.0, 1.4
norm   = Normalize(vmin=VMIN, vmax=VMAX)
levels = np.linspace(VMIN, VMAX, 9)
ticks  = [1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4]
clev   = np.array([1.05, 1.10, 1.20, 1.30])
clev   = clev[(clev >= VMIN) & (clev <= VMAX)]

fig, ax = plt.subplots(1, 1, figsize=(8, 5.2))
data_c = np.clip(ratio, VMIN, VMAX)
pm = ax.contourf(inj_lat, inj_pres, data_c, levels=levels,
                 cmap=WARM_NOWHITE, norm=norm, extend='max')
cs = ax.contour(inj_lat, inj_pres, data_c, levels=clev,
                colors='k', linewidths=0.6, alpha=0.85)
ax.clabel(cs, inline=True, fmt='%.2f', fontsize=8)

# Standard formatting
ax.set_yscale('log'); ax.invert_yaxis()
ax.set_ylim(P_MAX * 1.02, inj_pres.min() * 0.98)
ax.set_xlim(inj_lat.min() - 3, inj_lat.max() + 3)
ax.set_xticks(np.arange(-60, 61, 20))
ax.set_xlabel('Injection latitude (deg)')
ax.set_ylabel('Injection pressure (hPa)')
ax.yaxis.set_major_locator(FixedLocator([50, 70, 100, 120]))
ax.yaxis.set_minor_locator(FixedLocator([60, 80, 90, 110]))
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{int(v)}'))
ax.yaxis.set_minor_formatter(mticker.NullFormatter())
ax.plot(tp_lat, tp_pres, 'w-', lw=2.2)
ax.plot(tp_lat, tp_pres, 'k-', lw=1.0)

# Altitude (km) secondary axis
ax2 = ax.secondary_yaxis('right', functions=(p_to_z, z_to_p))
ax2.set_ylabel('Altitude (km)', labelpad=12)
ax2.set_yticks(ALT_TICKS)
ax2.set_yticklabels([str(v) for v in ALT_TICKS])
ax2.minorticks_off()

ax.set_title(r'$\tau(D_f{=}1.6) / \tau(D_f{=}3.0)$, silica 0.5 $\mu$m, with coag',
             fontsize=11)

# Layout: panel right edge at 0.78 leaves room for altitude (km) label;
# colorbar shifted right so it (and its tick labels) clear the altitude label.
fig.subplots_adjust(left=0.10, right=0.78, top=0.92, bottom=0.13)
cax = fig.add_axes([0.92, 0.13, 0.018, 0.79])
cb  = fig.colorbar(pm, cax=cax, ticks=ticks)
cb.set_label(r'$\tau(D_f=1.6)/\tau(D_f=3.0)$')

savefig_both(fig, HERE)
