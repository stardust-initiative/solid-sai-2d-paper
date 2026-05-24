#!/usr/bin/env python3
"""Plot K_phi_phi seasonal comparison: ERA5 (2008-2017) vs GSFC 2-D (year 2000).

Diverging DIV_BWR (ACP f02(a) prototype) on a log-spaced half-decade grid,
centred on the boundary at 1.0 (= 10^5 m^2/s in displayed units). Although
K_phi_phi is strictly positive, the geometric midpoint at 10^5 m^2/s
separates the tropical-pipe regime (blue, suppressed mixing) from the
surf-zone mid-latitude mixing regime (red).

Six discrete bands at log-spaced boundaries [0.1, 0.316, 1.0, 3.16, 10, 31.6,
100]. Tick labels rounded to 0.1, 0.3, 1, 3, 10, 30, 100 for readability.

Self-contained: reads data.nc, writes fig.png.
"""
import os, sys, numpy as np, xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import FixedLocator

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, DIV_BWR, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))

def p_to_z(p): return 7.0 * np.log(1000.0 / p)
def z_to_p(z): return 1000.0 * np.exp(-z / 7.0)
ALT_TICKS = [14, 16, 18, 20, 22, 24, 26, 28, 30, 32]

ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
seasons   = [str(s) for s in ds.season.values]
era_K     = ds.D_phi_phi_ERA5.values
era_lat   = ds.era_lat.values
era_lev   = ds.era_level.values
gsfc_K    = ds.D_phi_phi_GSFC.values
gsfc_lat  = ds.gsfc_lat.values
gsfc_lev  = ds.gsfc_level.values
tp_lat    = ds.tp_lat.values
tp_pres   = ds.tp_pres_clim.values
ds.close()

BOUNDARIES = np.array([0.1, 0.316, 1.0, 3.162, 10.0, 31.62, 100.0])
TICK_LABELS = ['0.1', '0.3', '1', '3', '10', '30', '100']
LOG_LO, LOG_MID, LOG_HI = -1.0, 0.0, 2.0

def log_pos(v):
    lv = np.log10(v)
    if lv < LOG_MID:
        return 0.5 * (lv - LOG_LO) / (LOG_MID - LOG_LO)
    return 0.5 + 0.5 * (lv - LOG_MID) / (LOG_HI - LOG_MID)

GEO_CENTERS = np.sqrt(BOUNDARIES[:-1] * BOUNDARIES[1:])
COLORS = [DIV_BWR(log_pos(c)) for c in GEO_CENTERS]
contour_levels = np.array([0.3, 1, 3, 10, 30])

fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), sharex=True, sharey=True,
                         gridspec_kw=dict(left=0.07, right=0.87, top=0.90, bottom=0.09,
                                          wspace=0.10, hspace=0.20))

def panel(ax, lat, lev, K, title=None, alt_axis=False):
    Z = np.clip(K, BOUNDARIES[0], BOUNDARIES[-1])
    pm = ax.contourf(lat, lev, Z, levels=BOUNDARIES, colors=COLORS, extend='both')
    cs = ax.contour(lat, lev, Z, levels=contour_levels,
                    colors='k', linewidths=0.6, alpha=0.85)
    ax.clabel(cs, inline=True, fmt=lambda v: f'{v:g}', fontsize=8)
    ax.plot(tp_lat, tp_pres, 'k-', lw=2.2)
    ax.set_yscale('log'); ax.invert_yaxis()
    ax.set_ylim(150, 10)
    ax.set_xlim(-90, 90)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.yaxis.set_major_locator(FixedLocator([10, 30, 100, 150]))
    ax.yaxis.set_minor_locator(FixedLocator([20, 40, 50, 60, 70, 80, 90, 110, 120, 130, 140]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{int(v)}'))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    if title: ax.set_title(title, fontsize=11)
    if alt_axis:
        ax2 = ax.secondary_yaxis('right', functions=(p_to_z, z_to_p))
        ax2.set_ylabel('Altitude (km)', labelpad=10)
        ax2.set_yticks(ALT_TICKS)
        ax2.set_yticklabels([str(v) for v in ALT_TICKS])
        ax2.minorticks_off()
    return pm

pm = None
for col, s in enumerate(seasons):
    is_last = (col == len(seasons) - 1)
    pm = panel(axes[0, col], era_lat,  era_lev,  era_K[col],  title=s, alt_axis=is_last)
    panel(axes[1, col], gsfc_lat, gsfc_lev, gsfc_K[col], alt_axis=is_last)

for ax in axes[1, :]: ax.set_xlabel('Latitude (deg)')
for ax in axes[:, 0]: ax.set_ylabel('Pressure (hPa)')

fig.text(0.012, 0.70, 'ERA5\n(2008-2017 clim)', ha='center', va='center', rotation=90, fontsize=11)
fig.text(0.012, 0.29, 'GSFC 2-D\n(year 2000)',  ha='center', va='center', rotation=90, fontsize=11)

cax = fig.add_axes([0.93, 0.10, 0.013, 0.78])
cb  = fig.colorbar(pm, cax=cax, ticks=BOUNDARIES)
cb.ax.set_yticklabels(TICK_LABELS)
cb.set_label(r'$D_{\varphi\varphi}\ (10^{5}\ \mathrm{m}^{2}\ \mathrm{s}^{-1})$')

savefig_both(fig, HERE)
