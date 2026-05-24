#!/usr/bin/env python3
"""Fig 4 (fig:base_fn_examples) — three injection cases, MMR + column density."""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import xarray as xr
from matplotlib.ticker import FixedLocator, LogFormatterMathtext
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, WARM_SEQ, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))
ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
lag_lat=ds.lag_lat.values; lag_pres=ds.lag_pres.values
twod_lat=ds.twod_lat.values; twod_pres=ds.twod_pres.values
case_label=ds.case_label.values
inj_lat_v=ds.inj_lat.values; inj_pres_v=ds.inj_pres.values
B_LAG_v=ds.B_LAG.values; B_2D_v=ds.B_2D.values; ratio_v=ds.ratio.values
lag_mmr=ds.lag_mmr.values; twod_mmr=ds.twod_mmr.values
lag_cd=ds.lag_col_cd.values; twod_cd=ds.twod_col_cd.values
tp_lat=ds.tp_lat.values; tp_pres=ds.tp_pres_2010.values
ds.close()

LOGNORM_DECADES = 4.0
N_FILL_LEVELS = 9

def add_tropopause(ax):
    ax.plot(tp_lat, tp_pres, 'w-', lw=2.2, zorder=4)
    ax.plot(tp_lat, tp_pres, 'k-', lw=1.0, zorder=4)

fig, axes = plt.subplots(3, 3, figsize=(15, 10),
                         gridspec_kw={'width_ratios': [1, 1, 0.9]})

for r, lbl in enumerate(case_label):
    lag_p = lag_mmr[r]; twod_p = twod_mmr[r]
    peak = max(np.nanmax(lag_p), np.nanmax(twod_p))
    vmin = peak * 10.0**(-LOGNORM_DECADES); vmax = peak
    norm = LogNorm(vmin=vmin, vmax=vmax)
    levels = np.logspace(np.log10(vmin), np.log10(vmax), N_FILL_LEVELS)
    contour_levels = peak * np.array([1e-2, 3e-2, 1e-1, 3e-1])
    case_label_full = f'({"i"*(r+1)}) {lbl}'

    # (a) Lagranto
    ax = axes[r, 0]
    Z = np.clip(lag_p, vmin, vmax)
    pm_a = ax.contourf(lag_lat, lag_pres, Z, levels=levels,
                       cmap=WARM_SEQ, norm=norm, extend='both')
    ax.contour(lag_lat, lag_pres, Z, levels=contour_levels,
               colors='k', linewidths=0.5, alpha=0.7, norm=norm)
    add_tropopause(ax)
    ax.set_yscale('log'); ax.invert_yaxis()
    ax.set_xlim(-90, 90); ax.set_xticks([-90,-60,-30,0,30,60,90])
    ax.set_ylim(200, 10); ax.yaxis.set_major_locator(FixedLocator([10,30,100,200]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v,_: f'{int(v)}'))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.scatter([inj_lat_v[r]], [inj_pres_v[r]], marker='*', s=180, c='white',
               edgecolor='k', linewidth=0.8, zorder=5)
    ax.set_ylabel(f'{case_label_full}\nPressure (hPa)', fontsize=10)
    if r == 0:
        ax.set_title('(a) Lagranto', fontsize=11)
    if r == 2:
        ax.set_xlabel('Latitude (°)')
    cax = ax.inset_axes([1.02, 0.0, 0.04, 1.0])
    cb = fig.colorbar(pm_a, cax=cax)
    dlo = int(np.ceil(np.log10(vmin))); dhi = int(np.floor(np.log10(vmax)))
    cb.set_ticks(10.0**np.arange(dlo, dhi+1))
    cb.ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    cb.set_label('mmr (mg/kg)', fontsize=8, labelpad=2)
    cb.ax.tick_params(labelsize=7, pad=1.5)

    # (b) 2-D model
    ax = axes[r, 1]
    Z = np.clip(twod_p, vmin, vmax)
    ax.contourf(twod_lat, twod_pres, Z, levels=levels,
                cmap=WARM_SEQ, norm=norm, extend='both')
    ax.contour(twod_lat, twod_pres, Z, levels=contour_levels,
               colors='k', linewidths=0.5, alpha=0.7, norm=norm)
    add_tropopause(ax)
    ax.set_yscale('log'); ax.invert_yaxis()
    ax.set_xlim(-90, 90); ax.set_xticks([-90,-60,-30,0,30,60,90])
    ax.set_ylim(200, 10); ax.yaxis.set_major_locator(FixedLocator([10,30,100,200]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v,_: f'{int(v)}'))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.scatter([inj_lat_v[r]], [inj_pres_v[r]], marker='*', s=180, c='white',
               edgecolor='k', linewidth=0.8, zorder=5)
    if r == 0:
        ax.set_title('(b) 2-D model', fontsize=11)
    if r == 2:
        ax.set_xlabel('Latitude (°)')

    # (c) Column density
    ax = axes[r, 2]
    ax.plot(lag_lat,  lag_cd[r],  color='#1F77B4', lw=1.8, ls='-', label='Lagranto')
    ax.plot(twod_lat, twod_cd[r], color='#C44E52', lw=1.8, ls='-', label='2-D model')
    ax.set_xlim(-90, 90); ax.set_xticks([-90,-60,-30,0,30,60,90])
    ax.grid(alpha=0.3)
    if r == 0:
        ax.set_title('(c) Column density', fontsize=11)
        ax.legend(fontsize=9, loc='upper right')
    if r == 2:
        ax.set_xlabel('Latitude (°)')
    ax.set_ylabel('column density', fontsize=9)
    info = (f'$B_{{\\rm LAG}}$={B_LAG_v[r]:.2f} Tg\n'
            f'$B_{{\\rm 2D}}$={B_2D_v[r]:.2f} Tg\n'
            f'ratio={ratio_v[r]:.2f}')
    ax.text(0.04, 0.96, info, transform=ax.transAxes,
            fontsize=8, va='top', ha='left',
            bbox=dict(facecolor='white', edgecolor='grey', alpha=0.85, pad=2.5))

fig.subplots_adjust(left=0.07, right=0.97, top=0.94, bottom=0.06,
                    hspace=0.30, wspace=0.40)
savefig_both(fig, HERE)
