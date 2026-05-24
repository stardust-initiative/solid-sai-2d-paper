# %%
import sys
import os
sys.path.append(os.getcwd())
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import interpolate
from scipy.interpolate import interp1d

data5 = np.load("for_paper_2026_04_age_of_air_10years.npz")
lat = data5["lat"]
lev = data5["lev"]
aoa_2d = data5["age_of_air"]
lagranto = np.load("Lagranto_age_of_air_try3_years_10_20.npz")
lagranto_age = lagranto['mean_age'].copy()
lagranto_age[np.isnan(lagranto_age)] = 0
interp_lagranto = interpolate.RectBivariateSpline(lagranto['lat_centers'],lagranto['p_centers'],lagranto_age, kx=1, ky=1)
lagranto1 = interp_lagranto(lat,lev)

interp_2d = interpolate.RectBivariateSpline(lat,lev,aoa_2d, kx=1, ky=1)
aoa_at_50hpa_2d = interp_2d(lat,50.0)
interp_lag = interpolate.RectBivariateSpline(lat,lev,lagranto1, kx=1, ky=1)
aoa_at_50hpa_lag = interp_lag(lat,50.0)

lat_range = np.array([-10, 10])
lat_mask = (lat >= lat_range[0]) & (lat <= lat_range[1])
mean_age_trop_2d = (
    (aoa_2d[lat_mask, :] * np.cos(lat[lat_mask, None] * np.pi / 180)).sum(axis=0)
    / np.cos(lat[lat_mask] * np.pi / 180).sum()
)
mean_age_trop_lag = (
    (lagranto1[lat_mask, :] * np.cos(lat[lat_mask, None] * np.pi / 180)).sum(axis=0)
    / np.cos(lat[lat_mask] * np.pi / 180).sum()
)

SMALL_SIZE = 16
MEDIUM_SIZE = 16

def _make_legend(ax, **kw):
    leg = ax.legend(**kw)
    leg.get_frame().set_edgecolor('black')
    leg.get_frame().set_linewidth(2.0)
    return leg

def _add_legend_separator(fig, leg, sep_after=4):
    """Draw a horizontal line inside the legend box between entries sep_after and sep_after+1."""
    renderer = fig.canvas.get_renderer()
    texts = leg.get_texts()
    bb_above = texts[sep_after].get_window_extent(renderer)
    bb_below = texts[sep_after + 1].get_window_extent(renderer)
    y = (bb_above.y0 + bb_below.y1) / 2
    frame = leg.get_frame().get_window_extent(renderer)
    pad = 6
    inv = fig.transFigure.inverted()
    p0 = inv.transform([frame.x0 + pad, y])
    p1 = inv.transform([frame.x1 - pad, y])
    fig.add_artist(Line2D([p0[0], p1[0]], [p0[1], p1[1]],
                          transform=fig.transFigure, color='gray',
                          linewidth=2.0, clip_on=False))

# --- Load digitized reference data from Chabrillat et al. (2018) Fig. 4a ---
# CSV columns: Obs(X,Y), error_low(X,Y), error_high(X,Y), CFSR(X,Y),
#              JRA(X,Y), MERRA2(X,Y), ERAI(X,Y)  — two rows of headers, then data.
def _extract(raw, col):
    x = pd.to_numeric(raw.iloc[:, col],   errors='coerce')
    y = pd.to_numeric(raw.iloc[:, col+1], errors='coerce')
    mask = x.notna() & y.notna()
    x, y = x[mask].values, y[mask].values
    order = np.argsort(x)
    return x[order], y[order]

_raw = pd.read_csv("chabrillat_18_fig_5_27_a_datasets.csv", header=None, skiprows=2)
obs_lat,      obs_age      = _extract(_raw,  0)
err_low_lat,  err_low_age  = _extract(_raw,  2)
err_high_lat, err_high_age = _extract(_raw,  4)
cfsr_lat,     cfsr_age     = _extract(_raw,  6)
jra_lat,      jra_age      = _extract(_raw,  8)
merra2_lat,   merra2_age   = _extract(_raw, 10)
erai_lat,     erai_age     = _extract(_raw, 12)

# Interpolate error bounds onto a shared fine grid for fill_between
lat_fill = np.linspace(max(err_low_lat[0], err_high_lat[0]),
                       min(err_low_lat[-1], err_high_lat[-1]), 500)
err_low_fill  = interp1d(err_low_lat,  err_low_age,  kind='linear')(lat_fill)
err_high_fill = interp1d(err_high_lat, err_high_age, kind='linear')(lat_fill)

# --- Combined figure: panel (a) left, panel (b) right ---
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(18, 7))

# --- Panel (a): mean age at 50 hPa vs latitude ---
ax.fill_between(lat_fill, err_low_fill, err_high_fill,
                color='grey', alpha=0.35, label='_nolegend_')
ax.plot(obs_lat,    obs_age,    color='black',   marker='D', markersize=9, linestyle='none', label='Obs')
ax.plot(cfsr_lat,  cfsr_age,  color='#0072B2', linestyle='-',  linewidth=4, label='CFSR-CFSv2')
ax.plot(jra_lat,   jra_age,   color='#D55E00', linestyle='--', linewidth=4, label='JRA-55')
ax.plot(merra2_lat, merra2_age, color='#009E73', linestyle=':',  linewidth=4, label='MERRA-2')
ax.plot(erai_lat,  erai_age,  color='#CC79A7', linestyle='-.', linewidth=4, label='ERA-I')
ax.plot(lat, aoa_at_50hpa_2d,  color='#E69F00', linestyle='-',  linewidth=4, label='2D model')
ax.plot(lat, aoa_at_50hpa_lag, color='#56B4E9', linestyle='--', linewidth=4, label='Lagranto')

ax.set_xlim([-90, 90])
ax.set_ylim([0, 6])
ax.set_xlabel('Latitude', fontsize=MEDIUM_SIZE)
ax.set_ylabel('Age of Air [years]', fontsize=MEDIUM_SIZE)
ax.tick_params(axis='both', labelsize=SMALL_SIZE)
ax.set_xticks(np.arange(-90, 100, 30))
ax.grid(True, alpha=0.3, linestyle='--')
leg_a = _make_legend(ax, fontsize=SMALL_SIZE, framealpha=0.9, loc='lower right')
ax.text(0.02, 0.98, '(a)', transform=ax.transAxes,
        fontsize=MEDIUM_SIZE+6, fontweight='bold', va='top', ha='left')

# --- Panel (b): tropical AoA vertical profile ---
# X = AoA (years), Y = pressure (hPa), log-inverted y-axis.

_raw_b = pd.read_csv("chabrillat_18_fig_5_27_b_datasets.csv", header=None, skiprows=2)

def _extract_b(raw, col):
    """Extract (aoa, pressure) pair, sorted by ascending pressure."""
    x = pd.to_numeric(raw.iloc[:, col],   errors='coerce')
    y = pd.to_numeric(raw.iloc[:, col+1], errors='coerce')
    mask = x.notna() & y.notna()
    x, y = x[mask].values, y[mask].values
    order = np.argsort(y)
    return x[order], y[order]

obs_aoa_b,       obs_pres_b       = _extract_b(_raw_b,  0)
err_low_aoa_b,   err_low_pres_b   = _extract_b(_raw_b,  2)
err_high_aoa_b,  err_high_pres_b  = _extract_b(_raw_b,  4)
cfsr_aoa_b,      cfsr_pres_b      = _extract_b(_raw_b,  6)
jra_aoa_b,       jra_pres_b       = _extract_b(_raw_b,  8)
merra2_aoa_b,    merra2_pres_b    = _extract_b(_raw_b, 10)
erai_aoa_b,      erai_pres_b      = _extract_b(_raw_b, 12)

# Per-point rectangle bounds: geometric midpoints between consecutive pressure levels.
# obs, error_low, error_high share the same pressure levels (same count, same order).
_p = obs_pres_b  # sorted ascending (low pressure first = top of inverted axis)
_mid = np.sqrt(_p[:-1] * _p[1:])
_y_lo = np.empty(len(_p))
_y_hi = np.empty(len(_p))
_y_lo[0]   = _p[0]  * (_p[0]  / _p[1])  ** 0.5   # extend upward from first point
_y_lo[1:]  = _mid
_y_hi[:-1] = _mid
_y_hi[-1]  = _p[-1] * (_p[-1] / _p[-2]) ** 0.5   # extend downward from last point
_y_lo = np.clip(_y_lo, 5, 100)
_y_hi = np.clip(_y_hi, 5, 100)

# Model data filtered to plot range
_pm = (lev >= 5) & (lev <= 110)

# Grey rectangles: one per obs point, spanning error_low to error_high in AoA
for i in range(len(_p)):
    ax2.fill_betweenx([_y_lo[i], _y_hi[i]], err_low_aoa_b[i], err_high_aoa_b[i],
                      color='grey', alpha=0.35, linewidth=0)

ax2.plot(obs_aoa_b,    obs_pres_b,    color='black',   marker='D', markersize=9, linestyle='none', label='Obs')
ax2.plot(cfsr_aoa_b,   cfsr_pres_b,   color='#0072B2', linestyle='-',  linewidth=4, label='CFSR-CFSv2')
ax2.plot(jra_aoa_b,    jra_pres_b,    color='#D55E00', linestyle='--', linewidth=4, label='JRA-55')
ax2.plot(merra2_aoa_b, merra2_pres_b, color='#009E73', linestyle=':',  linewidth=4, label='MERRA-2')
ax2.plot(erai_aoa_b,   erai_pres_b,   color='#CC79A7', linestyle='-.', linewidth=4, label='ERA-I')
ax2.plot(mean_age_trop_2d[_pm],  lev[_pm], color='#E69F00', linestyle='-',  linewidth=4, label='2D model')
ax2.plot(mean_age_trop_lag[_pm], lev[_pm], color='#56B4E9', linestyle='--', linewidth=4, label='Lagranto')

ax2.set_yscale('log')
ax2.invert_yaxis()
ax2.set_ylim([100, 5])
ax2.set_xlim([0, 6])

# All 15 values become major ticks → gridlines at every one of them.
# Only a subset gets a visible label to avoid crowding.
_all_yticks  = [5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
_label_set   = {5, 10, 20, 30, 50, 70, 100}
ax2.set_yticks(_all_yticks)
ax2.set_yticklabels([str(v) if v in _label_set else '' for v in _all_yticks])

ax2.set_xlabel('Age of Air [years]', fontsize=MEDIUM_SIZE)
ax2.set_ylabel('Pressure [hPa]', fontsize=MEDIUM_SIZE)
ax2.tick_params(axis='both', labelsize=SMALL_SIZE)
ax2.grid(True, alpha=0.3, linestyle='--')
leg_b = _make_legend(ax2, fontsize=SMALL_SIZE, framealpha=0.9, loc='lower right')
ax2.text(0.02, 0.98, '(b)', transform=ax2.transAxes,
         fontsize=MEDIUM_SIZE+6, fontweight='bold', va='top', ha='left')

plt.tight_layout()
fig.canvas.draw()
_add_legend_separator(fig, leg_a)
_add_legend_separator(fig, leg_b)
plt.savefig('for_paper_2026_04_aoa_combined.png', dpi=300, bbox_inches='tight')
plt.show()

# %%
