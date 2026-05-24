#!/usr/bin/env python3
"""Fig 10/11 (fig:heating_materials) — silica + calcite, symmetric injection.

Layout: 2 rows x 2 cols
  Row 1 (a, b): SARF efficacy  epsilon       (W m^-2 per Tg yr^-1)   <= 0
  Row 2 (c, d): Heating cost   eta_heating   (K per W m^-2)          >= 0
  Col 1: Silica  d_0 = 0.5 um, with coagulation
  Col 2: Calcite d_0 = 0.3 um, with coagulation

Colormap (post ACP-audit pass):
  Row 1 (epsilon, all <= 0): COOL_HALF sequential -- dark navy = most cooling.
                              Both materials share the same scale.
  Row 2 (eta, all >= 0):     WARM_HALF sequential -- dark red  = most heating.
                              Silica and calcite use per-material scales because
                              their ranges differ by ~order of magnitude.

Outputs fig.pdf (vector, ACP submission) and fig.png (raster, preview).
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

ds = xr.open_dataset(os.path.join(HERE, "data.nc"))
inj_lat  = ds.inj_lat.values
inj_pres = ds.inj_pres.values
eps      = ds.eps.values
eta      = ds.eta_heating.values
tp_lat   = ds.tp_lat.values
tp_pres  = ds.tp_pres_clim.values
ds.close()

P_MAX = 120
ALT_TICKS = [15, 17, 19, 21]
TITLES = [r"Silica, $d_0 = 0.5\,\mu$m, with coag",
          r"Calcite, $d_0 = 0.3\,\mu$m, with coag"]


def p_to_z(p): return 7.0 * np.log(1000.0 / p)
def z_to_p(z): return 1000.0 * np.exp(-z / 7.0)


def fmt_axes(ax, ylabel=False, altaxis=False, xlabel=True):
    ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_ylim(P_MAX * 1.02, inj_pres.min() * 0.98)
    ax.set_xlim(-1, inj_lat.max() + 3)
    ax.set_xticks([0, 20, 40, 60])
    if xlabel:
        ax.set_xlabel("Injection latitude (°)")
    ax.yaxis.set_major_locator(FixedLocator([50, 70, 100, 120]))
    ax.yaxis.set_minor_locator(FixedLocator([60, 80, 90, 110]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v)}"))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.plot(tp_lat, tp_pres, "w-", lw=2.2)
    ax.plot(tp_lat, tp_pres, "k-", lw=1.0)
    if ylabel:
        ax.set_ylabel("Injection pressure (hPa)")
    if altaxis:
        ax2 = ax.secondary_yaxis("right", functions=(p_to_z, z_to_p))
        ax2.set_ylabel("Altitude (km)", labelpad=12)
        ax2.set_yticks(ALT_TICKS)
        ax2.set_yticklabels([str(v) for v in ALT_TICKS])
        ax2.minorticks_off()


# Sequential norms. ε vmax shifted to -0.05 (not 0) so the lightest in-scale
# band reads pale blue rather than near-white; extend='max' below caps the
# small-magnitude cells in (-0.05, 0] to off-scale white.
levels_eps    = np.arange(-0.40, -0.05 + 0.001, 0.05)
labels_eps    = [-0.35, -0.30, -0.25, -0.20, -0.15, -0.10]
norm_eps      = Normalize(vmin=-0.40, vmax=-0.05)

levels_eta_si = np.arange(0.0, 2.401, 0.10)
labels_eta_si = [0.4, 0.8, 1.2, 1.6, 2.0]
norm_eta_si   = Normalize(vmin=0.0, vmax=2.4)

# Calcite eta scale extended to 0.4 so the ~0.36 peak at the high-lat lower-strat
# corner is no longer saturated (audit pass fix).
levels_eta_ca = np.arange(0.0, 0.401, 0.025)
labels_eta_ca = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]
norm_eta_ca   = Normalize(vmin=0.0, vmax=0.40)

fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)

pms_eps = []
pms_eta = []
for c_idx in range(eps.shape[0]):
    eps_g = np.clip(eps[c_idx], levels_eps[0], levels_eps[-1])
    if c_idx == 0:    # silica
        norm_eta, levels_eta, labels_eta = norm_eta_si, levels_eta_si, labels_eta_si
    else:             # calcite
        norm_eta, levels_eta, labels_eta = norm_eta_ca, levels_eta_ca, labels_eta_ca
    eta_g = np.clip(eta[c_idx], levels_eta[0], levels_eta[-1])

    ax = axes[0, c_idx]
    pm = ax.contourf(inj_lat, inj_pres, eps_g,
                     levels=levels_eps, cmap=COOL_SEQ, norm=norm_eps, extend='max')
    cs = ax.contour(inj_lat, inj_pres, eps_g, levels=labels_eps,
                    colors="k", linewidths=0.6, alpha=0.7)
    ax.clabel(cs, inline=True, fmt="%.2f", fontsize=8)
    fmt_axes(ax, ylabel=(c_idx == 0), altaxis=(c_idx == 1), xlabel=False)
    ax.set_title(TITLES[c_idx], fontsize=11)
    panel_label(ax, ["a", "b"][c_idx])
    pms_eps.append(pm)

    ax = axes[1, c_idx]
    pm = ax.contourf(inj_lat, inj_pres, eta_g,
                     levels=levels_eta, cmap=WARM_SEQ, norm=norm_eta, extend="max")
    cs = ax.contour(inj_lat, inj_pres, eta_g, levels=labels_eta,
                    colors="k", linewidths=0.6, alpha=0.7)
    ax.clabel(cs, inline=True, fmt="%.2f", fontsize=8)
    fmt_axes(ax, ylabel=(c_idx == 0), altaxis=(c_idx == 1), xlabel=True)
    panel_label(ax, ["c", "d"][c_idx])
    pms_eta.append(pm)

fig.subplots_adjust(left=0.07, right=0.84, top=0.94, bottom=0.08,
                    hspace=0.22, wspace=0.45)

# Top-row shared epsilon colorbar (far right of top half)
cax_eps = fig.add_axes([0.92, 0.55, 0.012, 0.39])
cb_eps  = fig.colorbar(pms_eps[1], cax=cax_eps,
                       ticks=[-0.40, -0.30, -0.20, -0.10, -0.05])
cb_eps.set_label(r"$\epsilon$ (W m$^{-2}$ per Tg yr$^{-1}$)")

# Bottom-row eta: silica cbar in inter-panel gap (label on LEFT to clear panel d);
# calcite cbar at far right of bottom half (after altitude axis on panel d).
cax_eta_si = fig.add_axes([0.46, 0.08, 0.012, 0.39])
cb_eta_si  = fig.colorbar(pms_eta[0], cax=cax_eta_si,
                          ticks=[0.0, 0.6, 1.2, 1.8, 2.4])
cb_eta_si.ax.yaxis.set_label_position("left")
cb_eta_si.set_label(r"$\eta_{\rm heating}$, silica  (K W$^{-1}$ m$^{2}$)")

cax_eta_ca = fig.add_axes([0.92, 0.08, 0.012, 0.39])
cb_eta_ca  = fig.colorbar(pms_eta[1], cax=cax_eta_ca,
                          ticks=[0.0, 0.10, 0.20, 0.30, 0.40])
cb_eta_ca.set_label(r"$\eta_{\rm heating}$, calcite  (K W$^{-1}$ m$^{2}$)")

savefig_both(fig, HERE)
