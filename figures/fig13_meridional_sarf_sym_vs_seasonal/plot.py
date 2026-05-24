#!/usr/bin/env python3
"""Fig 12 (fig:meridional_sarf_profiles) — §3.4 Role of injection season."""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))
ds = xr.open_dataset(os.path.join(HERE, 'data.nc'))
lat = ds.lat.values
sym = ds.sym_profile.values
alt = ds.alt_profile.values
f_sym = ds.f_sym.values
f_alt = ds.f_alt.values
tau_sym = ds.tau_sym.values
tau_alt = ds.tau_alt.values
case_label = ds.case_label.values
ds.close()

n_cases = sym.shape[0]
fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True, sharey=True)
panels = axes.flatten()
tags = ['(a)', '(b)', '(c)', '(d)']

for ci in range(n_cases):
    ax = panels[ci]
    lab_sym = (r'symmetric  ($\tau=%.2f$ yr,  $f=%+.2f$ W m$^{-2}$/Tg)'
               % (tau_sym[ci], f_sym[ci]))
    lab_alt = (r'seasonal   ($\tau=%.2f$ yr,  $f=%+.2f$ W m$^{-2}$/Tg)'
               % (tau_alt[ci], f_alt[ci]))
    ax.plot(lat, sym[ci], color='#2E5984', lw=2.0, label=lab_sym)
    ax.plot(lat, alt[ci], color='#C44E52', lw=2.0, ls='--', label=lab_alt)
    ax.axhline(0,  color='grey', lw=0.6, ls=':')
    ax.axhline(-1, color='grey', lw=0.6, ls=':')
    ax.set_xlim(-90, 90)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    if ci >= 2:
        ax.set_xlabel('Latitude (degree)', fontsize=12)
    if ci % 2 == 0:
        ax.set_ylabel('normalised SARF anomaly\nper response lat (W m$^{-2}$)',
                      fontsize=11)
    ax.set_title(f'injection at {case_label[ci]}', fontsize=12)
    ax.text(0.02, 0.97, tags[ci], transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=2.0))
    ax.grid(alpha=0.3)
    ax.tick_params(axis='both', labelsize=10)
    ax.legend(loc='lower left', fontsize=9, framealpha=0.95)

fig.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.09,
                    wspace=0.10, hspace=0.25)
savefig_both(fig, HERE)
