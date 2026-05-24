#!/usr/bin/env python3
"""Fig 2 (fig:weisenstein_comparison) — coagulation model validation."""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))
ds   = xr.open_dataset(os.path.join(HERE, 'data.nc'))
N    = ds.N.values
m1m  = ds.mass_1tg_model.values
m8m  = ds.mass_8tg_model.values
w1   = ds.mass_1tg_w15.values
w8   = ds.mass_8tg_w15.values
ds.close()

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(N, m1m, 'o-', color='#1F77B4', mfc='#1F77B4',
        markersize=9, linewidth=2.2, label='this work, 1 Tg yr$^{-1}$')
ax.plot(N, m8m, 's-', color='#C44E52', mfc='#C44E52',
        markersize=9, linewidth=2.2, label='this work, 8 Tg yr$^{-1}$')
ax.plot(N, w1, 'o--', color='#1F77B4', mfc='white',
        markersize=9, linewidth=1.8, mew=1.6,
        label='Weisenstein et al. 2015, 1 Tg yr$^{-1}$')
ax.plot(N, w8, 's--', color='#C44E52', mfc='white',
        markersize=9, linewidth=1.8, mew=1.6,
        label='Weisenstein et al. 2015, 8 Tg yr$^{-1}$')

ax.set_xscale('log', base=2)
ax.set_xticks(list(N))
ax.set_xticklabels([str(int(v)) for v in N])
ax.set_xlim([0.72, 22])
ax.set_ylim([0.0, 1.0])
ax.tick_params(axis='both', labelsize=14)
ax.set_xlabel('No. of Monomers in Agglomerate', fontsize=18)
ax.set_ylabel('Mass Fraction', fontsize=18)
ax.grid(True, alpha=0.4, linestyle='--')
ax.legend(fontsize=12, loc='upper right', framealpha=0.95)

savefig_both(fig, HERE)
