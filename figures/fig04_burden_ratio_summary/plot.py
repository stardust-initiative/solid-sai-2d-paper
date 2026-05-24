#!/usr/bin/env python3
"""Fig 3 (fig:base_fn_burden_ratio) — §3.1 base distribution functions."""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _style import apply, savefig_both
apply()

HERE = os.path.dirname(os.path.abspath(__file__))
ds   = xr.open_dataset(os.path.join(HERE, 'data.nc'))
inj_lat  = ds.inj_lat.values
inj_pres = ds.inj_pres.values
ratio    = ds.ratio.values
ds.close()

PRES_GROUPS = [
    ( 40, 'o', '#C44E52', '40 hPa injection'),
    ( 65, 's', '#E59C2D', '65 hPa injection'),
    ( 95, '^', '#984EA3', '95 hPa injection'),
    (120, 'D', '#3673BB', '120 hPa injection'),
]

fig, ax = plt.subplots(1, 1, figsize=(11, 5.6))
for pres, marker, colour, label in PRES_GROUPS:
    mask = np.isclose(inj_pres, pres)
    if not mask.any():
        continue
    ax.scatter(inj_lat[mask], ratio[mask],
               marker=marker, c=colour, s=80,
               edgecolor='k', linewidth=0.6, label=label)

ax.axhline(1.0, color='k', linestyle=':', alpha=0.5)
ax.set_xlim(-70, 70)
ax.set_xticks([-60, -45, -30, -15, 0, 15, 30, 45, 60])
ax.set_xlabel('Injection latitude (°)')
ax.set_ylabel('Global burden ratio  (2-D / LAGRANTO)')
ax.set_ylim(0.30, 1.40)
ax.grid(alpha=0.3)
ax.legend(title='injection altitude', loc='lower right',
          fontsize=10, title_fontsize=10, framealpha=0.95)

fig.tight_layout()
savefig_both(fig, HERE)
savefig_both(fig, HERE)
