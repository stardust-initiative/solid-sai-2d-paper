# Radiation-comparison figure — self-contained bundle

This bundle lets you regenerate and modify the radiation-transfer comparison
figure (MODTRAN vs original / modified RRTMG, SW radiative forcing per unit
burden vs latitude) without any of the heavy machinery — no climlab, no
MODTRAN, no CESM. Everything the plotting script needs is included as a small
set of cached NumPy archives.

## Contents

```
radiation_figure_bundle/
├── README.md                          (this file)
├── requirements.txt                   (numpy + matplotlib, that's it)
├── scripts/
│   └── plot_appendix_modtran_rrtmg.py (the plotting script)
├── data/
│   └── silica/
│       └── rf_results_day091.npz      (the day shown in the paper)
└── figures/                           (output directory; created on first run)
```

> **Filename note:** `plot_appendix_modtran_rrtmg.py` — the "appendix" is
> historical. This figure originated in the supplementary material and was
> later promoted to the main text. The name is kept to match the source
> repository so cross-references stay unambiguous.

## Quick start

```bash
# 1. (recommended) make a clean Python environment
python -m venv venv && source venv/bin/activate     # or use conda
pip install -r requirements.txt

# 2. generate the paper figure (day 91, silica)
python scripts/plot_appendix_modtran_rrtmg.py

# 3. the figure is written to figures/appendix_modtran_rrtmg_day091.png
```

That's it. The script reads `data/silica/rf_results_day091.npz` and writes a
300-dpi PNG into `figures/`.

## Options

```
python scripts/plot_appendix_modtran_rrtmg.py [--day 91] [--material silica]
                                              [--ymin -0.85] [--ymax 0.05]
                                              [--outdir DIR]
```

- `--day` — any of the cached days: 21, 52, 81, 91, 111, 141, 172, 202, 233,
  264, 294, 326, 356. The figure used in the paper is **day 91** (a near-equinox
  date, chosen so the whole latitude range is illuminated). The qualitative
  picture is the same at the other days.
- `--material` — `silica` is the only material included in this bundle.
- `--ymin` / `--ymax` — y-axis limits, in W m⁻² Tg⁻¹.
- `--outdir` — write the PNG somewhere else (default `figures/`).

## Modifying the figure

The script is short, plain matplotlib, and meant to be edited. The key spot is
`build_figure()`:

- **Curves plotted.** Four are drawn by default: MODTRAN (black solid),
  original RRTMG "3D" (green solid), original RRTMG "2D" (green dashed), and
  modified RRTMG "2D" (purple dashed). A fifth — modified RRTMG "3D" (purple
  solid) — is computed (`rf_new_inst`) but commented out, because it sits almost
  exactly on top of the MODTRAN curve and adds clutter. Un-comment that line if
  you want it back.
- **What "3D" and "2D" mean.** "3D" = the radiation code is called at every hour
  of the day and the result is integrated diurnally — the way RRTMG is naturally
  driven inside a 3D GCM. "2D" = the radiation code is called once per latitude
  with the daily-averaged cosine of the solar zenith angle as the only
  solar-geometry input — the way a zonally-averaged 2D transport model is forced
  to drive it. The modification to RRTMG closes most of the gap between the "2D"
  and "3D" behaviours.
- **Labels, colours, legend position, axes** — all in `build_figure()`. Change
  freely.
- The figure refreshes in well under a second, so iterate quickly: edit, rerun,
  look at `figures/`.

## What's in each `rf_results_day{DDD}.npz`

Each archive stores RF fields on a `(lat × hour)` grid (extended symmetrically
around local noon) plus daily-averaged summary vectors. Keys used by the
plotting script:

| key | shape | meaning |
|---|---|---|
| `lat_vect` | `(91,)` | latitudes [deg] |
| `hours_mat_extended` | `(91, 481)` | valid-hour mask per latitude (NaN outside daylight) |
| `rf_mat_modtran_extended` | `(91, 481)` | MODTRAN layer RF at (lat, hour) [W m⁻²] |
| `model_rf_mat_old_rrtmg_extended` | `(91, 481)` | original RRTMG, hourly, at (lat, hour) |
| `model_rf_mat_extended` | `(91, 481)` | modified RRTMG, hourly, at (lat, hour) |
| `model_rf_mat_old_rrtmg_daily_avg` | `(91,)` | original RRTMG, single call with ⟨μ⟩_daily |
| `model_rf_mat_new_rrtmg_daily_avg` | `(91,)` | modified RRTMG, single call with ⟨μ⟩_daily |

(Other keys are present — `coszen_extended`, `coszen_daily_avg_vect`, etc. — and
are not needed by this script.) Divide any RF field by `M_TG = 10` to get
per-Tg units; the script does this.

## Provenance / regenerating the data from scratch

The `.npz` archives are cached outputs of `comparison_to_modtran.py`, which:
1. loads a reference climate state (climlab),
2. builds a silica aerosol layer and runs the semi-analytic layer model,
3. constructs original- and modified-RRTMG radiation columns and computes ASR,
4. reads precomputed MODTRAN TAPE6 sweeps over (latitude × cos-zenith),
5. projects everything back onto the (lat, hour) grid and writes the `.npz`.

That full pipeline (with climlab, the MODTRAN run outputs, and the modified
RRTMG fork) lives in the transport-paper fork repository. You do **not** need
it to regenerate or modify the figure — the cached `.npz` archives are
sufficient. Ask the sender if you ever want the full regeneration pipeline.

## Size note

The 13 `.npz` files total ~26 MB. If you only need the paper figure (day 91),
you can delete all of `data/silica/rf_results_day*.npz` except
`rf_results_day091.npz` (~2 MB).
