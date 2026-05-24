# `figures/` — figure data + rendering scripts

For every machine-generated figure in the paper, this directory holds:

- the **data** that figure plots (`data.nc`, `*.npz`, or a small `data/` subdir), and
- the **rendering script** (`plot.py`) that turns it into the panel in the paper.

The rendered figures themselves (`.pdf` / `.png`) are not redistributed here —
they're in the paper. The upstream production code that *generated* the
data files lives in this repository's `scripts/`, in the component repos
([`stardust-climate`](https://github.com/stardust-initiative/stardust-climate),
[`climlab-stardust-extension`](https://github.com/stardust-initiative/climlab-stardust-extension),
[`stardust-2d-inputs`](https://github.com/stardust-initiative/stardust-2d-inputs)),
and in the Zenodo deposit
([`10.5281/zenodo.20271742`](https://doi.org/10.5281/zenodo.20271742)).

## Layout

Folder numbering matches the paper's figure numbering one-to-one (main text
**Figs 1–13**, appendix **Figs A1–A4**):

```
figures/
├── README.md                                ← this file
├── _style.py                                ← shared matplotlib style + cmaps
├── fig01_validation_aoa/                    ← Fig 1
├── fig02_weisenstein/                       ← Fig 2
├── fig03_compare_to_modtran/                ← Fig 3
├── fig04_burden_ratio_summary/              ← Fig 4
├── fig05_three_examples_lag_vs_2D/          ← Fig 5
├── fig06_lifetime_05_coag_nocoag_deff/      ← Fig 6
├── fig07_ratio_Df16_Df30_lifetime/          ← Fig 7
├── fig08_ratio_03_07_lifetime_and_deff_coag/← Fig 8
├── fig09_swlw_split/                        ← Fig 9
├── fig10_silica_d0_optimum/                 ← Fig 10
├── fig11_heating_silica05_sym/              ← Fig 11
├── fig12_heating_silica05_seasonal/         ← Fig 12
├── fig13_meridional_sarf_sym_vs_seasonal/   ← Fig 13
├── figA1_kphiphi_gsfc_vs_era5/              ← Fig A1
├── figA2_seasonal_psi_res_tem/              ← Fig A2
├── figA3_material_panels/                   ← Fig A3
└── figA4_morphology_gallery/                ← Fig A4
```

A few figures in the paper are intentionally *not* here because they're not
script-regenerable from a small data file:

- The **Fig 1 schematic** (`figure1_diagram.png`) — produced from a slide
  deck; no `data.nc` to extract. The validation panel that the script
  composites alongside it is in `fig01_validation_aoa/`.

## How to regenerate a figure

```bash
cd figures/figNN_<name>
python plot.py
```

`plot.py` is fast (seconds) and reads only its local data file(s). It
expects `_style.py` to be one directory up and adds that to `sys.path`
itself — no installation needed.

### Dependencies for the plotting scripts

The standard scientific-Python stack:

```bash
pip install numpy matplotlib xarray
```

A few figures need extras and document them in their own folder (e.g.
`fig03_compare_to_modtran/requirements.txt`).

## Per-figure catalog

| Paper Fig | Folder | LaTeX label | Production data source |
|---|---|---|---|
| **Fig 1** | `fig01_validation_aoa` | `fig:validation_aoa` | Stardust 2-D transport AoA vs Chabrillat 2018 Fig 5.27 a/b + Lagranto reference (npz + CSV inputs; bespoke composer script) |
| **Fig 2** | `fig02_weisenstein` | `fig:weisenstein_comparison` | Stardust coagulation model + Weisenstein 2015 Fig 4c |
| **Fig 3** | `fig03_compare_to_modtran` | `fig:compare_to_modtran` | climlab + MODTRAN day-091 SW radiative-forcing comparison (silica) |
| **Fig 4** | `fig04_burden_ratio_summary` | `fig:base_fn_burden_ratio` | `results/compare_real2010/` Lagranto-vs-2D burden-ratio scatter |
| **Fig 5** | `fig05_three_examples_lag_vs_2D` | `fig:base_fn_examples` | `compare_runs_seasonal_for_paper.nc`, 2010 tropopause |
| **Fig 6** | `fig06_lifetime_05_coag_nocoag_deff` | `fig:lifetime_coag_nocoag_deff` | Silica d₀=0.5 µm transport runs + tropopause climatology (2008–2017) |
| **Fig 7** | `fig07_ratio_Df16_Df30_lifetime` | `fig:ratio_Df16_Df30_lifetime` | Silica coag runs at fractal dimensions Df=1.6 vs 3.0 |
| **Fig 8** | `fig08_ratio_03_07_lifetime_and_deff_coag` | `fig:ratio_03_07_lifetime_coag` | Silica d₀=0.3 vs 0.7 µm coag runs |
| **Fig 9** | `fig09_swlw_split` | `fig:swlw_split` | Seasonally-aggregated SARF — SW/LW split per material/size |
| **Fig 10** | `fig10_silica_d0_optimum` | `fig:silica_d0_optimum` | Yearly-aggregated SARF — silica d₀ optimum (nocoag vs coag) |
| **Fig 11** | `fig11_heating_silica05_sym` | `fig:heating_materials` | SARF maps — silica d₀=0.5 µm symmetric injection |
| **Fig 12** | `fig12_heating_silica05_seasonal` | `fig:lifetime_seasonal_ratio` | SARF maps — silica d₀=0.5 µm seasonal injection |
| **Fig 13** | `fig13_meridional_sarf_sym_vs_seasonal` | `fig:meridional_sarf_profiles` | Meridional SARF profiles — symmetric vs seasonal injection |
| **Fig A1** | `figA1_kphiphi_gsfc_vs_era5` | `fig:seasonal_Dphiphi` | GSFC-derived vs ERA5-derived horizontal-diffusivity Kφφ |
| **Fig A2** | `figA2_seasonal_psi_res_tem` | `fig:seasonal_Psi_resTEM` | TEM residual stream function Ψ_res,TEM — seasonal climatology |
| **Fig A3** | `figA3_material_panels` | `fig:mie_compare_dsweep` | Mie scattering, d₀=0.5 µm, silica + calcite |
| **Fig A4** | `figA4_morphology_gallery` | `fig:dda_all_morphologies` | Mie morphology gallery at R=250 nm (silica) |

## Conventions in the data files

- **`data.nc`** carries only the points the figure plots — typically tens of
  kilobytes per figure. Variable names mirror the figure's legend
  (`mass_1tg_model`, `eps`, `eta_heating`, `ratio_lifetime`, …); coordinate
  names mirror the physical axes (`inj_lat`, `inj_pres`, `lat`, `pres`, …).
- Each variable carries a human-readable `description` attribute and a
  `units` attribute where physically meaningful.
- Datasets carry `title`, `figure_label` (mirroring the LaTeX
  `\label{fig:...}`), `section`, and `source_files` global attributes.
- A handful of `figure_label` attrs are out of sync with the final
  manuscript (drafted earlier in the project). The folder name and the
  table above are the authoritative mapping; the LaTeX label column is the
  paper's `\label{fig:...}`.
