# Solid-particle stratospheric aerosol injection — a 2-D modeling exploration of the design space

Companion repository for the paper *"Solid-particle stratospheric aerosol
injection: a 2-D modeling exploration of the design space"* (Lederer et al.,
submitted to *Atmospheric Chemistry and Physics*, 2026). This umbrella
repository is the entry point for reproducing the paper's results: it hosts
the paper-specific driver scripts and per-figure data + rendering code, and
links out to the simulation libraries and input-data engine in the component
repositories.

## Status

**`v0.2.x`.** All five component repositories are public, the input-data
deposit is published on Zenodo, and the SARF + 2-D transport pipelines in
`scripts/` run end-to-end against the public chain with no credentials. The
`v1.0` release will be tagged to match the paper's version of record.

## Paper

> Yoav Lederer, Nahliel Wygoda, Dorri Halbertal, and Brian E. J. Rose.
> **Solid-particle stratospheric aerosol injection: a 2-D modeling exploration
> of the design space.** Submitted to *Atmospheric Chemistry and Physics*,
> 2026.
>
> Preprint DOI: [`10.5194/egusphere-2026-2772`](https://doi.org/10.5194/egusphere-2026-2772)

Author affiliations:
- Yoav Lederer (corresponding) — Stardust Labs ltd, Ness Ziona, Israel
- Nahliel Wygoda — Stardust Labs ltd, Ness Ziona, Israel
- Dorri Halbertal — Stardust Labs ltd, Ness Ziona, Israel
- Brian E. J. Rose — University at Albany, State University of New York, USA

## Repository layout

```
.
├── scripts/
│   ├── build_sarf_ref.py         build the SARF reference state
│   ├── preprocess_profiles.py    transport runs → SARF-ready profiles
│   ├── run_sarf_case.py          compute SARF for one preprocessed case
│   ├── run_sarf_sweep.py         parallel SARF sweep over a manifest
│   └── run_transport.py          2-D aerosol-transport runner
├── figures/                      per-figure data + plot.py (17 figures)
│   └── README.md                 figure → folder → label → data-source map
├── README.md / CITATION.cff / LICENSE / CHANGELOG.md / …
└── component repositories (separate GitHub repos under stardust-initiative):
```

| Component | Role |
|---|---|
| [`stardust-climate`](https://github.com/stardust-initiative/stardust-climate) | Climate-model runner library — the model-factory, reference-model, post-processing primitives the pipeline scripts here import. |
| [`climlab-stardust-extension`](https://github.com/stardust-initiative/climlab-stardust-extension) | Extension of `climlab` providing the aerosol-layer radiation coupling, the modified convection, the 2-D meridional diffusion–advection transport, and aerosol optical-depth tables. |
| [`climlab-rrtmg_stardust`](https://github.com/stardust-initiative/climlab-rrtmg_stardust) | Stardust fork of `climlab`'s RRTMG bindings — OpenMP column parallelization, ensemble cloud sampling, spectral SW flux exposure, aerosol-layer (`r_mu`/`t_mu`/`r_bar`/`t_bar`) hooks. |
| [`climlab-sbm-convection_stardust`](https://github.com/stardust-initiative/climlab-sbm-convection_stardust) | Stardust fork of the Simplified Betts–Miller convection scheme, with the explicit surface-pressure variant (`betts_miller_pstar`) needed for the 2-D model's non-uniform pressure grid. |
| [`stardust-2d-inputs`](https://github.com/stardust-initiative/stardust-2d-inputs) | Input-data engine — lean runtime loader + provenance registry, Zenodo-backed for the public release set. |

## Reproducing the paper's results

The pipeline runs in a fresh conda environment with no credentials. Every
Stardust-modified dependency is pinned by version tag in the runner's
`pyproject.toml`, so a single `pip install` resolves the whole chain.

### 1. Set up the environment

```bash
conda create -n stardust_env python=3.11 -y
conda activate stardust_env
conda install -c conda-forge climlab compilers meson meson-python -y
pip install git+https://github.com/stardust-initiative/stardust-climate.git@v0.1.2
```

The Stardust forks of `climlab-rrtmg` and `climlab-sbm-convection` are
compiled from source during `pip install`; that's why the conda environment
includes a compiler toolchain. The `stardust-2d-inputs` engine is pulled in
transitively and ships with a built-in default config — the loader fetches
input data from the Zenodo deposit
([10.5281/zenodo.20271742](https://doi.org/10.5281/zenodo.20271742)) on
demand, content-addressed and SHA-256-verified.

### 2. Build the SARF reference state

```bash
cd scripts/
python build_sarf_ref.py --grid-npz <a-transport-run.npz> --out-dir sarf-ref/
```

The reference is an ERA5-pinned, transport-free radiative-convective column
model integrated to a radiative fixed point on the lat/lev grid extracted
from a 2-D transport run. See `--help` for the spin-up knobs.

### 3. Run the 2-D aerosol-transport simulation(s)

```bash
python run_transport.py \
    --injection-geometry symmetric --inject-lat 30 \
    --months 60 \
    --output transport_run.npz
```

This is the offline 2-D aerosol-transport calculation (residual circulation
+ eddy diffusion + sedimentation + coagulation across size bins). For
paper-scale runs use `--months 60` or longer; the default is a short proof
configuration.

### 4. Preprocess transport profiles for SARF

```bash
python preprocess_profiles.py --mode annual --input-root <dir-of-transport-runs>
```

Reduces the per-injection-point transport `.npz` outputs into the flat
per-bin mass-mixing-ratio NetCDFs that the SARF driver consumes.

### 5. Compute SARF — one case or a parallel sweep

One case:

```bash
python run_sarf_case.py \
    --preprocessed-nc <file.nc> \
    --ref-nc sarf-ref/model_ref_era5ref_minimal_notransp.nc \
    --out-nc result.nc
```

Parallel sweep over a preprocessed-profile manifest:

```bash
python run_sarf_sweep.py --pre-root <preprocessed-tree> --ref-dir sarf-ref/
```

### 6. Regenerate the figures

```bash
cd figures/figNN_<name>
python plot.py
```

Every script-regenerable manuscript figure is in its own folder under
`figures/`, with the data it plots (`data.nc` or local `*.npz`) and a small
`plot.py` that reads only that local data. See
[`figures/README.md`](figures/README.md) for the full figure ↔ folder map.

## Scope of what's reproducible from this repository

This repository hosts the **paper-specific driver scripts and per-figure
data**. Some figures cannot be reproduced from a small data file:

- The schematic of the modelling scheme (paper Fig 1, panel a — the diagram)
  is produced from a slide deck; only the validation-AoA half of Fig 1 is
  script-regenerable and ships in `figures/fig01_validation_aoa/`.

What is reproducible end-to-end from this repository:

- The 2-D aerosol-transport simulations underlying §§2–3.
- The SARF computations underlying §§4–5.
- Every machine-generated figure in the paper (main text and appendix), from
  the cached per-figure `data.nc` via the local `plot.py`.

## Data availability

The reference atmospheric state files (ERA5-derived zonal-mean climatologies),
the 2-D transport drivers (TEM-decomposed residual circulation + eddy
diffusivity + tropopause climatology), and the RRTMG-banded aerosol
optical-property tables ship via the public **Zenodo deposit**:

> **Stardust climate-model input database, v0.1.0.**
> DOI: [`10.5281/zenodo.20271742`](https://doi.org/10.5281/zenodo.20271742)
> (concept DOI: [`10.5281/zenodo.20271741`](https://doi.org/10.5281/zenodo.20271741))

The deposit's README + `deposit_manifest.json` list every file with its
SHA-256 hash, upstream source, and the generator that produced it. The
`stardust-2d-inputs` engine pulls from this deposit on demand and verifies
each file against its hash.

## Contact

- **Scientific / paper correspondence** — corresponding author **Yoav
  Lederer**, [y.lederer@stardust-initiative.com](mailto:y.lederer@stardust-initiative.com)
  (Stardust Labs ltd, Ness Ziona, Israel).
- **Code / repository correspondence** — repository maintainer **Dorri
  Halbertal**, [d.halbertal@stardust-initiative.com](mailto:d.halbertal@stardust-initiative.com)
  (Stardust Labs ltd).

For routine bug reports, documentation problems, or reproducibility
questions, please open a GitHub issue on this repository rather than
emailing directly. For matters that should be reported privately, see
[`SECURITY.md`](SECURITY.md).

## License

This umbrella repository's content (driver scripts, figure code + data,
documentation, citation metadata) is released under the MIT License — see
[`LICENSE`](LICENSE). The individual component repositories carry their own
licenses; please refer to each component for details.

## Citation

If you use any of the materials linked from this repository, please cite the
paper (see [`CITATION.cff`](CITATION.cff) or use the "Cite this repository"
button on GitHub) and, where appropriate, the Zenodo deposit. Final citation
details (ACP DOI, paper volume/page numbers) will land here on acceptance.
