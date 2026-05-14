# Solid-particle stratospheric aerosol injection — a 2-D modeling exploration of the design space

Companion repository for the paper *"Solid-particle stratospheric aerosol
injection: a 2-D modeling exploration of the design space"* (Lederer et al.,
submitted to *Atmospheric Chemistry and Physics*, 2026). This umbrella
repository is the entry point for reproducing the paper's results: it links to
the simulation code, configuration files, and the paper-specific scripts
that drive the component-repo simulations and post-analyze their output
into the results discussed in the manuscript.

## Status

**Submission-state skeleton (v0.1).** The paper is currently under review at
*Atmospheric Chemistry and Physics*. This repository exists at submission time
to (a) provide a stable public URL we can cite in the paper, and (b) document
the structure of the planned full release. The component repositories listed
below are presently private; we are committed to making the full set publicly
available — together with a Zenodo-archived snapshot of inputs and pre-computed
outputs — before the paper is published in its final form.

The repository will be updated in place as the component repositories are
opened and the reproduction workflow is finalized, and a `v1.0` release will be
tagged to match the paper's version of record.

## Paper

> Yoav Lederer, Nahliel Wygoda, Dorri Halbertal, and Brian E. J. Rose.
> **Solid-particle stratospheric aerosol injection: a 2-D modeling exploration
> of the design space.** Submitted to *Atmospheric Chemistry and Physics*,
> 2026.
>
> DOI: <!-- TODO: fill in ACPD discussion-paper DOI when assigned, then the
> final ACP DOI on acceptance -->

Author affiliations:
- Yoav Lederer (corresponding) — Stardust-Labs, Ness Ziona, Israel
- Nahliel Wygoda — Stardust-Labs, Ness Ziona, Israel
- Dorri Halbertal — Stardust-Labs, Ness Ziona, Israel
- Brian E. J. Rose — University at Albany, State University of New York, USA

## Repository structure

In its current skeleton state (`v0.x`) this repository contains only
documentation and citation metadata; by `v1.0` (paper acceptance) it will
additionally host the paper-specific scripts — those that drive the
component-repo simulations and those that post-analyze their output into
the results discussed in the manuscript. The *generic* simulation machinery
it depends on lives in separate **component repositories** — all hosted (or
to be hosted) under the
[`stardust-initiative`](https://github.com/stardust-initiative) GitHub
organisation:

| Component | Role |
|---|---|
| [`stardust-climate`](https://github.com/stardust-initiative/stardust-climate) | Main analysis driver: SARF / ERF sweeps, post-processing pipeline, paper plot generation, plus the entry-point scripts for the offline 2-D aerosol-transport simulations. |
| [`climlab-stardust-extension`](https://github.com/stardust-initiative/climlab-stardust-extension) | Extension of `climlab` providing the aerosol-layer radiation coupling, the modified convection, and the 2-D meridional diffusion–advection transport used both for atmospheric moisture in the reference state and for the offline aerosol-transport simulations that produce the SARF inputs. |
| [`climlab-rrtmg_stardust`](https://github.com/stardust-initiative/climlab-rrtmg_stardust) | Fork of `climlab`'s RRTMG bindings with spectral diagnostics and the aerosol-layer (`r_mu` / `t_mu` / `r_bar` / `t_bar`) hooks. |
| [`climlab-sbm-convection_stardust`](https://github.com/stardust-initiative/climlab-sbm-convection_stardust) | Fork of `climlab`'s Simplified Betts–Miller convection scheme used by the reference state. |
| [`optical-tables-generator`](https://github.com/stardust-initiative/optical-tables-generator) | Refractive-index data, Mie code, and per-RRTMG-band optical tables for the candidate aerosol materials. <!-- TODO: repo does not yet exist; placeholder name reserved. --> |

The offline 2-D aerosol-transport machinery currently lives on a development
branch of a `climlab_stardust` working repository maintained by one of the
co-authors; for the published version it will be folded into
[`climlab-stardust-extension`](https://github.com/stardust-initiative/climlab-stardust-extension)
(library / numerics) and [`stardust-climate`](https://github.com/stardust-initiative/stardust-climate)
(entry-point scripts and configurations), so there will be no separate
"transport" repository.

**Note on availability.** All component repositories above are *private at
submission time*. They will be opened — and any final renaming applied — before
the paper is published in its final form. Links above resolve once each
component is made public; before that, they will return 404 to external
visitors. The corresponding author can grant pre-publication access to
reviewers on request.

## Reproducing the results

The intended end-to-end workflow, from a clean clone to the paper's figures:

1. **Clone** this umbrella repository together with the component repositories
   listed above. <!-- TODO: provide a single `git clone --recurse-submodules`
   instruction once the component URLs are fixed (and we have decided whether
   to use git submodules or a top-level setup script). -->
2. **Install** the modified `climlab` stack (`climlab-stardust-extension` and
   its `climlab-rrtmg-stardust` / `climlab-sbm-convection-stardust`
   dependencies) into a fresh Python environment. <!-- TODO: pin the Python
   version, list the conda/pip command, and document the (small number of)
   compiled-extension prerequisites. -->
3. **Download** the reference climate-database files and the material optical
   tables from the Zenodo deposit (see *Data availability* below).
4. **(Optional)** Re-run the 2-D transport simulations to regenerate the
   aerosol-input profiles, or skip directly to the SARF step using the
   pre-computed profiles included in the Zenodo deposit. <!-- TODO: command(s)
   for the transport step. -->
5. **Run the SARF/ERF sweep** using the scripts under `stardust-climate/`.
   <!-- TODO: list the entry-point script(s) and the configuration files that
   parametrise the published sweeps. -->
6. **Generate the figures** using the post-processing scripts under
   `stardust-climate/transport_paper/scripts/`. <!-- TODO: list the
   plot-generation entry points and the figures they produce. -->

The complete reproduction workflow — environment specification, exact
command-line invocations, expected runtimes and disk footprint — will be
documented here before final publication.

## Data availability

Two public data sources will support reproduction:

1. **Reference-atmospheric-state input database** — a public dataset of
   reference atmospheric state files (zonal-mean climatologies and related
   fields derived from CDS / ERA5 reanalyses) which the analysis code reads
   directly. This database is not yet set up; its location and access details
   will be announced here before final publication.
   <!-- TODO: fill in the public URL / DOI of the CDS-ERA5-derived input
   database once the source is set up. -->
2. **Zenodo deposit** — a versioned archive holding the refractive-index data,
   the Mie-derived optical tables, the pre-computed 2-D aerosol-transport
   profiles, and the per-figure plot data, sufficient for re-creating every
   number and figure in the paper without re-running the upstream simulations.
   <!-- TODO: Zenodo DOI to be added after the first deposit is created. -->

Both will be cited in the paper's *Data availability* statement.

## Contact

- **Scientific / paper correspondence** — corresponding author **Yoav
  Lederer**, [y.lederer@stardust-initiative.com](mailto:y.lederer@stardust-initiative.com)
  (Stardust-Labs, Ness Ziona, Israel).
- **Code / repository correspondence** — repository maintainer **Dorri
  Halbertal**, [d.halbertal@stardust-initiative.com](mailto:d.halbertal@stardust-initiative.com)
  (Stardust-Labs).

For routine bug reports, documentation problems, or reproducibility
questions, please open a GitHub issue on this repository rather than emailing
directly. For matters that should be reported privately, see
[`SECURITY.md`](SECURITY.md).

## License

This umbrella repository's content (documentation, citation metadata,
configuration stubs) is released under the MIT License — see
[`LICENSE`](LICENSE). The individual component repositories will carry their
own licenses; please refer to each component for details.

## Citation

If you use any of the materials linked from this repository, please cite the
paper (see [`CITATION.cff`](CITATION.cff) or use the "Cite this repository"
button on GitHub) and, where appropriate, the Zenodo deposit. Citation details
will be finalized on acceptance.
