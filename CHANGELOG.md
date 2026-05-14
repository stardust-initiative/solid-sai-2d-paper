# Changelog

All notable changes to this repository are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project loosely
adheres to [Semantic Versioning](https://semver.org/), interpreted as:

- **major (`v1.0`)** — paper accepted; repository state matches the final
  version of record. Subsequent major bumps are not anticipated.
- **minor (`v0.x`)** — substantive updates during the open-review period
  (e.g. component repositories opened, DOIs added, reproduction workflow
  fleshed out). A post-publication minor bump (`v1.1`, `v1.2`, …) would
  only accompany a published corrigendum or addendum.
- **patch (`v0.x.y`, `v1.0.y`, …)** — corrections, typo fixes, and other
  changes that do not affect the science (including post-publication errata).

## [Unreleased]

## [v0.1] — 2026-05-13

Initial submission-state skeleton accompanying the manuscript *"Solid-particle
stratospheric aerosol injection: a 2-D modeling exploration of the design
space"* (Lederer, Wygoda, Halbertal, and Rose, submitted to *Atmospheric
Chemistry and Physics*, 2026).

### Added — documentation and metadata

- `README.md` — planned umbrella structure, component repositories
  (`stardust-climate`, `climlab-stardust-extension`, `climlab-rrtmg_stardust`,
  `climlab-sbm-convection_stardust`, `optical-tables-generator`), the two
  planned public data sources (CDS / ERA5-derived reference-state input
  database, and a Zenodo deposit for refractive-index data, optical tables,
  transport profiles, and per-figure plot data), and the split-routing
  contact section (paper → Yoav Lederer; code / repository → Dorri Halbertal).
- `CITATION.cff` — author metadata and a `preferred-citation` block for the
  paper (DOI to be added on acceptance).
- `LICENSE` — MIT.
- `CONTRIBUTING.md` — slim contribution rules for humans.
- `SECURITY.md` — private-disclosure channel.
- `AGENTS.md` — single source of truth for repository-level conventions
  (scope, layout, contribution rules, what does and does not belong in this
  repository).
- `CLAUDE.md` — one-paragraph redirect to `AGENTS.md`.
- `CHANGELOG.md` — this file.

### Added — CI and repository hygiene

- `.github/CODEOWNERS` — review routing to the repository maintainer
  (`@dorrih-stardust`).
- `.github/pull_request_template.md` — PR template flagging the
  paper-citation implications of changes.
- `.github/workflows/cff-validation.yml` — runs `cffconvert --validate` on
  every push and pull request against `main`; serves as a required status
  check for the branch's protection ruleset.
- `.gitignore` — scientific-Python defaults plus a safety net against
  accidental data commits.

### Notes

- Component repositories are private at this stage and will be made public,
  archived on Zenodo, and assigned permanent DOIs prior to the publication of
  the final-revised version.
