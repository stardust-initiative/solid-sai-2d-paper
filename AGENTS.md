# AGENTS

> *AI coding agents (Cursor, Claude Code, GitHub Copilot, etc.): start here.
> This file is the single source of truth for repository-level conventions in
> `stardust-initiative/solid-sai-2d-paper`. See also [CLAUDE.md](CLAUDE.md)
> (which redirects here) and [CONTRIBUTING.md](CONTRIBUTING.md) (which targets
> human collaborators).*

## What this repository is

The umbrella repository accompanying the paper *"Solid-particle
stratospheric aerosol injection: a 2-D modeling exploration of the design
space"* (Lederer et al., submitted to *Atmospheric Chemistry and Physics*,
2026).

At present (skeleton state, `v0.x`) it holds only documentation, citation
metadata, and CI configuration. By `v1.0` (paper acceptance) it will
additionally host the **paper-specific scripts**: those that drive the
component-repo simulations and those that post-analyze their output into
the results discussed in the manuscript.

The *generic* simulation machinery (radiation, transport, optical tables,
modified `climlab` schemes) lives in the **component repositories** listed
in [README.md](README.md) — those are reusable libraries the paper-specific
scripts call into.

## Layout

| path | role |
|---|---|
| `README.md` | primary entry point — status, paper citation, component list, reproduction-workflow placeholders, contacts. |
| `LICENSE` | MIT. |
| `CITATION.cff` | machine-readable citation metadata; validated in CI on every PR. |
| `CHANGELOG.md` | release history; *Keep a Changelog* format. |
| `CONTRIBUTING.md` | contribution rules for humans. |
| `SECURITY.md` | private-disclosure channel. |
| `AGENTS.md` | this file. |
| `CLAUDE.md` | redirect to this file. |
| `.github/CODEOWNERS` | review routing to the repository maintainer. |
| `.github/pull_request_template.md` | PR template — flags paper-citation implications of changes. |
| `.github/workflows/cff-validation.yml` | `cffconvert --validate` against `CITATION.cff`. |
| `.gitignore` | scientific-Python defaults + safety net against accidental data commits. |

## Setup / build / test

In skeleton state, no setup is required. The only automated check is
`cffconvert --validate` against `CITATION.cff`, which runs on every push and
PR to `main` that touches the citation file.

Once the paper-specific scripts are committed, this section will document
the Python environment (version pin, dependencies, install command), how
the scripts are invoked, and any data-fetching steps required from the
planned public input database and Zenodo deposit.

## How to contribute

1. **Fork** the repository and create a feature branch off `main`:
   ```
   git checkout -b <type>/<short-desc>
   ```
   where `<type>` ∈ {`feat`, `fix`, `docs`, `chore`, `refactor`, `test`,
   `perf`, `ci`, `build`}.
2. Use **[Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/)**
   for commit messages: `<type>[optional-scope]: <description>`.
   Examples: `docs(readme): correct transport-2d framing`,
   `chore(citation): set date-released to v0.2 tag date`,
   `ci(cff-validation): pin action to v2.0.1`.
3. Open a **pull request** against `main`. Direct pushes to `main` are
   blocked by repository rules.
4. Merge strategy is **rebase-only** (no merge commits, no squash); `main`
   keeps a linear history.
5. The PR must pass the `cff-validation` check (when applicable) before it
   can be merged.

## Scope — what does and does not belong here

- **Belongs here**: documentation, citation metadata, CHANGELOG, and the
  paper-specific scripts — those that orchestrate the component-repo
  simulations, and those that post-analyze their outputs into the results
  discussed in the manuscript.
- **Does not belong here**: the *generic* simulation machinery (lives in the
  component repositories), the *raw* or *pre-computed* simulation outputs
  (will live in the planned Zenodo deposit), and the *reference-atmospheric-
  state input data* (will live in the planned CDS / ERA5-derived public
  input database). See README's *Data availability* section.
- Once `v1.0` is tagged at paper acceptance, the contents are
  reproducibility-frozen: further changes must be either an erratum (with a
  `CHANGELOG.md` entry under a new patch tag) or a metadata-only update
  (e.g. inserting a DOI).

## Security / sensitive content

This repository must never contain credentials, private paths, internal
URLs, or unpublished data. Push-protection (secret-scanning) is enabled on
the repository as a safety net. To report a sensitive issue privately, see
[SECURITY.md](SECURITY.md).

## Org-wide conventions referenced here

- **Conventional Commits 1.0.0** — <https://www.conventionalcommits.org/en/v1.0.0/>
- **Semantic Versioning** — <https://semver.org/>
  (interpreted for paper-companion repos as: major = paper acceptance,
  minor = open-review updates, patch = corrections; see `CHANGELOG.md`).
- **AGENTS.md convention** — <https://agents.md>
