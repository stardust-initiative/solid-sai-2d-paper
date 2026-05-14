# Contributing

This repository is the umbrella for a published (or, at the time of writing,
under-review) paper. Its primary purpose is to provide a stable public entry
point to the code and data used in that paper.

## What can change

Once a `v1.0` release is tagged to match the paper's version of record, the
**contents of the umbrella repository are frozen** in the sense that any
change which would alter the reproducibility of the published results will
**not** be merged into `main`; such changes belong in follow-up work and in
the relevant component repositories' own development branches.

What *will* be updated post-`v1.0` here is purely descriptive:

- Errata and clarifications in `README.md`.
- DOIs (paper, Zenodo) once assigned.
- Links to component repositories as they are opened.

## Issues and discussion are welcome

Please open a GitHub issue if you:

- find a problem reproducing a result described in the paper,
- spot a bug or inconsistency in the documentation,
- have a question about the methodology that the paper does not answer, or
- want to flag an erratum that we should record here.

For substantive discussion that is not a bug report, the GitHub *Discussions*
tab on this repository is the preferred venue.

For correspondence outside GitHub, please contact the corresponding author
(see [`README.md`](README.md)).

## Component-repository contributions

If your contribution is to one of the *component* repositories
(`stardust-climate`, `climlab-stardust-extension`, etc.), please open the
issue or pull request on that component's own repository, not here.
