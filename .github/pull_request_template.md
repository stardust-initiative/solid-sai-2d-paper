<!--
Thanks for proposing a change to this repository.

This repository is the umbrella for a published (or under-review) paper, so
its content is more conservative than a normal development project: changes
that affect what the paper cites should be flagged explicitly so we can
record them in CHANGELOG.md and decide whether they need a release tag or
an erratum.
-->

## Summary

<!-- One or two sentences describing what this PR changes and why. -->

## Type of change

- [ ] Documentation / README / metadata
- [ ] Citation metadata (`CITATION.cff`)
- [ ] Repository configuration (`.gitignore`, `.github/`, branch / tag
      protection, hooks)
- [ ] Changelog entry only
- [ ] Other (please describe)

## Checklist

- [ ] The change is consistent with the manuscript's *Code availability* and
      *Data availability* statements (or those statements have been updated
      accordingly in a separate PR / will be updated at the next revision).
- [ ] If touching `CITATION.cff`, the file still passes CFF validation
      (the workflow at `.github/workflows/cff-validation.yml` will check
      this automatically on PRs).
- [ ] If touching `README.md`, internal anchors and external links still
      resolve.
- [ ] An entry has been added under `[Unreleased]` in `CHANGELOG.md` if the
      change is user-visible.
- [ ] If this PR alters anything cited from the paper's version of record,
      the README's *Status* section / an erratum has been added to make
      that visible.
