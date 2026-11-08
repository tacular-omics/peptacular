## What this changes

<!-- One or two sentences. If it fixes an open issue, write "Fixes #123" so GitHub links them. -->

## Why

<!-- The behaviour that was wrong, or the gap this fills. If the change follows a specification,
     a controlled vocabulary or a published method, name it so a reviewer can check it. -->

## How it was tested

<!-- Which tests you ran and what they reported. New behaviour needs a new test. -->

- [ ] `just check` passes locally, or CI on this pull request is green
- [ ] New or changed behaviour has a test that fails without this change
- [ ] Public API changes are reflected in the docstrings and the docs
- [ ] `CHANGELOG.md` has an entry under `## [Unreleased]`

## Notes for the maintainer

<!-- Anything you were unsure about, or a decision you would like checked. -->

---

Two things that are the maintainer's job, not yours, so please leave them alone:

- **Do not change the version** in `pyproject.toml`, `__init__.py` or `CITATION.cff`, and do
  not add a new version heading to `CHANGELOG.md`. Put your entry under `## [Unreleased]` and
  the release will pick it up.
- **Do not edit `uv.lock`** unless your change adds or removes a dependency. If it does, say so
  here and a maintainer will refresh the lock.

If this is your first pull request here, welcome, and thank you. Setting up the development
environment is `uv sync --locked` from a clone, which installs everything the test suite needs,
including the optional readers. `CONTRIBUTING.md` has the rest.
