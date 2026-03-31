# Contributing to Peptacular

Thank you for your interest in contributing to peptacular, a ProForma 2.1 compliant peptide sequence library.
Whether you are fixing a bug, adding a feature, or improving documentation, your help is welcome.
This project is currently under JOSS review, so clear and well-tested contributions are especially valuable.

All participants are expected to follow our [Code of Conduct](CODE_OF_CONDUCT.md) (Contributor Covenant 2.0).

## Reporting Bugs and Requesting Features

Please use [GitHub Issues](https://github.com/tacular-omics/peptacular/issues) for bug reports and feature requests.

**Bug reports** should include:

- Python version and peptacular version (`python --version`, `pip show peptacular`)
- A minimal reproducing example, ideally with a ProForma sequence string
- Expected vs. actual behavior

**Feature requests** should describe the scientific or practical use case, not just the desired API.

## Development Setup

Prerequisites: **Python 3.12+**, [uv](https://docs.astral.sh/uv/), and [just](https://just.systems/).

```bash
# Fork and clone the repository
git clone https://github.com/<your-username>/peptacular.git
cd peptacular

# Install all dependencies (dev + extras)
just install

# Verify everything works
just check   # runs lint, tests, and type checking
```

There are no pre-commit hooks configured; please run `just check` before pushing.

## Making Changes

1. Create a branch from `main` with a descriptive name: `fix/mass-rounding`, `feature/new-ion-type`, `docs/isotope-example`.
2. Keep each pull request focused on a single change.
3. Write clear commit messages that explain *why*, not just *what*.

## Code Style

Formatting and linting are handled by [ruff](https://docs.astral.sh/ruff/) (line length 160, Python 3.12 target):

```bash
just format   # auto-format code
just lint     # check for lint errors
```

**Docstrings** use Sphinx style:

```python
def mass(sequence: str, monoisotopic: bool = True) -> float:
    """Calculate the neutral mass of a peptide sequence.

    :param sequence: A ProForma-formatted peptide string.
    :type sequence: str
    :param monoisotopic: Use monoisotopic masses if True, average if False.
    :type monoisotopic: bool
    :return: The neutral mass in Daltons.
    :rtype: float
    :raises ValueError: If the sequence cannot be parsed.
    """
```

**Type annotations** are required on all public functions. Run `just ty` to check.

## Testing

- All new code needs tests. Tests live in `tests/`.
- Use `pytest.approx` for floating-point comparisons.
- Use the `tmp_path` fixture for any file I/O; do not mock library internals.
- Aim to maintain at least **83% branch coverage**.

```bash
just test       # run tests
just test-cov   # run tests with coverage report
```

## Documentation

- Add docstrings to all public functions and classes.
- Build the docs locally with `just docs` and verify your changes render correctly.
- When introducing new functionality, consider adding a worked example to the `examples/` directory.

## Submitting a Pull Request

Before opening a PR, confirm:

- [ ] `just check` passes (lint + tests + type checking)
- [ ] New code has tests and docstrings
- [ ] PR description explains the change and its motivation
- [ ] Any related issue is referenced (e.g., "Closes #42")

CI (GitHub Actions) will run lint, type checking, and tests automatically when you push. A maintainer will review your PR and may request changes.

## Questions?

Open an [issue](https://github.com/tacular-omics/peptacular/issues) or reach out to the maintainers. We are happy to help.
