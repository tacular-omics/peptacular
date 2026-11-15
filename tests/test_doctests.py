"""Execute every docstring example in the package."""

import doctest
import importlib
import importlib.util
import pkgutil

import pytest

import peptacular

# Modules that cannot be imported without an optional extra, each mapped to the
# third-party package it needs. The skip is keyed to that package being absent, not to the
# exception, so an ImportError raised for any other reason still fails the test. Measured
# 2026-11-08 in an environment with no extras installed: these four are the only imports
# that fail, and `peptacular.mcp` and `peptacular.mcp.cli` import without the extra, so
# they are not listed here and must import. `server` needs mcp itself, the other three
# reach pydantic through `contracts`.
NEEDS_EXTRA = {
    "peptacular.mcp.contracts": "pydantic",
    "peptacular.mcp.operations": "pydantic",
    "peptacular.mcp.outputs": "pydantic",
    "peptacular.mcp.server": "mcp",
}

# Measured 2026-11-08 at 774f136 with `uv sync --all-extras`. The module count is the names
# the walk yields, which is the same with and without extras, because `peptacular.mcp` has
# an empty __init__ and the walk enters it either way. Raise these when you add modules or
# examples, and never lower them without saying why in the commit message.
MIN_MODULES = 74
# The 17 property examples add an import and a calculation each.
MIN_EXAMPLES = 251
MIN_IMPORTED = MIN_MODULES - len(NEEDS_EXTRA)


def _walk_module_names() -> list[str]:
    """Every module name in the package, from a walk rather than a hand-list.

    :return: Module names, the package itself first.
    :rtype: list[str]
    :raises AssertionError: If walking into a subpackage fails, which must never be silent.
    """

    def _fail(name: str) -> None:
        raise AssertionError(f"walking into {name} failed, so its submodules were not seen")

    names = [peptacular.__name__]
    names.extend(name for _finder, name, _ispkg in pkgutil.walk_packages(peptacular.__path__, prefix=f"{peptacular.__name__}.", onerror=_fail))
    return names


def test_docstring_examples() -> None:
    """Run every docstring example in every module and hold the three measured floors.

    :return: None
    :rtype: None
    """
    names = _walk_module_names()
    assert len(names) >= MIN_MODULES, f"module walk found only {len(names)}, expected >= {MIN_MODULES}"

    stale = sorted(set(NEEDS_EXTRA) - set(names))
    assert not stale, f"NEEDS_EXTRA names modules the walk does not yield, so those entries filter nothing: {stale}"

    attempted = 0
    imported = 0
    failures: list[str] = []
    skipped: list[str] = []
    for name in names:
        required = NEEDS_EXTRA.get(name)
        if required is not None and importlib.util.find_spec(required) is None:
            skipped.append(f"{name} (no {required})")
            continue
        module = importlib.import_module(name)
        imported += 1
        result = doctest.testmod(module, verbose=False)
        attempted += result.attempted
        if result.failed:
            failures.append(f"{name} ({result.failed} of {result.attempted})")

    assert not failures, f"docstring examples failed in {len(failures)} modules: {', '.join(failures)}"
    assert imported >= MIN_IMPORTED, f"only {imported} modules imported, expected >= {MIN_IMPORTED} (skipped: {skipped or 'none'})"
    assert attempted >= MIN_EXAMPLES, f"only {attempted} examples ran, expected >= {MIN_EXAMPLES} (skipped: {skipped or 'none'})"


def test_needs_extra_gate_agrees_with_what_imports() -> None:
    """Every `NEEDS_EXTRA` entry must predict, per module, whether that module imports.

    The floors cannot see a wrong root, because `MIN_IMPORTED` is computed by subtraction:
    a verdict stuck at absent skips a module that imports fine, and a verdict stuck at
    present hides one that does not. This asserts the gate per module, which also catches
    an over-listed root, a module skipped in a lean install for a package it never needed.

    The invariant this rests on: each entry's root must be one the module genuinely needs,
    directly or through its own relative imports, which is what makes "root absent implies
    the import fails" true. That is not the same as the module's own module-scope import
    list, and here it is demonstrably not: `peptacular.mcp.operations` is keyed to pydantic
    while its own roots name only tacular, because it reaches pydantic through a relative
    import of `.contracts`. The strict form is still sound for it, since pydantic absent
    makes `.contracts` fail and `operations` with it. The shape that breaks this form is one
    entry shared by several modules whose roots differ, where a half-installed environment
    leaves some of them importable, so each module here keeps its own entry.

    :return: None
    :rtype: None
    """
    assert NEEDS_EXTRA, "NEEDS_EXTRA is empty, so the gate is inert and nothing is checked"
    rootless = sorted(name for name, root in NEEDS_EXTRA.items() if not root)
    assert not rootless, f"these NEEDS_EXTRA entries name no root, so they read as present in every environment: {rootless}"

    for name, root in sorted(NEEDS_EXTRA.items()):
        absent = importlib.util.find_spec(root) is None
        try:
            importlib.import_module(name)
        except pytest.skip.Exception as exc:
            pytest.fail(f"importing {name} raised Skipped, which pytest reports as a skip rather than a failure: {exc}")
        except ImportError as exc:
            failed_root = (getattr(exc, "name", None) or "").split(".")[0]
            assert absent, (
                f"{name} failed to import although its key's root {root} is present, so it is broken rather than ungated: {type(exc).__name__}: {exc}"
            )
            assert failed_root == root, (
                f"{name} is keyed to {root} but failed on {failed_root or 'an unnamed root'}, "
                f"so the key does not name what actually went missing: {type(exc).__name__}: {exc}"
            )
        else:
            assert not absent, (
                f"{name} imported although {root} is absent, so its key names a package the module does not need "
                f"and a lean install skips examples it could have run"
            )
