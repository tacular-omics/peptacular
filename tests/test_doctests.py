"""Execute every docstring example in the package."""

import doctest
import importlib
import importlib.util
import pkgutil

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
MIN_EXAMPLES = 217
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
    """The gate's verdict must match what actually fails to import.

    `MIN_IMPORTED` is computed by subtraction, so it cannot see a wrong `NEEDS_EXTRA` root:
    a verdict stuck at absent skips a module that imports fine, and a verdict stuck at
    present hides a module that does not. This imports every gated module and compares.
    Nothing here inspects the modules' own import statements, because a gated module may
    reach its root through a relative import and still be keyed correctly:
    `peptacular.mcp.operations` is keyed to pydantic and reaches it through `.contracts`.

    :return: None
    :rtype: None
    """
    absent_roots = sorted({root for root in NEEDS_EXTRA.values() if importlib.util.find_spec(root) is None})
    predicted_absent = bool(absent_roots)

    failures: list[tuple[str, str, str]] = []
    for name in sorted(NEEDS_EXTRA):
        try:
            importlib.import_module(name)
        except ImportError as exc:
            failures.append((name, (getattr(exc, "name", None) or "").split(".")[0], f"{type(exc).__name__}: {exc}"))

    named = [f"{module} (root {root or 'unnamed'})" for module, root, _message in failures]
    assert predicted_absent == bool(failures), (
        f"the gate says the extra is {'absent' if predicted_absent else 'present'} "
        f"(absent roots: {absent_roots or 'none'}), but {len(failures)} of {len(NEEDS_EXTRA)} gated modules failed to import: {named or 'none'}"
    )

    unkeyed = [(module, root, message) for module, root, message in failures if root != NEEDS_EXTRA[module]]
    assert not unkeyed, "a gated module failed on a root its key does not name, so it is broken rather than ungated: " + ", ".join(
        f"{module} is keyed to {NEEDS_EXTRA[module]} but failed on {root or 'an unnamed root'} ({message})" for module, root, message in unkeyed
    )

    if predicted_absent:
        assert failures, f"the gate says {absent_roots} absent, which skips every module keyed to them, yet all {len(NEEDS_EXTRA)} imported"
