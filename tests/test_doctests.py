"""Execute every docstring example in the package."""

import doctest
import importlib
import pkgutil

import peptacular

# Modules that need an optional extra and may therefore be absent. A module is skipped only
# if it is named here AND its import raises ImportError: anything else propagates, and a
# module that is not named here must import. Measured 2026-11-08 in an environment with no
# extras installed: all 74 modules the walk yields import except these four, and
# `peptacular.mcp` itself imports (its __init__ is empty), so the walk still reaches them.
OPTIONAL_MODULES = frozenset(
    {
        # The `mcp` extra: mcp and pydantic.
        "peptacular.mcp.contracts",
        "peptacular.mcp.operations",
        "peptacular.mcp.outputs",
        "peptacular.mcp.server",
    }
)

# Measured 2026-11-08 at 0a4d366 with `uv sync --all-extras`. The module count is the names
# the walk yields, which is the same with and without extras. Raise these when you add
# modules or examples, and never lower them without saying why in the commit message.
MIN_MODULES = 74
MIN_EXAMPLES = 217


def _walk_module_names() -> list[str]:
    """Every importable module name in the package, from a walk rather than a hand-list.

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
    """Run every docstring example in every module and hold both measured floors.

    :return: None
    :rtype: None
    """
    names = _walk_module_names()
    assert len(names) >= MIN_MODULES, f"module walk found only {len(names)}, expected >= {MIN_MODULES}"

    attempted = 0
    failures: list[str] = []
    skipped: list[str] = []
    for name in names:
        try:
            module = importlib.import_module(name)
        except ImportError:
            if name not in OPTIONAL_MODULES:
                raise
            skipped.append(name)
            continue
        result = doctest.testmod(module, verbose=False)
        attempted += result.attempted
        if result.failed:
            failures.append(f"{name} ({result.failed} of {result.attempted})")

    assert not failures, f"docstring examples failed in {len(failures)} modules: {', '.join(failures)}"
    assert attempted >= MIN_EXAMPLES, f"only {attempted} examples ran, expected >= {MIN_EXAMPLES} (skipped: {skipped or 'none'})"
