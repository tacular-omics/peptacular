"""Lazy optional-dependency loading."""

from importlib import import_module
from types import ModuleType

from ._errors import MissingOptionalDependencyError


def require_dependency(module_name: str, extra_name: str) -> ModuleType:
    """Import an optional dependency or raise an actionable error.

    :param module_name: Importable top-level module name.
    :type module_name: str
    :param extra_name: Peptacular extra that installs the dependency.
    :type extra_name: str
    :return: Imported module.
    :rtype: ModuleType
    :raises MissingOptionalDependencyError: If the requested module is absent.
    """
    try:
        return import_module(module_name)
    except ModuleNotFoundError as exc:
        if exc.name not in {module_name, module_name.partition(".")[0]}:
            raise
        raise MissingOptionalDependencyError(
            f"The '{module_name.partition('.')[0]}' package is required for this conversion. Install it with 'pip install peptacular[{extra_name}]'."
        ) from exc
