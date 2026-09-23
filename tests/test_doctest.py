import doctest
import importlib
import pkgutil

import pytest

import peptacular as pt


def get_all_modules(package, prefix=""):
    """Recursively get all modules and submodules from a package."""
    modules = []

    # Add the package itself if it has a __file__ (is a module)
    if hasattr(package, "__file__"):
        modules.append(package)

    # Get the package path
    if hasattr(package, "__path__"):
        # Iterate through all submodules
        for _importer, modname, _ispkg in pkgutil.walk_packages(
            path=package.__path__,
            prefix=f"{package.__name__}.",
            onerror=lambda x: None,  # Skip modules that fail to import
        ):
            try:
                module = importlib.import_module(modname)
                modules.append(module)
            except Exception as e:
                # Skip modules that can't be imported
                print(f"Skipping {modname}: {e}")
                continue

    return modules


# Build module list
modules = [pt.isotope, pt.regex_utils]

# Add all sequence submodules
modules.extend(get_all_modules(pt.sequence))
modules.extend(get_all_modules(pt.annotation))
modules.extend(get_all_modules(pt.digestion))


@pytest.mark.parametrize("module", modules)
def test_doctests(module):
    """Test doctests for each module."""
    result = doctest.testmod(module, verbose=False)
    assert result.failed == 0, f"Doctests failed in {module.__name__}: {result.failed} failures out of {result.attempted} tests"
