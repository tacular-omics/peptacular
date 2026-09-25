"""
pytest configuration file for peptacular tests.

Fast by default. Two switches bring back the full run that CI uses:

- ``--run-slow`` (or ``RUN_SLOW=1``) runs tests marked ``slow``: MCP stdio subprocess
  launches and the AlphaBase interop tests (importing AlphaBase JIT-compiles numba code).
- ``HYPOTHESIS_PROFILE=ci`` runs every property test with 150 examples (CI);
  ``HYPOTHESIS_PROFILE=thorough`` runs 2000. The default profile runs 25.
"""

import os
import sys
from pathlib import Path

import pytest
from hypothesis import HealthCheck, settings

# Add src directory to Python path
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

_HYPOTHESIS_COMMON = {"deadline": None, "suppress_health_check": [HealthCheck.too_slow]}
settings.register_profile("default", max_examples=25, **_HYPOTHESIS_COMMON)
settings.register_profile("ci", max_examples=150, **_HYPOTHESIS_COMMON)
settings.register_profile("thorough", max_examples=2000, **_HYPOTHESIS_COMMON)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "default"))


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption("--run-slow", action="store_true", default=False, help="also run tests marked slow (or set RUN_SLOW=1)")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if config.getoption("--run-slow") or os.environ.get("RUN_SLOW", "") not in ("", "0"):
        return
    skip_slow = pytest.mark.skip(reason="slow: run with --run-slow or RUN_SLOW=1")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
