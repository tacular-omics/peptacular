import pytest

pytest.importorskip("mcp")
pytest.importorskip("pytest_asyncio")

from peptacular.mcp.config import Config
from peptacular.mcp.storage import Store


@pytest.fixture
def config(tmp_path):
    return Config(tmp_path, cache=tmp_path / "cache")


@pytest.fixture
def store(config):
    return Store(config)
