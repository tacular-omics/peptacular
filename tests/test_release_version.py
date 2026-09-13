import io
import json
import runpy
import tarfile
import zipfile
from pathlib import Path

import pytest

release = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts/release_version.py"))


@pytest.fixture
def project(tmp_path):
    source = tmp_path / "src/peptacular/__init__.py"
    source.parent.mkdir(parents=True)
    source.write_text('raise RuntimeError("Must not import the package")\n__version__ = "3.3.0"\n')
    (tmp_path / "pyproject.toml").write_text(
        '[project]\ndynamic = ["version"]\n[tool.hatch.version]\npath = "src/peptacular/__init__.py"\n'
    )
    (tmp_path / "CITATION.cff").write_text('title: Keep this title\nversion: "3.3.0"\nauthors:\n  - name: Test Author\n')
    (tmp_path / ".zenodo.json").write_text(json.dumps({"version": "3.3.0", "creators": [{"name": "Test Author"}]}))
    (tmp_path / "CHANGELOG.md").write_text("# Changelog\n\n## [4.0.0] (Unreleased)\n\nKeep these release notes.\n")
    return tmp_path


def test_one_command_sets_version_and_preserves_other_metadata(project):
    release["sync"](project, "4.0.0")
    assert release["check"](project) == "4.0.0"
    assert 'raise RuntimeError("Must not import the package")' in (project / "src/peptacular/__init__.py").read_text()
    assert "name: Test Author" in (project / "CITATION.cff").read_text()
    assert json.loads((project / ".zenodo.json").read_text())["creators"] == [{"name": "Test Author"}]
    assert "Keep these release notes." in (project / "CHANGELOG.md").read_text()


def test_manual_version_edit_is_detected_and_synchronized(project):
    (project / "src/peptacular/__init__.py").write_text('__version__ = "4.0.0rc1"\n')
    with pytest.raises(ValueError, match="Stale version metadata"):
        release["check"](project)
    release["sync"](project)
    assert release["check"](project) == "4.0.0rc1"


@pytest.mark.parametrize("version", ["v4.0.0", "04.0.0", "4.0", "4.0.0\n", "4.0.0+bad space"])
def test_invalid_version_does_not_modify_files(project, version):
    before = {p: p.read_bytes() for p in project.rglob("*") if p.is_file()}
    with pytest.raises(ValueError, match="Invalid version"):
        release["sync"](project, version)
    assert before == {p: p.read_bytes() for p in before}


def test_release_requires_matching_tag_and_dated_current_changelog(project):
    release["sync"](project, "4.0.0")
    with pytest.raises(ValueError, match="does not match"):
        release["check"](project, "v3.3.0")
    with pytest.raises(ValueError, match="Date the"):
        release["check"](project, "v4.0.0")
    changelog = project / "CHANGELOG.md"
    changelog.write_text("## [4.0.0] (2026-09-13)\n")
    assert release["check"](project, "v4.0.0") == "4.0.0"
    changelog.write_text("## [3.3.0] (2026-09-13)\n")
    with pytest.raises(ValueError, match="first changelog release"):
        release["check"](project, "v4.0.0")


def make_artifacts(directory, wheel_version="4.0.0", source_version="4.0.0", runtime="4.0.0"):
    wheel = directory / "peptacular-4.0.0-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("peptacular-4.0.0.dist-info/METADATA", f"Name: peptacular\nVersion: {wheel_version}\n")
        archive.writestr("peptacular/__init__.py", f'__version__ = "{runtime}"\n')
    with tarfile.open(directory / "peptacular-4.0.0.tar.gz", "w:gz") as archive:
        for name, text in {
            "PKG-INFO": f"Name: peptacular\nVersion: {source_version}\n",
            "src/peptacular/__init__.py": f'__version__ = "{runtime}"\n',
        }.items():
            info = tarfile.TarInfo("peptacular-4.0.0/" + name)
            info.size = len(text.encode())
            archive.addfile(info, io.BytesIO(text.encode()))
    return wheel


def test_built_artifacts_agree_and_stale_wheel_is_rejected(tmp_path):
    wheel = make_artifacts(tmp_path)
    release["check_artifacts"](tmp_path, "4.0.0")
    (tmp_path / "peptacular-3.3.0-py3-none-any.whl").write_bytes(wheel.read_bytes())
    with pytest.raises(ValueError, match="exactly one wheel"):
        release["check_artifacts"](tmp_path, "4.0.0")


@pytest.mark.parametrize("field", ["wheel_version", "source_version", "runtime"])
def test_artifact_version_mismatch_blocks_release(tmp_path, field):
    make_artifacts(tmp_path, **{field: "3.3.0"})
    with pytest.raises(ValueError, match="metadata|runtime version"):
        release["check_artifacts"](tmp_path, "4.0.0")


def test_artifact_filename_must_match_release(tmp_path):
    wheel = make_artifacts(tmp_path)
    wheel.rename(tmp_path / "peptacular-3.3.0-py3-none-any.whl")
    with pytest.raises(ValueError, match="filenames differ"):
        release["check_artifacts"](tmp_path, "4.0.0")
