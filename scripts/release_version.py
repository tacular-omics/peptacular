"""Synchronize and verify release versions without importing the package.

Shared across the tacular-omics packages. The canonical copy lives in the
tacular-omics workspace at templates/scripts/release_version.py; edit it there
and re-sync, do not edit a package's copy by hand.

The version lives only in ``__version__`` in the file named by
``[tool.hatch.version] path``. ``sync`` copies it to CITATION.cff (and
.zenodo.json when present). ``sync --set X.Y.Z`` also turns the changelog's
``## [Unreleased]`` section into ``## [X.Y.Z] (today)``.
"""

import argparse
import ast
import json
import re
import tarfile
import tomllib
import zipfile
from datetime import date
from email.parser import Parser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
NAME = CONFIG["project"]["name"]
DIST = NAME.replace("-", "_")
SOURCE = Path(CONFIG["tool"]["hatch"]["version"]["path"])
VERSION = r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)(?:(?:a|b|rc)[0-9]+)?(?:\.post[0-9]+)?(?:\.dev[0-9]+)?"


def validate_version(value: str) -> str:
    if not re.fullmatch(VERSION, value):
        raise ValueError(f"Invalid version {value!r}. Use X.Y.Z, optionally with a PEP 440 a, b, rc, post, or dev suffix.")
    return value


def source_version(text: str) -> str:
    values = [
        ast.literal_eval(node.value)
        for node in ast.parse(text).body
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets)
    ]
    if len(values) != 1 or not isinstance(values[0], str):
        raise ValueError("Expected exactly one literal __version__ assignment")
    return validate_version(values[0])


def replace_field(text: str, pattern: str, value: str) -> str:
    updated, count = re.subn(pattern, lambda match: match[1] + json.dumps(value), text, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"Expected exactly one version field matching {pattern!r}")
    return updated


def metadata(root: Path) -> tuple[str, dict[Path, str]]:
    source = (root / SOURCE).read_text(encoding="utf-8")
    version = source_version(source)
    config = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    if config["tool"]["hatch"]["version"]["path"] != SOURCE.as_posix() or "version" not in config["project"].get("dynamic", []):
        raise ValueError(f"Package builds must read the dynamic version from {SOURCE}")
    if "version" in config["project"]:
        raise ValueError("Remove the duplicated static project version from pyproject.toml")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    updates = {root / "CITATION.cff": replace_field(citation, r"^(version: *)[^\n]+$", version)}
    if (root / ".zenodo.json").exists():
        zenodo = json.loads((root / ".zenodo.json").read_text(encoding="utf-8"))
        zenodo["version"] = version
        updates[root / ".zenodo.json"] = json.dumps(zenodo, indent=2, ensure_ascii=False) + "\n"
    return version, updates


def sync(root: Path, new_version: str | None = None) -> str:
    version, updates = metadata(root)
    if new_version is not None:
        version = validate_version(new_version)
        updates[root / SOURCE] = replace_field((root / SOURCE).read_text(encoding="utf-8"), r"^(__version__ = )[^\n]+$", version)
        updates[root / "CITATION.cff"] = replace_field(updates[root / "CITATION.cff"], r"^(version: *)[^\n]+$", version)
        if root / ".zenodo.json" in updates:
            zenodo = json.loads(updates[root / ".zenodo.json"])
            zenodo["version"] = version
            updates[root / ".zenodo.json"] = json.dumps(zenodo, indent=2, ensure_ascii=False) + "\n"
        changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
        if f"## [{version}]" not in changelog:
            heading = f"## [Unreleased]\n\n## [{version}] ({date.today().isoformat()})"
            changelog, count = re.subn(r"^## \[Unreleased\][^\n]*$", heading, changelog, count=1, flags=re.MULTILINE)
            if count != 1:
                raise ValueError("CHANGELOG.md needs a '## [Unreleased]' section to release from")
            updates[root / "CHANGELOG.md"] = changelog
    for path, content in updates.items():
        if path.read_text(encoding="utf-8") != content:
            path.write_text(content, encoding="utf-8")
    return version


def check(root: Path, tag: str | None = None) -> str:
    version, _ = metadata(root)
    citation = re.findall(r"^version: *([^\n]+)$", (root / "CITATION.cff").read_text(encoding="utf-8"), re.MULTILINE)
    actual = {"CITATION.cff": citation[0].strip().strip("\"'")}
    if (root / ".zenodo.json").exists():
        actual[".zenodo.json"] = json.loads((root / ".zenodo.json").read_text(encoding="utf-8"))["version"]
    stale = [path for path, value in actual.items() if value != version]
    if stale:
        raise ValueError(f"Stale version metadata: {', '.join(stale)}. Run just sync-version.")
    if tag is not None:
        if tag != f"v{version}":
            raise ValueError(f"Release tag {tag!r} does not match v{version}")
        headings = re.findall(r"^## \[([^\]]+)\]([^\n]*)", (root / "CHANGELOG.md").read_text(encoding="utf-8"), re.MULTILINE)
        headings = [h for h in headings if h[0] != "Unreleased"]
        if not headings or headings[0][0] != version:
            raise ValueError(f"The first changelog release must be [{version}]")
        dated = re.fullmatch(r" \((\d{4}-\d{2}-\d{2})\)", headings[0][1])
        if dated is None:
            raise ValueError(f"Date the [{version}] changelog entry before publishing, using (YYYY-MM-DD)")
        date.fromisoformat(dated[1])
    return version


def check_metadata(text: str, version: str, artifact: Path) -> None:
    fields = Parser().parsestr(text)
    if fields.get_all("Name") not in ([NAME], [DIST]) or fields.get_all("Version") != [version]:
        raise ValueError(f"{artifact.name}: package metadata must identify {NAME} {version}")


def check_artifacts(directory: Path, version: str) -> None:
    wheels = list(directory.glob("*.whl"))
    sdists = list(directory.glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise ValueError("Expected exactly one wheel and one source distribution. Build into an empty output directory.")
    if not wheels[0].name.startswith(f"{DIST}-{version}-") or sdists[0].name != f"{DIST}-{version}.tar.gz":
        raise ValueError("Distribution filenames differ from the release version")
    with zipfile.ZipFile(wheels[0]) as archive:
        names = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(names) != 1:
            raise ValueError("Wheel must contain exactly one METADATA file")
        check_metadata(archive.read(names[0]).decode(), version, wheels[0])
        if source_version(archive.read(SOURCE.relative_to("src").as_posix()).decode()) != version:
            raise ValueError("Wheel runtime version differs from the release version")
    with tarfile.open(sdists[0]) as archive:
        names = [name for name in archive.getnames() if name.count("/") == 1 and name.endswith("/PKG-INFO")]
        if len(names) != 1:
            raise ValueError("Source distribution must contain exactly one root PKG-INFO file")
        info = archive.extractfile(names[0])
        source = archive.extractfile(names[0].split("/")[0] + "/" + SOURCE.as_posix())
        if info is None or source is None:
            raise ValueError("Source distribution is missing release metadata or package source")
        check_metadata(info.read().decode(), version, sdists[0])
        if source_version(source.read().decode()) != version:
            raise ValueError("Source distribution runtime version differs from the release version")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    update = commands.add_parser("sync", help="Copy the package version to citation and Zenodo metadata")
    update.add_argument("--set", dest="new_version", help="Set the package version and synchronize metadata in one command")
    verify = commands.add_parser("check", help="Fail if release versions disagree")
    verify.add_argument("--tag", help="Also require an exact vX.Y.Z tag and dated changelog entry")
    verify.add_argument("--dist", type=Path, help="Also verify built wheel and source distribution versions")
    args = parser.parse_args()
    try:
        if args.command == "sync":
            version = sync(ROOT, args.new_version)
            print(f"Synchronized version {version}")
        else:
            version = check(ROOT, args.tag)
            if args.dist is not None:
                check_artifacts(args.dist, version)
            print(f"Version checks passed for {version}")
    except (ValueError, OSError, KeyError, SyntaxError, tarfile.TarError, zipfile.BadZipFile) as error:
        parser.exit(1, f"Version check failed: {error}\n")


if __name__ == "__main__":
    main()
