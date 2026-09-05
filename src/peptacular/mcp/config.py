"""Explicit local roots and conservative service budgets."""

import hashlib
import os
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Config:
    workspace: Path
    cache: Path | None = None
    read_roots: tuple[Path, ...] = ()
    output_root: Path | None = None
    workers: int = 2
    ttl_seconds: int = 86400
    storage_bytes: int = 256 * 1024 * 1024
    input_bytes: int = 16 * 1024 * 1024
    max_records: int = 50000
    max_residues: int = 2000000
    max_sequence_length: int = 100000
    page_bytes: int = 256000
    preview_bytes: int = 16000
    inline_seconds: int = 5

    def __post_init__(self):
        workspace = self.workspace.expanduser().resolve(strict=True)
        if not workspace.is_dir():
            raise ValueError("Workspace must be a directory")
        object.__setattr__(self, "workspace", workspace)
        roots = tuple(root.expanduser().resolve(strict=True) for root in self.read_roots)
        if any(not root.is_dir() for root in roots):
            raise ValueError("Read roots must be directories")
        object.__setattr__(self, "read_roots", (workspace, *roots))
        output = (self.output_root or workspace).expanduser().resolve(strict=True)
        if not output.is_dir():
            raise ValueError("Output root must be a directory")
        object.__setattr__(self, "output_root", output)
        if self.cache is None:
            base = (
                Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData/Local"))
                if os.name == "nt"
                else Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
            )
            digest = hashlib.sha256(str(workspace).encode()).hexdigest()[:24]
            object.__setattr__(self, "cache", base / "peptacular" / "mcp" / digest)
        else:
            object.__setattr__(self, "cache", self.cache.expanduser().resolve())
        for name in (
            "workers",
            "ttl_seconds",
            "storage_bytes",
            "input_bytes",
            "max_records",
            "max_residues",
            "max_sequence_length",
            "page_bytes",
            "preview_bytes",
            "inline_seconds",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        object.__setattr__(self, "workers", min(self.workers, 2, os.cpu_count() or 1))

    def limits(self):
        return {key: value for key, value in asdict(self).items() if isinstance(value, int)}

    def input_path(self, supplied):
        path = Path(supplied).expanduser()
        path = (path if path.is_absolute() else self.workspace / path).resolve(strict=True)
        if not any(path.is_relative_to(root) for root in self.read_roots) or not path.is_file():
            raise ValueError("Input must be a regular file inside a configured read root")
        return path

    def destination(self, supplied):
        relative = Path(supplied)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Destination must be relative to the configured output root")
        assert self.output_root is not None
        path = (self.output_root / relative).resolve()
        if not path.is_relative_to(self.output_root) or not path.parent.is_dir():
            raise ValueError("Destination must have an existing parent inside the output root")
        return path
