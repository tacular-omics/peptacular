"""Optional console entry point with no SDK import on base-package startup."""

import argparse
import json
import logging
import sys
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description="Peptacular local stdio MCP server")
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--read-root", type=Path, action="append", default=[])
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--ttl-seconds", type=int, default=86400)
    parser.add_argument("--storage-mib", type=int, default=256)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("command", nargs="?", choices=["cache"])
    parser.add_argument("action", nargs="?", choices=["clean"])
    args = parser.parse_args(argv)
    if bool(args.command) != bool(args.action):
        parser.error("Use cache clean together")
    try:
        from .config import Config
        from .operations import versions
        from .server import create_server
        from .storage import Store
    except ImportError:
        print('MCP support is unavailable. Install with: pip install "peptacular[mcp]"', file=sys.stderr)
        return 2
    try:
        config = Config(
            workspace=args.workspace,
            cache=args.cache,
            read_roots=tuple(args.read_root),
            output_root=args.output_root,
            workers=args.workers,
            ttl_seconds=args.ttl_seconds,
            storage_bytes=args.storage_mib * 1024 * 1024,
        )
        if args.check:
            print(
                json.dumps(
                    {
                        "ok": True,
                        "workspace": str(config.workspace),
                        "cache": str(config.cache),
                        "output_root": str(config.output_root),
                        "read_roots": [str(p) for p in config.read_roots],
                        "versions": versions(),
                        "limits": config.limits(),
                    }
                )
            )
            return 0
        if args.command:
            print(json.dumps({"removed_expired_objects": Store(config).clean()}))
            return 0
        logging.basicConfig(stream=sys.stderr, level=logging.WARNING)
        logging.captureWarnings(True)
        create_server(config).run(transport="stdio")
        return 0
    except (ValueError, OSError) as exc:
        print(f"Peptacular MCP: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
