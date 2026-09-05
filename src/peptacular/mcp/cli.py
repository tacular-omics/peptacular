"""Optional console entry point with no SDK import on base-package startup."""

import argparse
import json
import logging
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description="Peptacular stateless stdio MCP server")
    parser.add_argument("--check", action="store_true", help="Report installed capabilities and exit")
    args = parser.parse_args(argv)
    try:
        from .contracts import REQUESTS
        from .operations import LIMITS, versions
        from .server import create_server
    except ImportError:
        print('MCP support is unavailable. Install with: pip install "peptacular[mcp]"', file=sys.stderr)
        return 2
    if args.check:
        print(json.dumps({"ok": True, "tools": list(REQUESTS), "versions": versions(), "limits": LIMITS}))
        return 0
    logging.basicConfig(stream=sys.stderr, level=logging.WARNING)
    logging.captureWarnings(True)
    create_server().run(transport="stdio")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
