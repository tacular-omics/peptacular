"""Install the freshly built wheel into a temporary environment and smoke test it."""

import argparse
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--extra", choices=["pyteomics", "psm-utils", "alphabase", "mcp", "mcp,pyteomics", "mcp,psm-utils", "mcp,alphabase"])
args = parser.parse_args()
wheels = sorted((root / "dist").glob("peptacular-*.whl"), key=lambda path: path.stat().st_mtime)
if not wheels:
    raise SystemExit("Build a wheel with uv build before running this check.")
with TemporaryDirectory(prefix="peptacular-wheel-") as temporary:
    env = Path(temporary) / "env"
    subprocess.run(["uv", "venv", "--python", sys.executable, str(env)], check=True)
    python = env / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    requirement = f"{wheels[-1]}[{args.extra}]" if args.extra else str(wheels[-1])
    subprocess.run(["uv", "pip", "install", "--python", str(python), requirement], check=True)
    smoke_args = ["--extra", args.extra] if args.extra else []
    subprocess.run([str(python), "-I", str(root / "scripts/wheel_smoke.py"), *smoke_args], cwd=temporary, check=True)
