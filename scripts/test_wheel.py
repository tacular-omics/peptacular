"""Install the freshly built wheel into a temporary environment and smoke test it."""

import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

root = Path(__file__).resolve().parent.parent
wheels = sorted((root / "dist").glob("peptacular-*.whl"), key=lambda path: path.stat().st_mtime)
if not wheels:
    raise SystemExit("Build a wheel with uv build before running this check.")
with TemporaryDirectory(prefix="peptacular-wheel-") as temporary:
    env = Path(temporary) / "env"
    subprocess.run(["uv", "venv", "--python", sys.executable, str(env)], check=True)
    python = env / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    subprocess.run(["uv", "pip", "install", "--python", str(python), str(wheels[-1])], check=True)
    subprocess.run([str(python), "-I", str(root / "scripts/wheel_smoke.py")], cwd=temporary, check=True)
