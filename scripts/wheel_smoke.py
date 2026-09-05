"""Exercise the installed distribution without importing from the source tree."""

import argparse
import importlib.util
import io
import sys
from pathlib import Path

import peptacular as pt
from peptacular.interop import (
    MissingOptionalDependencyError,
    from_alphabase_dataframe,
    from_psm_utils,
    from_pyteomics,
    to_alphabase_dataframe,
    to_alphabase_row,
    to_psm_utils,
    to_pyteomics,
    to_pyteomics_composition,
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--extra", choices=["pyteomics", "psm-utils", "alphabase", "mcp", "mcp,pyteomics", "mcp,psm-utils", "mcp,alphabase"])
args = parser.parse_args()
package = Path(pt.__file__).resolve()
assert "site-packages" in package.parts, package
assert package.with_name("py.typed").is_file()
assert abs(pt.mass("PEPTIDE") - 799.3599643488) < 1e-5
assert pt.parse("PEM[Oxidation]TIDE").serialize() == "PEM[Oxidation]TIDE"
assert pt.fragment("PEPTIDE", ion_types=["b", "y"], charges=[1])
assert pt.digest("PEPTIDEKS", pt.Proteases.TRYPSIN, missed_cleavages=0)
assert list(pt.iter_fasta(io.StringIO(">protein\nPEPTIDE")))[0].sequence == "PEPTIDE"
results = pt.batch("mass", ["PEPTIDE", "PEP[UnknownModification]TIDE"], errors="collect")
assert results[0].ok and results[1].error.code == "unresolved_modification"
assert pt.diagnose("PEP[+42]TIDE", "comp").code == "unavailable_composition"
annotation = pt.parse("PEM[Oxidation]TIDE/2")
assert pt.ProFormaAnnotation.from_json(annotation.to_json()).to_dict() == annotation.to_dict()
assert pt.get_proforma_json_schema()["$id"] == pt.PROFORMA_JSON_SCHEMA_ID

extras = set((args.extra or "").split(","))
assert "mcp" not in sys.modules

if "pyteomics" in extras:
    assert from_pyteomics(to_pyteomics(annotation)) == annotation
    assert dict(to_pyteomics_composition({"C": 2, "13C": 1})) == {"C": 2, "C[13]": 1}
elif "psm-utils" in extras:
    assert from_psm_utils(to_psm_utils(annotation)) == annotation
elif "alphabase" in extras:
    assert from_alphabase_dataframe(to_alphabase_dataframe([annotation])) == [annotation]
elif "mcp" not in extras:
    for optional in ("alphabase", "psm_utils", "pyteomics", "pandas", "mcp", "pydantic"):
        assert importlib.util.find_spec(optional) is None, optional
        assert optional not in sys.modules, optional
    try:
        to_alphabase_row(annotation)
    except MissingOptionalDependencyError:
        pass
    else:
        raise AssertionError("A core installation must not include optional integrations")
print(f"Installed wheel smoke checks passed ({args.extra or 'core'}): {package}")

if "mcp" in extras:
    import asyncio

    from mcp import Client
    from mcp.client.stdio import StdioServerParameters

    async def check_mcp():
        parameters = StdioServerParameters(
            command=sys.executable,
            args=["-I", "-m", "peptacular.mcp", "--workspace", str(Path.cwd()), "--cache", str(Path.cwd() / "cache")],
        )
        async with Client(parameters) as client:
            tools = await client.list_tools()
            assert len(tools.tools) == 18
            result = await client.call_tool(
                "analyze_peptides",
                {
                    "request": {"inputs": {"kind": "inline", "records": [{"annotation": "PEPTIDE/2"}]}, "measurements": ["mz"]},
                },
            )
            assert not result.is_error
            assert result.structured_content["records"][0]["mz"] > 400

    asyncio.run(check_mcp())
    print("Installed MCP stdio smoke check passed")
