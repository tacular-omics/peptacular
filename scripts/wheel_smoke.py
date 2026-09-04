"""Exercise the installed distribution without importing from the source tree."""

import io
from pathlib import Path

import peptacular as pt

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
print(f"Installed wheel smoke checks passed: {package}")
