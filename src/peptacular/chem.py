from collections import Counter
from collections.abc import Mapping, Sequence
from typing import overload

from tacular import ELEMENT_LOOKUP, ElementInfo

from .proforma_components import ChargedFormula
from .sequence.parallel import (
    parallel_apply_internal,
    parallelMethod,
    parallelMethodLiteral,
)


def _parse_formula_single(formula: str | Mapping[ElementInfo | str, int], sep: str = "") -> Counter[ElementInfo]:
    if isinstance(formula, str):
        return ChargedFormula.from_string(formula, allow_zero=False, require_formula_prefix=False, sep=sep).get_composition()

    composition: Counter[ElementInfo] = Counter[ElementInfo]()
    for key, value in formula.items():
        if isinstance(key, ElementInfo):
            composition[key] += value
        else:
            composition[ELEMENT_LOOKUP[key]] += value

    return composition


@overload
def parse_formula(
    formula: str | Mapping[ElementInfo | str, int],
    sep: str = "",
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> Counter[ElementInfo]: ...


@overload
def parse_formula(
    formula: Sequence[str | Mapping[ElementInfo | str, int]],
    sep: str = "",
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[Counter[ElementInfo]]: ...


def parse_formula(
    formula: str | Mapping[ElementInfo | str, int] | Sequence[str | Mapping[ElementInfo | str, int]],
    sep: str = "",
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> Counter[ElementInfo] | list[Counter[ElementInfo]]:
    """Parse a chemical formula string or list of formulas into elemental composition."""
    if isinstance(formula, Sequence) and not isinstance(formula, str):
        return parallel_apply_internal(
            _parse_formula_single,
            formula,
            sep=sep,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            reuse_pool=reuse_pool,
        )
    else:
        return _parse_formula_single(formula, sep=sep)


def _chem_comp_single(
    formula: str | Mapping[ElementInfo | str, int],
) -> Counter[ElementInfo]:
    return _parse_formula_single(formula)


@overload
def chem_comp(
    formula: str | Mapping[ElementInfo | str, int],
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> Counter[ElementInfo]: ...


@overload
def chem_comp(
    formula: Sequence[str | Mapping[ElementInfo | str, int]],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[Counter[ElementInfo]]: ...


def chem_comp(
    formula: str | Mapping[ElementInfo | str, int] | Sequence[str | Mapping[ElementInfo | str, int]],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> Counter[ElementInfo] | list[Counter[ElementInfo]]:
    """Get the elemental composition of a chemical formula or list of formulas."""
    if isinstance(formula, Sequence) and not isinstance(formula, str):
        return parallel_apply_internal(
            _chem_comp_single,
            formula,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            reuse_pool=reuse_pool,
        )
    else:
        return _chem_comp_single(formula)


def _chem_mass_single(
    formula: str | Mapping[ElementInfo | str, int],
    monoisotopic: bool = True,
) -> float:
    comp = _parse_formula_single(formula)
    mass = 0.0
    for element, count in comp.items():
        mass += element.get_mass(monoisotopic=monoisotopic) * count
    return mass


@overload
def chem_mass(
    formula: str | Mapping[ElementInfo | str, int],
    monoisotopic: bool = True,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> float: ...


@overload
def chem_mass(
    formula: Sequence[str | Mapping[ElementInfo | str, int]],
    monoisotopic: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[float]: ...


def chem_mass(
    formula: str | Mapping[ElementInfo | str, int] | Sequence[str | Mapping[ElementInfo | str, int]],
    monoisotopic: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> float | list[float]:
    """Calculate the mass of a chemical formula or list of formulas."""
    if isinstance(formula, Sequence) and not isinstance(formula, str):
        return parallel_apply_internal(
            _chem_mass_single,
            formula,
            monoisotopic=monoisotopic,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            reuse_pool=reuse_pool,
        )
    else:
        return _chem_mass_single(formula, monoisotopic=monoisotopic)


def _chem_formula_single(
    comp: str | Mapping[ElementInfo | str, int],
    hill_order: bool = True,
    sep: str = "",
    include_formula_prefix: bool = False,
) -> str:
    """Generate a chemical formula string from an elemental composition."""
    if isinstance(comp, str):
        comp: Counter[ElementInfo] = _parse_formula_single(comp, sep=sep)
    return ChargedFormula.from_composition(comp, charge=None).serialize(hill_order=hill_order, sep=sep, include_formula_prefix=include_formula_prefix)


@overload
def chem_formula(
    comp: str | Mapping[ElementInfo | str, int],
    hill_order: bool = True,
    sep: str = "",
    include_formula_prefix: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> str: ...


@overload
def chem_formula(
    comp: str | Sequence[Mapping[ElementInfo | str, int]],
    hill_order: bool = True,
    sep: str = "",
    include_formula_prefix: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[str]: ...


def chem_formula(
    comp: str | Mapping[ElementInfo | str, int] | Sequence[Mapping[ElementInfo | str, int]],
    hill_order: bool = True,
    sep: str = "",
    include_formula_prefix: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> str | list[str]:
    """Generate a chemical formula string from an elemental composition or list of compositions."""
    if isinstance(comp, Sequence) and not isinstance(comp, str):
        return parallel_apply_internal(
            _chem_formula_single,
            comp,
            hill_order=hill_order,
            sep=sep,
            include_formula_prefix=include_formula_prefix,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            reuse_pool=reuse_pool,
        )
    else:
        return _chem_formula_single(
            comp,
            hill_order=hill_order,
            sep=sep,
            include_formula_prefix=include_formula_prefix,
        )
