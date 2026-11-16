"""Fragmentation presets: the ions an activation method usually produces.

A preset is the keyword arguments :func:`peptacular.fragment` already takes, so naming an
activation method replaces writing the ion series, the neutral losses and the charge rule by
hand:

.. code-block:: python

    import peptacular as pt

    pt.fragment("PEPTIDEK", **pt.fragmentation_preset("hcd", charge_state=2))

**A preset is a default, not a measurement.** It names the series a method usually gives for
tryptic-sized peptides, taken from the literature cited for each entry. The real series depend on
the instrument, the activation energy, the precursor charge and the peptide itself, so a preset is
a starting point to override, never evidence about a particular spectrum.

Every preset sets ``ion_types``, ``neutral_deltas`` and ``max_ndeltas``, and
:func:`fragmentation_preset` adds ``charges`` when it is given a ``charge_state``. The charge rule
is the same in every preset, ``1`` to ``max(1, charge_state - 1)``, because the fragments of a
``z+`` precursor cannot each carry ``z`` charges.

=========  ===============================  ================  ==========================================
preset     ion_types                        neutral_deltas    primary literature
=========  ===============================  ================  ==========================================
cid        b, y                             H2O, NH3          Roepstorff and Fohlman 1984, Biol Mass
                                                              Spectrom 11:601 (10.1002/bms.1200111109),
                                                              the letter that names a, b, c, x, y, z;
                                                              Paizs and Suhai 2005, Mass Spectrom Rev
                                                              24:508-548 (10.1002/mas.20024) for b and y
                                                              dominating low-energy CID and for the
                                                              water and ammonia losses
hcd        b, y                             H2O, NH3          Olsen et al. 2007, Nat Methods 4:709-712
                                                              (10.1038/nmeth1060)
etd        c, z.                            H2O, NH3          Syka et al. 2004, PNAS 101:9528-9533
                                                              (10.1073/pnas.0402700101)
ecd        c, z.                            H2O, NH3          Zubarev et al. 1998, JACS 120:3265-3266
                                                              (10.1021/ja973478k)
ethcd      b, y, c, z.                      H2O, NH3          Frese et al. 2012, Anal Chem 84:9668-9673
                                                              (10.1021/ac3025366); Swaney et al. 2007,
                                                              Anal Chem 79:477-485 (10.1021/ac061457f)
                                                              for the b and y ions that collisional
                                                              activation of the electron-transfer
                                                              survivor adds
ead        a, b, c, x, y, z.                H2O, NH3          Baba et al. 2021, J Am Soc Mass Spectrom
                                                              32:1964-1975 (10.1021/jasms.0c00425)
uvpd       a, b, c, x, y, z                 H2O               Shaw et al. 2013, JACS 135:12646-12651
                                                              (10.1021/ja4029654); Brodbelt 2020, Chem
                                                              Rev 120:3328-3380
                                                              (10.1021/acs.chemrev.9b00440)
=========  ===============================  ================  ==========================================

``etd`` holds no ``b`` or ``y`` series. The b and y ions seen in an electron-transfer experiment
come from collisional activation of the charge-reduced survivor (Swaney et al. 2007, and Han et al.
2007, Rapid Commun Mass Spectrom 21:1567-1573, 10.1002/rcm.2994), which is what ``ethcd`` is for.

A preset for a method not listed here is a dict: write it into :data:`FRAGMENTATION_PRESETS` and
:func:`fragmentation_preset` will find it.

``fast_fragment`` takes only ``ion_types``, ``charges`` and ``monoisotopic``, so a preset cannot be
splatted into it. Filter the preset first:

.. code-block:: python

    preset = pt.fragmentation_preset("hcd", charge_state=2)
    fast = {k: v for k, v in preset.items() if k in ("ion_types", "charges", "monoisotopic")}
    pt.fast_fragment("PEPTIDEK", **fast)
"""

from __future__ import annotations

from typing import Any

from .diagnostics import PeptacularError

__all__ = ["FRAGMENTATION_PRESETS", "FRAGMENTATION_PRESET_ALIASES", "FRAGMENTATION_PRESET_CITATIONS", "fragmentation_preset"]

_WATER_AND_AMMONIA = ("H2O", "NH3")

FRAGMENTATION_PRESETS: dict[str, dict[str, Any]] = {
    "cid": {"ion_types": ("b", "y"), "neutral_deltas": _WATER_AND_AMMONIA, "max_ndeltas": 1},
    "hcd": {"ion_types": ("b", "y"), "neutral_deltas": _WATER_AND_AMMONIA, "max_ndeltas": 1},
    "etd": {"ion_types": ("c", "z."), "neutral_deltas": _WATER_AND_AMMONIA, "max_ndeltas": 1},
    "ecd": {"ion_types": ("c", "z."), "neutral_deltas": _WATER_AND_AMMONIA, "max_ndeltas": 1},
    "ethcd": {"ion_types": ("b", "y", "c", "z."), "neutral_deltas": _WATER_AND_AMMONIA, "max_ndeltas": 1},
    "ead": {"ion_types": ("a", "b", "c", "x", "y", "z."), "neutral_deltas": _WATER_AND_AMMONIA, "max_ndeltas": 1},
    "uvpd": {"ion_types": ("a", "b", "c", "x", "y", "z"), "neutral_deltas": ("H2O",), "max_ndeltas": 1},
}
"""The keyword arguments of :func:`peptacular.fragment` for each activation method, except
``charges``, which :func:`fragmentation_preset` adds from the precursor charge. Writable: a new
method is a new key."""

FRAGMENTATION_PRESET_ALIASES: dict[str, str] = {
    "cad": "cid",
    "beam-type cid": "hcd",
    "pasef": "hcd",
    "etcid": "ethcd",
    "etcad": "ethcd",
}
"""Other names for a preset. ``pasef`` fragments by beam-type collisional dissociation, so it maps
to ``hcd``; ``etcid`` and ``etcad`` are electron transfer with supplemental collisional activation,
the same series as ``ethcd``."""

FRAGMENTATION_PRESET_CITATIONS: dict[str, str] = {
    "cid": (
        "Roepstorff and Fohlman 1984, Biol Mass Spectrom 11:601 (10.1002/bms.1200111109); "
        "Paizs and Suhai 2005, Mass Spectrom Rev 24:508-548 (10.1002/mas.20024)"
    ),
    "hcd": "Olsen et al. 2007, Nat Methods 4:709-712 (10.1038/nmeth1060)",
    "etd": "Syka et al. 2004, PNAS 101:9528-9533 (10.1073/pnas.0402700101)",
    "ecd": "Zubarev et al. 1998, JACS 120:3265-3266 (10.1021/ja973478k)",
    "ethcd": ("Frese et al. 2012, Anal Chem 84:9668-9673 (10.1021/ac3025366); Swaney et al. 2007, Anal Chem 79:477-485 (10.1021/ac061457f)"),
    "ead": "Baba et al. 2021, J Am Soc Mass Spectrom 32:1964-1975 (10.1021/jasms.0c00425)",
    "uvpd": ("Shaw et al. 2013, JACS 135:12646-12651 (10.1021/ja4029654); Brodbelt 2020, Chem Rev 120:3328-3380 (10.1021/acs.chemrev.9b00440)"),
}
"""The primary literature each preset rests on, one entry per key of
:data:`FRAGMENTATION_PRESETS`."""


def fragmentation_preset(method: str, *, charge_state: int | None = None) -> dict[str, Any]:
    """The :func:`peptacular.fragment` keyword arguments for an activation method.

    A preset is a default, not a measurement: it names the series the method usually gives, from
    the literature in :data:`FRAGMENTATION_PRESET_CITATIONS`, and the real series depend on the
    instrument, the energy and the peptide.

    The returned dict is a fresh copy, so overriding one entry is a dict literal and the caller's
    key wins.

    :param method: A key of :data:`FRAGMENTATION_PRESETS` or of
        :data:`FRAGMENTATION_PRESET_ALIASES`, in any case and with surrounding spaces ignored.
    :type method: str
    :param charge_state: The precursor charge. When given, the preset carries
        ``charges=(1, ..., max(1, abs(charge_state) - 1))``. When None, no ``charges`` key is set
        and :func:`peptacular.fragment` uses the annotation's own charge.
    :type charge_state: int | None
    :return: The keyword arguments, a new dict on every call.
    :rtype: dict[str, Any]
    :raises PeptacularError: If ``method`` names no preset. The message lists the valid names.

    >>> import peptacular as pt
    >>> pt.fragmentation_preset("hcd")
    {'ion_types': ('b', 'y'), 'neutral_deltas': ('H2O', 'NH3'), 'max_ndeltas': 1}
    >>> pt.fragmentation_preset("etd", charge_state=3)["charges"]
    (1, 2)
    >>> len(pt.fragment("PEPTIDEK", **pt.fragmentation_preset("hcd", charge_state=2)))
    39
    >>> hcd = pt.fragmentation_preset("hcd", charge_state=2)
    >>> len(pt.fragment("PEPS[Phospho]TIDEK", **{**hcd, "neutral_deltas": ("H3PO4",)}))
    28
    """

    key = method.strip().lower()
    key = FRAGMENTATION_PRESET_ALIASES.get(key, key)
    if key not in FRAGMENTATION_PRESETS:
        names = ", ".join(sorted(FRAGMENTATION_PRESETS))
        aliases = ", ".join(sorted(FRAGMENTATION_PRESET_ALIASES))
        raise PeptacularError(f"Unknown activation method: {method!r}. Presets: {names}. Aliases: {aliases}")
    preset = dict(FRAGMENTATION_PRESETS[key])
    if charge_state is not None:
        preset["charges"] = tuple(range(1, max(1, abs(charge_state) - 1) + 1))
    return preset
