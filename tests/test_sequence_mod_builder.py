"""Coverage for the functional modification-builder API (peptacular.sequence.mod_builder)."""

import pytest

import peptacular as pt
from peptacular.sequence.mod_builder import (
    append_mods,
    condense_static_mods,
    condense_to_peptidoform,
    extend_mods,
    filter_mods,
    from_ms2_pip,
    get_mods,
    modify,
    pop_mods,
    remove_mods,
    set_mods,
    strip_mods,
    to_ms2_pip,
)

SEQ = "PEPTIDE"


class TestModify:
    def test_scalar(self):
        result = modify(SEQ, internal_variable={"P": [79.966]}, max_variable_mods=1)
        assert "PEPTIDE" in result
        assert any("79.966" in r for r in result)

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert modify(a, internal_variable={"P": [79.966]}, max_variable_mods=1) == modify(SEQ, internal_variable={"P": [79.966]}, max_variable_mods=1)

    def test_batch(self):
        result = modify([SEQ, "PROTEIN"], internal_variable={"P": [79.966]}, max_variable_mods=1)
        assert len(result) == 2

    def test_batch_with_parallel_kwargs(self):
        result = modify([SEQ, "PROTEIN"], internal_variable={"P": [79.966]}, max_variable_mods=1, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_unique_peptidoforms(self):
        result = modify(SEQ, internal_variable={"P": [79.966]}, max_variable_mods=1, unique_peptidoforms=True)
        assert len(result) == len(set(result))

    def test_static_mods_kwargs(self):
        result = modify(SEQ, internal_static={"P": [10.0]})
        assert len(result) == 1

    def test_terminal_mod_kwargs(self):
        result = modify(SEQ, nterm_static={None: [10.0]}, cterm_variable={None: [20.0]}, max_variable_mods=1)
        assert len(result) >= 1


class TestGetMods:
    def test_returns_dict(self):
        result = get_mods("PEP[Oxidation]TIDE")
        assert isinstance(result, dict)

    def test_annotation_input(self):
        a = pt.parse("PEP[Oxidation]TIDE")
        assert get_mods(a) == get_mods("PEP[Oxidation]TIDE")


class TestSetMods:
    def test_sets_internal_mod(self):
        assert set_mods(SEQ, {3: "Oxidation"}) == "PEPT[Oxidation]IDE"


class TestAppendMods:
    def test_appends_to_existing(self):
        assert append_mods("PEP[Oxidation]TIDE", {3: "Amide"}) == "PEP[Oxidation]T[Amide]IDE"


class TestExtendMods:
    def test_extends_with_list(self):
        assert extend_mods(SEQ, {3: ["Oxidation", "Amide"]}) == "PEPT[Oxidation][Amide]IDE"


class TestCondenseStaticMods:
    def test_condenses(self):
        assert condense_static_mods("<13C><[100]@P>PEPTIDE") == "<13C>P[100]EP[100]TIDE"

    def test_unmodified_passthrough(self):
        assert condense_static_mods(SEQ) == SEQ


class TestPopMods:
    def test_returns_seq_and_dict(self):
        seq, mod_dict = pop_mods("PEP[phospho]TIDE")
        assert seq == SEQ
        assert isinstance(mod_dict, dict)


class TestRemoveMods:
    def test_strips_mods(self):
        assert remove_mods("PEP[phospho]TIDE") == SEQ


class TestStripMods:
    def test_scalar(self):
        assert strip_mods("PEP[phospho]TIDE") == SEQ

    def test_scalar_annotation(self):
        a = pt.parse("PEP[phospho]TIDE")
        assert strip_mods(a) == SEQ

    def test_batch(self):
        assert strip_mods(["PEP[phospho]TIDE", "PEP[phospho]TIDE"]) == [SEQ, SEQ]

    def test_mods_filter_kwarg(self):
        assert strip_mods("PEP[phospho]TIDE", mods="internal") == SEQ


class TestFilterMods:
    def test_keeps_specified(self):
        assert filter_mods("PEP[phospho]TIDE", mods="internal") == "PEP[phospho]TIDE"


class TestToMs2Pip:
    def test_scalar(self):
        assert to_ms2_pip("PEP[Phospho]TIDE") == (SEQ, "3|Phospho")

    def test_scalar_annotation(self):
        a = pt.parse("PEP[Phospho]TIDE")
        assert to_ms2_pip(a) == to_ms2_pip("PEP[Phospho]TIDE")

    def test_batch(self):
        result = to_ms2_pip(["PEP[Phospho]TIDE", "PROT[Oxidation]EIN"])
        assert result == [(SEQ, "3|Phospho"), ("PROTEIN", "4|Oxidation")]

    def test_batch_with_parallel_kwargs(self):
        result = to_ms2_pip(["PEP[Phospho]TIDE", "PEP[Phospho]TIDE"], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2


class TestFromMs2Pip:
    def test_scalar(self):
        assert from_ms2_pip((SEQ, "3|Phospho")) == "PEP[Phospho]TIDE"

    def test_scalar_empty_mods(self):
        assert from_ms2_pip((SEQ, "")) == SEQ

    def test_batch(self):
        items = [(SEQ, "3|Phospho"), ("PROTEIN", "4|Oxidation")]
        assert from_ms2_pip(items) == ["PEP[Phospho]TIDE", "PROT[Oxidation]EIN"]

    def test_batch_with_parallel_kwargs(self):
        items = [(SEQ, "3|Phospho"), (SEQ, "3|Phospho")]
        result = from_ms2_pip(items, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_static_mods_kwarg(self):
        result = from_ms2_pip((SEQ, ""), static_mods={"C": 57.021})
        assert result == SEQ

    def test_batch_with_non_tuple_item_raises(self):
        with pytest.raises(ValueError, match="must be tuples"):
            from_ms2_pip([(SEQ, "3|Phospho"), "not_a_tuple"])

    def test_scalar_non_tuple_string_raises(self):
        # A bare string is itself a Sequence, so it is treated as a batch of
        # single characters and fails the "all items are tuples" check.
        with pytest.raises(ValueError, match="must be tuples"):
            from_ms2_pip("PEPTIDE")

    def test_scalar_wrong_length_tuple_raises(self):
        with pytest.raises(ValueError, match="for single processing"):
            from_ms2_pip(("a", "b", "c"))


class TestCondenseToPeptidoform:
    def test_scalar(self):
        result = condense_to_peptidoform("PEP[Oxidation]TIDE")
        assert result.startswith("[Oxidation]?")

    def test_scalar_annotation(self):
        a = pt.parse("PEP[Oxidation]TIDE")
        assert condense_to_peptidoform(a) == condense_to_peptidoform("PEP[Oxidation]TIDE")

    def test_batch(self):
        result = condense_to_peptidoform(["PEP[Oxidation]TIDE", "PEP[Oxidation]TIDE"])
        assert len(result) == 2
