"""Coverage for the functional sequence-basic API (peptacular.sequence.basic):
parse_chimeric, serialize_chimeric, parse, serialize, sequence_length, is_ambiguous,
is_modified, count_residues, percent_residues, annotate_ambiguity, validate, generate_random.
"""

import pytest

import peptacular as pt
from peptacular.sequence.basic import (
    annotate_ambiguity,
    count_residues,
    generate_random,
    is_ambiguous,
    is_modified,
    parse,
    parse_chimeric,
    percent_residues,
    sequence_length,
    serialize,
    serialize_chimeric,
    validate,
)

SEQ = "PEPTIDE"


class TestParseChimeric:
    def test_scalar(self):
        result = parse_chimeric("PEPTIDE+SEQUENCE")
        assert len(result) == 2
        assert all(isinstance(a, pt.ProFormaAnnotation) for a in result)

    def test_batch_matches_scalar(self):
        scalar = parse_chimeric("PEPTIDE+SEQUENCE")
        batch = parse_chimeric(["PEPTIDE+SEQUENCE", "PEPTIDE+SEQUENCE"])
        assert len(batch) == 2
        assert [a.serialize() for a in batch[0]] == [a.serialize() for a in scalar]

    def test_batch_with_parallel_kwargs(self):
        result = parse_chimeric(["PEPTIDE+SEQUENCE", "PEPTIDE+SEQUENCE"], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_validate_kwarg(self):
        result = parse_chimeric("PEPTIDE+SEQUENCE", validate=True)
        assert len(result) == 2


class TestSerializeChimeric:
    def test_scalar_str_list(self):
        assert serialize_chimeric(["PEPTIDE", "SEQUENCE"]) == "PEPTIDE+SEQUENCE"

    def test_scalar_annotation_list(self):
        annots = parse_chimeric("PEPTIDE+SEQUENCE")
        assert serialize_chimeric(annots) == "PEPTIDE+SEQUENCE"

    def test_batch(self):
        result = serialize_chimeric([["PEPTIDE", "SEQUENCE"], ["PEPTIDE", "SEQUENCE"]])
        assert result == ["PEPTIDE+SEQUENCE", "PEPTIDE+SEQUENCE"]

    def test_batch_with_parallel_kwargs(self):
        result = serialize_chimeric(
            [["PEPTIDE", "SEQUENCE"], ["PEPTIDE", "SEQUENCE"]], n_workers=1, chunksize=1, method="sequential"
        )
        assert result == ["PEPTIDE+SEQUENCE", "PEPTIDE+SEQUENCE"]

    def test_compound_name_mismatch_raises(self):
        a = pt.parse("PEPTIDE").copy()
        a.compound_name = "foo"
        b = pt.parse("PEPTIDE").copy()
        b.compound_name = "bar"
        with pytest.raises(ValueError, match="compound name"):
            serialize_chimeric([a, b])

    def test_static_mods_mismatch_raises(self):
        a = pt.parse("<[Carbamidomethyl]@C>PEPTIDE")
        b = pt.parse("PEPTIDE")
        with pytest.raises(ValueError, match="static modifications"):
            serialize_chimeric([a, b])

    def test_isotope_mods_mismatch_raises(self):
        a = pt.parse("<13C>PEPTIDE")
        b = pt.parse("PEPTIDE")
        with pytest.raises(ValueError, match="isotopic modifications"):
            serialize_chimeric([a, b])


class TestParse:
    def test_scalar(self):
        result = parse(SEQ)
        assert isinstance(result, pt.ProFormaAnnotation)
        assert result.serialize() == SEQ

    def test_batch(self):
        result = parse([SEQ, SEQ])
        assert len(result) == 2

    def test_batch_with_parallel_kwargs(self):
        result = parse([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_validate_kwarg(self):
        assert parse(SEQ, validate=True).serialize() == SEQ


class TestSerialize:
    def test_scalar_string(self):
        assert serialize(SEQ) == SEQ

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert serialize(a) == SEQ

    def test_batch(self):
        assert serialize([SEQ, SEQ]) == [SEQ, SEQ]

    def test_batch_with_parallel_kwargs(self):
        assert serialize([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential") == [SEQ, SEQ]


class TestSequenceLength:
    def test_scalar_string(self):
        assert sequence_length(SEQ) == 7

    def test_scalar_annotation(self):
        assert sequence_length(pt.parse(SEQ)) == 7

    def test_batch(self):
        assert sequence_length([SEQ, SEQ]) == [7, 7]

    def test_batch_with_parallel_kwargs(self):
        assert sequence_length([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential") == [7, 7]


class TestIsAmbiguous:
    def test_scalar_false(self):
        assert is_ambiguous(SEQ) is False

    def test_scalar_annotation(self):
        assert is_ambiguous(pt.parse(SEQ)) is False

    def test_scalar_true(self):
        assert is_ambiguous("(?PE)PTIDE") is True

    def test_batch(self):
        assert is_ambiguous([SEQ, "(?PE)PTIDE"]) == [False, True]

    def test_batch_with_parallel_kwargs(self):
        assert is_ambiguous([SEQ, "(?PE)PTIDE"], n_workers=1, chunksize=1, method="sequential") == [False, True]


class TestIsModified:
    def test_scalar_false(self):
        assert is_modified(SEQ) is False

    def test_scalar_annotation(self):
        assert is_modified(pt.parse(SEQ)) is False

    def test_scalar_true(self):
        assert is_modified("PEP[Oxidation]TIDE") is True

    def test_batch(self):
        assert is_modified([SEQ, "PEP[Oxidation]TIDE"]) == [False, True]

    def test_batch_with_parallel_kwargs(self):
        assert is_modified([SEQ, "PEP[Oxidation]TIDE"], n_workers=1, chunksize=1, method="sequential") == [False, True]


class TestCountResidues:
    def test_scalar(self):
        assert count_residues(SEQ) == {"P": 2, "E": 2, "T": 1, "I": 1, "D": 1}

    def test_scalar_annotation(self):
        assert count_residues(pt.parse(SEQ)) == {"P": 2, "E": 2, "T": 1, "I": 1, "D": 1}

    def test_include_mods_false(self):
        result = count_residues("PEP[Oxidation]TIDE", include_mods=False)
        assert result == {"P": 2, "E": 2, "T": 1, "I": 1, "D": 1}

    def test_batch(self):
        assert count_residues([SEQ, SEQ]) == [{"P": 2, "E": 2, "T": 1, "I": 1, "D": 1}] * 2

    def test_batch_with_parallel_kwargs(self):
        result = count_residues([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert result == [{"P": 2, "E": 2, "T": 1, "I": 1, "D": 1}] * 2

    def test_does_not_mutate_input_annotation(self):
        a = pt.parse("<[Carbamidomethyl]@C>PEPTIDE")
        before = a.serialize()
        count_residues(a)
        assert a.serialize() == before


class TestPercentResidues:
    def test_scalar(self):
        result = percent_residues(SEQ)
        assert round(result["P"], 2) == 28.57

    def test_scalar_annotation(self):
        result = percent_residues(pt.parse(SEQ))
        assert round(result["P"], 2) == 28.57

    def test_include_mods_false(self):
        result = percent_residues("PEP[Oxidation]TIDE", include_mods=False)
        assert "P[Oxidation]" not in result

    def test_batch(self):
        result = percent_residues([SEQ, SEQ])
        assert len(result) == 2

    def test_batch_with_parallel_kwargs(self):
        result = percent_residues([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2


class TestAnnotateAmbiguity:
    def test_basic(self):
        result = annotate_ambiguity(SEQ, [0, 1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0])
        assert result == "(?PE)PTI(?DE)"

    def test_with_mass_shift(self):
        result = annotate_ambiguity(SEQ, [1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1], 79.966)
        assert result == "PEPT[+79.966]IDE"

    def test_condense_to_xnotation(self):
        result = annotate_ambiguity(SEQ, [0, 1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0], condense_to_xnotation=True)
        assert result.startswith("X[")
        assert "?" not in result


class TestValidate:
    def test_scalar_valid(self):
        assert validate(SEQ) is True

    def test_scalar_invalid(self):
        assert validate("PEP[Oxidation") is False

    def test_scalar_annotation(self):
        assert validate(pt.parse(SEQ)) is True

    def test_batch(self):
        assert validate([SEQ, "PEP[Oxidation"]) == [True, False]

    def test_batch_with_parallel_kwargs(self):
        assert validate([SEQ, "PEP[Oxidation"], n_workers=1, chunksize=1, method="sequential") == [True, False]


class TestGenerateRandom:
    def test_single(self):
        result = generate_random()
        assert isinstance(result, pt.ProFormaAnnotation)

    def test_count(self):
        result = generate_random(count=3)
        assert len(result) == 3

    def test_count_with_parallel_kwargs(self):
        result = generate_random(count=3, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 3

    def test_feature_flags(self):
        result = generate_random(
            mod_probability=0.0,
            include_internal_mods=False,
            include_nterm_mods=False,
            include_cterm_mods=False,
            include_labile_mods=False,
            include_unknown_mods=False,
            include_isotopic_mods=False,
            include_static_mods=False,
            include_intervals=False,
            include_charge=False,
            require_composition=False,
        )
        assert isinstance(result, pt.ProFormaAnnotation)
