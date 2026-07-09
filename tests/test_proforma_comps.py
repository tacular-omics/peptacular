"""Coverage for peptacular.proforma_components.comps: the ProForma component dataclasses.

These are the frozen dataclasses that back the ProForma 2.1 object model (formulas,
modification tags, sequence elements, peptidoforms, etc). Tests exercise real lookups
(UNIMOD/PSI-MOD/RESID/GNO/XLMOD via `tacular`) rather than mocking them, per project
convention.
"""

from collections import Counter

import pytest
from tacular import AA_LOOKUP, ELEMENT_LOOKUP, AminoAcid, Element, Monosaccharide

from peptacular.constants import CV, Terminal
from peptacular.proforma_components.comps import (
    ChargedFormula,
    ComkpTag,
    CompoundPeptidoformIon,
    ComupTag,
    FixedModification,
    FormulaElement,
    GlobalChargeCarrier,
    GlycanComponent,
    GlycanTag,
    IsotopeReplacement,
    LimitTag,
    ModificationAmbiguousPrimary,
    ModificationAmbiguousSecondary,
    ModificationCrossLinker,
    ModificationTags,
    Peptidoform,
    PeptidoformIon,
    PositionRule,
    PositionScore,
    PositionTag,
    SequenceElement,
    SequenceRegion,
    TagAccession,
    TagCustom,
    TagInfo,
    TagMass,
    TagName,
)

# Real accessions/names with no defined elemental composition in their CVs, used to
# exercise the "found but missing composition" defensive branches with real data.
PSIMOD_NO_COMPOSITION_ID = "02028"
PSIMOD_NO_COMPOSITION_NAME = "iTRAQ4plex reporter+balance reagent acylated residue, average mass modification"


class TestMassPropertyMixin:
    """FormulaElement is a concrete MassPropertyMixin subclass; use it to exercise the mixin."""

    def test_monoisotopic_mass_property(self) -> None:
        fe = FormulaElement(element=Element.C, occurance=2)
        assert fe.monoisotopic_mass == fe.get_mass(monoisotopic=True)

    def test_average_mass_property(self) -> None:
        fe = FormulaElement(element=Element.C, occurance=2)
        assert fe.average_mass == fe.get_mass(monoisotopic=False)


class TestPositionScoreMixin:
    """Exercised through TagAccession, a concrete PositionScoreMixin subclass."""

    def test_both_position_and_score(self) -> None:
        tag = TagAccession(accession="35", cv=CV.UNIMOD, position_id="1", score=0.9)
        assert tag.serialize_position_score() == "#1(0.9)"

    def test_position_only(self) -> None:
        tag = TagAccession(accession="35", cv=CV.UNIMOD, position_id="1")
        assert tag.serialize_position_score() == "#1"

    def test_score_only(self) -> None:
        tag = TagAccession(accession="35", cv=CV.UNIMOD, score=0.5)
        assert tag.serialize_position_score() == "(0.5)"

    def test_neither(self) -> None:
        tag = TagAccession(accession="35", cv=CV.UNIMOD)
        assert tag.serialize_position_score() == ""


class TestFormulaElement:
    def test_is_valid_true_for_known_isotope(self) -> None:
        assert FormulaElement(element=Element.C, occurance=2).is_valid is True

    def test_is_valid_false_for_unknown_isotope(self) -> None:
        assert FormulaElement(element=Element.C, occurance=2, isotope=999).is_valid is False

    def test_validate_returns_none_when_valid(self) -> None:
        assert FormulaElement(element=Element.C, occurance=2).validate() is None

    def test_validate_returns_error_message_when_invalid(self) -> None:
        error = FormulaElement(element=Element.C, occurance=2, isotope=999).validate()
        assert error is not None
        assert "999" in error

    def test_from_element_info_wraps_isotopes_in_brackets(self) -> None:
        info = ELEMENT_LOOKUP[(Element.C, 13)]
        fe = FormulaElement.from_element_info(info, 2)
        assert fe.serialize() == "[13C2]"

    def test_from_element_info_no_brackets_for_natural_element(self) -> None:
        info = ELEMENT_LOOKUP[(Element.C, None)]
        fe = FormulaElement.from_element_info(info, 1)
        assert fe.serialize() == "C"

    def test_get_composition_returns_element_and_count(self) -> None:
        fe = FormulaElement(element=Element.C, occurance=3)
        elem_info, count = fe.get_element_count()
        assert fe.get_composition() == Counter({elem_info: count})


class TestChargedFormula:
    def test_is_valid_true(self) -> None:
        cf = ChargedFormula.from_string("H2O", require_formula_prefix=False)
        assert cf.is_valid is True

    def test_is_valid_false_for_bad_element(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1, isotope=999),))
        assert cf.is_valid is False
        assert cf.validate() is not None

    def test_formula_dict(self) -> None:
        cf = ChargedFormula.from_string("H2O", require_formula_prefix=False)
        assert cf.formula_dict() == {"H": 2, "O": 1}

    def test_get_dict_composition(self) -> None:
        cf = ChargedFormula.from_string("H2O", require_formula_prefix=False)
        assert cf.get_dict_composition() == {"H": 2, "O": 1}

    def test_to_mz_paf_positive(self) -> None:
        # mzPAF's chemical-formula notation reuses ProForma's own notation:
        # atom then count ("C2"), not a reversed count-before-atom form.
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=2),))
        assert cf.to_mz_paf() == "+C2"

    def test_to_mz_paf_negative(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=-2),))
        assert cf.to_mz_paf() == "-C2"

    def test_to_mz_paf_omits_count_of_one(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1),))
        assert cf.to_mz_paf() == "+C"

    def test_to_mz_paf_multi_element_matches_spec_water_loss_example(self) -> None:
        # mzPAF spec section 4.5 gives "-H2O" as the canonical water-loss example.
        cf = ChargedFormula(
            formula=(
                FormulaElement(element=Element.H, occurance=2),
                FormulaElement(element=Element.O, occurance=1),
            )
        )
        assert cf.to_mz_paf() == "+H2O"

    def test_to_mz_paf_raises_on_mixed_signs(self) -> None:
        cf = ChargedFormula(
            formula=(
                FormulaElement(element=Element.C, occurance=2),
                FormulaElement(element=Element.H, occurance=-2),
            )
        )
        with pytest.raises(ValueError, match="both positive and negative"):
            cf.to_mz_paf()

    def test_to_mz_paf_raises_when_empty(self) -> None:
        with pytest.raises(ValueError, match="no elements present"):
            ChargedFormula(formula=()).to_mz_paf()

    def test_from_mz_paf_positive(self) -> None:
        # Real mzPAF chemical formulas are never "Formula:"-prefixed.
        result = ChargedFormula.from_mz_paf("+C2H2")
        assert result.formula_dict() == {"C": 2, "H": 2}

    def test_from_mz_paf_positive_rejects_negative_occurance(self) -> None:
        with pytest.raises(ValueError, match="negative occurance in positive part"):
            ChargedFormula.from_mz_paf("+C2H-2")

    def test_from_mz_paf_negative(self) -> None:
        result = ChargedFormula.from_mz_paf("-H2O")
        assert result.formula_dict() == {"H": -2, "O": -1}

    def test_from_mz_paf_negative_rejects_negative_occurance(self) -> None:
        with pytest.raises(ValueError, match="negative occurance in negative part"):
            ChargedFormula.from_mz_paf("-C2H-2")

    def test_to_mz_paf_round_trips_through_from_mz_paf(self) -> None:
        for cf in (
            ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=2),)),
            ChargedFormula(
                formula=(
                    FormulaElement(element=Element.H, occurance=2),
                    FormulaElement(element=Element.O, occurance=1),
                )
            ),
            ChargedFormula(
                formula=(
                    FormulaElement(element=Element.H, occurance=-2),
                    FormulaElement(element=Element.O, occurance=-1),
                )
            ),
        ):
            assert ChargedFormula.from_mz_paf(cf.to_mz_paf()).get_composition() == cf.get_composition()

    def test_from_mz_paf_requires_leading_sign(self) -> None:
        with pytest.raises(ValueError):
            ChargedFormula.from_mz_paf("C2H2")

    def test_from_composition(self) -> None:
        cf = ChargedFormula.from_composition({"C": 2, "H": 4}, charge=1)
        assert cf.formula_dict() == {"C": 2, "H": 4}
        assert cf.charge == 1

    def test_serialize_and_str(self) -> None:
        cf = ChargedFormula.from_string("H2O", require_formula_prefix=False)
        assert cf.serialize() == "Formula:H2O"
        assert str(cf) == "Formula:H2O"

    def test_add_combines_composition_and_charge(self) -> None:
        cf1 = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=2),), charge=1)
        cf2 = ChargedFormula(formula=(FormulaElement(element=Element.H, occurance=4),), charge=1)
        combined = cf1 + cf2
        assert combined.formula_dict() == {"C": 2, "H": 4}
        assert combined.charge == 2

    def test_add_charge_is_none_if_either_operand_has_none_charge(self) -> None:
        cf1 = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=2),), charge=1)
        cf2 = ChargedFormula(formula=(FormulaElement(element=Element.H, occurance=4),), charge=None)
        assert (cf1 + cf2).charge is None

    def test_sub_combines_composition_and_charge(self) -> None:
        cf1 = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=2),), charge=2)
        cf2 = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1),), charge=1)
        combined = cf1 - cf2
        assert combined.formula_dict() == {"C": 1}
        assert combined.charge == 1

    def test_is_neutral_true_for_none_charge(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1),), charge=None)
        assert cf.is_neutral is True

    def test_is_neutral_true_for_zero_charge(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1),), charge=0)
        assert cf.is_neutral is True

    def test_is_neutral_false_for_nonzero_charge(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1),), charge=1)
        assert cf.is_neutral is False

    def test_is_charged_true_for_nonzero_charge(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1),), charge=1)
        assert cf.is_charged is True

    def test_is_charged_false_for_none_charge(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1),), charge=None)
        assert cf.is_charged is False


class TestTagAccession:
    def test_is_valid_false_for_unknown_accession(self) -> None:
        tag = TagAccession(accession="UNKNOWN123", cv=CV.UNIMOD)
        assert tag.is_valid is False
        error = tag.validate()
        assert error is not None
        assert "UNKNOWN123" in error

    def test_get_mass_raises_for_unknown_accession(self) -> None:
        tag = TagAccession(accession="UNKNOWN123", cv=CV.UNIMOD)
        with pytest.raises(ValueError, match="Unknown modification accession"):
            tag.get_mass()

    def test_get_composition_raises_for_unknown_accession(self) -> None:
        tag = TagAccession(accession="UNKNOWN123", cv=CV.UNIMOD)
        with pytest.raises(ValueError, match="Unknown modification accession"):
            tag.get_composition()

    def test_get_charge_is_always_none(self) -> None:
        assert TagAccession(accession="35", cv=CV.UNIMOD).get_charge() is None

    @pytest.mark.parametrize(
        ("cv", "accession"),
        [
            (CV.PSI_MOD, "00046"),
            (CV.RESID, "AA0037"),
            (CV.GNOME, "G00008BG"),
            (CV.XL_MOD, "02001"),
        ],
    )
    def test_lookup_by_accession_for_each_cv(self, cv: CV, accession: str) -> None:
        tag = TagAccession(accession=accession, cv=cv)
        assert tag.is_valid is True
        assert tag.get_mass() > 0
        assert len(tag.get_composition()) > 0

    def test_unsupported_cv_raises_on_validate(self) -> None:
        tag = TagAccession(accession="X", cv=CV.CUSTOM)
        error = tag.validate()
        assert error is not None
        assert "not implemented" in error

    def test_get_composition_raises_when_cv_entry_has_no_composition(self) -> None:
        tag = TagAccession(accession=PSIMOD_NO_COMPOSITION_ID, cv=CV.PSI_MOD)
        # The mass is defined even though the composition is not.
        assert tag.get_mass() > 0
        with pytest.raises(ValueError, match="no elemental composition"):
            tag.get_composition()

    def test_from_string_serialize_round_trip(self) -> None:
        tag = TagAccession.from_string("UNIMOD:35")
        assert tag.serialize() == "UNIMOD:35"
        assert str(tag) == "UNIMOD:35"


class TestTagMass:
    def test_post_init_adds_plus_sign(self) -> None:
        assert TagMass(mass_str="10.0").mass_str == "+10.0"

    def test_post_init_preserves_minus_sign(self) -> None:
        assert TagMass(mass_str="-10.0").mass_str == "-10.0"

    def test_post_init_accepts_float_input(self) -> None:
        assert TagMass(mass_str=10.0).mass_str == "+10.0"

    def test_is_valid_false_for_unparseable_mass(self) -> None:
        tag = TagMass(mass_str="abc")
        assert tag.is_valid is False
        assert tag.validate() is not None

    def test_is_valid_true_for_parseable_mass(self) -> None:
        tag = TagMass(mass_str="15.9949")
        assert tag.is_valid is True
        assert tag.validate() is None

    def test_from_string_serialize_round_trip(self) -> None:
        tag = TagMass.from_string("+15.9949")
        assert tag.serialize() == "+15.9949"
        assert str(tag) == "+15.9949"

    def test_get_charge_is_always_none(self) -> None:
        assert TagMass(mass_str="10.0").get_charge() is None

    def test_get_composition_unimod(self) -> None:
        tag = TagMass(mass_str="79.9663", cv=CV.UNIMOD)
        assert len(tag.get_composition()) > 0

    @pytest.mark.parametrize(
        ("cv", "mass_str"),
        [
            (CV.GNOME, "1118.37484932148"),
        ],
    )
    def test_get_composition_for_cv_with_single_match(self, cv: CV, mass_str: str) -> None:
        tag = TagMass(mass_str=mass_str, cv=cv)
        assert len(tag.get_composition()) > 0

    @pytest.mark.parametrize(
        ("cv", "mass_str"),
        [
            (CV.PSI_MOD, "79.966331"),
            (CV.RESID, "79.966331"),
            (CV.XL_MOD, "138.06807961"),
        ],
    )
    def test_get_composition_raises_for_ambiguous_mass(self, cv: CV, mass_str: str) -> None:
        tag = TagMass(mass_str=mass_str, cv=cv)
        with pytest.raises(ValueError, match="Multiple modifications found"):
            tag.get_composition()

    def test_get_composition_unsupported_cv_raises(self) -> None:
        tag = TagMass(mass_str="1.0", cv=CV.CUSTOM)
        with pytest.raises(ValueError, match="not implemented"):
            tag.get_composition()


class TestPositionScore:
    def test_get_mass_is_zero(self) -> None:
        assert PositionScore(position_id="1").get_mass() == 0.0

    def test_get_composition_is_empty(self) -> None:
        assert PositionScore(position_id="1").get_composition() == Counter()

    def test_validate_is_always_none(self) -> None:
        assert PositionScore(position_id="1").validate() is None


class TestTagName:
    def test_is_valid_false_for_unknown_name(self) -> None:
        tag = TagName(name="TotallyBogusModName12345")
        assert tag.is_valid is False
        error = tag.validate()
        assert error is not None
        assert "TotallyBogusModName12345" in error

    def test_get_mass_raises_for_unknown_name(self) -> None:
        tag = TagName(name="TotallyBogusModName12345")
        with pytest.raises(ValueError, match="Unknown modification name"):
            tag.get_mass()

    def test_get_charge_is_always_none(self) -> None:
        assert TagName(name="Oxidation", cv=CV.UNIMOD).get_charge() is None

    def test_lookup_with_explicit_unimod_cv(self) -> None:
        tag = TagName(name="Oxidation", cv=CV.UNIMOD)
        assert tag.get_mass() > 0

    def test_default_cv_falls_back_to_unimod_first(self) -> None:
        # "Oxidation" is defined in UNIMOD; default cv=None should resolve it there.
        tag = TagName(name="Oxidation")
        assert tag.get_mass() == TagName(name="Oxidation", cv=CV.UNIMOD).get_mass()

    def test_default_cv_falls_back_to_psimod_when_not_in_unimod(self) -> None:
        # "O-phospho-L-serine" is only defined in PSI-MOD, not UNIMOD.
        tag = TagName(name="O-phospho-L-serine")
        assert tag.get_mass() == TagName(name="O-phospho-L-serine", cv=CV.PSI_MOD).get_mass()

    def test_get_composition_with_explicit_unimod_cv(self) -> None:
        tag = TagName(name="Oxidation", cv=CV.UNIMOD)
        assert len(tag.get_composition()) > 0

    def test_from_string_serialize_round_trip(self) -> None:
        tag = TagName.from_string("Oxidation")
        assert tag.serialize() == "Oxidation"
        assert str(tag) == "Oxidation"

    @pytest.mark.parametrize(
        ("cv", "name"),
        [
            (CV.PSI_MOD, "O-phospho-L-serine"),
            (CV.RESID, "O-phospho-L-serine"),
            (CV.GNOME, "G00008BG"),
            (CV.XL_MOD, "DSS"),
        ],
    )
    def test_lookup_by_name_for_each_cv(self, cv: CV, name: str) -> None:
        tag = TagName(name=name, cv=cv)
        assert tag.is_valid is True
        assert tag.get_mass() > 0

    def test_unsupported_cv_raises_on_validate(self) -> None:
        tag = TagName(name="Oxidation", cv=CV.CUSTOM)
        error = tag.validate()
        assert error is not None
        assert "not implemented" in error

    def test_get_composition_raises_when_cv_entry_has_no_composition(self) -> None:
        tag = TagName(name=PSIMOD_NO_COMPOSITION_NAME, cv=CV.PSI_MOD)
        assert tag.get_mass() > 0
        with pytest.raises(ValueError, match="no elemental composition"):
            tag.get_composition()

    def test_get_composition_raises_for_unknown_name(self) -> None:
        tag = TagName(name="TotallyBogusModName12345")
        with pytest.raises(ValueError, match="Unknown modification name"):
            tag.get_composition()


class TestTagInfo:
    def test_is_valid_always_true(self) -> None:
        assert TagInfo(info="some info").is_valid is True

    def test_get_mass_is_zero(self) -> None:
        assert TagInfo(info="some info").get_mass() == 0.0

    def test_get_composition_is_empty(self) -> None:
        assert TagInfo(info="some info").get_composition() == Counter()

    def test_get_charge_is_none(self) -> None:
        assert TagInfo(info="some info").get_charge() is None

    def test_from_string_serialize_round_trip(self) -> None:
        tag = TagInfo.from_string("INFO:hello")
        assert tag.serialize() == "INFO:hello"
        assert str(tag) == "INFO:hello"


class TestTagCustom:
    def test_is_valid_always_true(self) -> None:
        assert TagCustom(name="MyCustomMod").is_valid is True

    def test_get_mass_is_zero(self) -> None:
        assert TagCustom(name="MyCustomMod").get_mass() == 0.0

    def test_get_composition_is_empty(self) -> None:
        assert TagCustom(name="MyCustomMod").get_composition() == Counter()

    def test_get_charge_is_none(self) -> None:
        assert TagCustom(name="MyCustomMod").get_charge() is None

    def test_from_string_serialize_round_trip(self) -> None:
        tag = TagCustom.from_string("C:MyCustomThing")
        assert tag.serialize() == "C:MyCustomThing"
        assert str(tag) == "C:MyCustomThing"


class TestGlycanComponent:
    def test_is_valid_true_for_known_monosaccharide(self) -> None:
        gc = GlycanComponent(monosaccharide=Monosaccharide.Hex, occurance=2)
        assert gc.is_valid is True
        assert gc.get_mass() > 0
        assert len(gc.get_composition()) > 0

    def test_is_valid_false_for_unknown_monosaccharide_name(self) -> None:
        gc = GlycanComponent(monosaccharide="NotAMonosaccharide", occurance=1)  # type: ignore[arg-type]
        assert gc.is_valid is False
        assert gc.validate() is not None

    def test_get_charge_is_none(self) -> None:
        gc = GlycanComponent(monosaccharide=Monosaccharide.Hex, occurance=1)
        assert gc.get_charge() is None

    def test_charged_formula_monosaccharide_is_rejected(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.C, occurance=1),))
        with pytest.raises(NotImplementedError):
            GlycanComponent(monosaccharide=cf, occurance=1)

    def test_from_string_serialize_round_trip(self) -> None:
        gc = GlycanComponent.from_string("Hex2")
        assert gc.serialize() == "Hex2"
        assert str(gc) == "Hex2"


class TestGlycanTag:
    def test_is_valid_true(self) -> None:
        tag = GlycanTag.from_string("Glycan:Hex2")
        assert tag.is_valid is True

    def test_is_valid_false_propagates_component_error(self) -> None:
        bad_component = GlycanComponent(monosaccharide="NotReal", occurance=1)  # type: ignore[arg-type]
        tag = GlycanTag(components=(bad_component,))
        assert tag.is_valid is False
        assert tag.validate() is not None

    def test_get_charge_is_none(self) -> None:
        tag = GlycanTag.from_string("Glycan:Hex2")
        assert tag.get_charge() is None

    def test_get_mass_and_composition(self) -> None:
        tag = GlycanTag.from_string("Glycan:Hex2")
        assert tag.get_mass() > 0
        assert len(tag.get_composition()) > 0

    def test_len_and_getitem(self) -> None:
        tag = GlycanTag.from_string("Glycan:Hex2HexNAc1")
        assert len(tag) == 2
        assert tag[0].monosaccharide == Monosaccharide.Hex
        assert tag[1].monosaccharide == Monosaccharide.HexNAc

    def test_from_string_serialize_round_trip(self) -> None:
        tag = GlycanTag.from_string("Glycan:Hex2")
        assert tag.serialize() == "Glycan:Hex2"
        assert str(tag) == "Glycan:Hex2"


class TestPlacementTagMixin:
    """Exercised through LimitTag, a concrete PlacementTagMixin subclass."""

    def test_get_mass_raises(self) -> None:
        with pytest.raises(ValueError, match="has no mass"):
            LimitTag(limit=1).get_mass()

    def test_get_composition_is_empty(self) -> None:
        assert LimitTag(limit=1).get_composition() == Counter()

    def test_get_charge_is_none(self) -> None:
        assert LimitTag(limit=1).get_charge() is None

    def test_validate_is_none(self) -> None:
        assert LimitTag(limit=1).validate() is None


class TestPositionTag:
    def test_from_string_round_trip(self) -> None:
        tag = PositionTag.from_string("Position:N-term")
        assert tag.serialize() == "Position:N-term"
        assert str(tag) == "Position:N-term"

    def test_from_string_requires_prefix(self) -> None:
        with pytest.raises(ValueError, match="must start with 'Position:'"):
            PositionTag.from_string("Bad:N-term")


class TestLimitTag:
    def test_from_string_round_trip(self) -> None:
        tag = LimitTag.from_string("Limit:2")
        assert tag.serialize() == "Limit:2"
        assert str(tag) == "Limit:2"

    def test_from_string_requires_prefix(self) -> None:
        with pytest.raises(ValueError, match="must start with 'Limit:'"):
            LimitTag.from_string("Bad:2")


class TestComkpTag:
    def test_from_string_round_trip(self) -> None:
        tag = ComkpTag.from_string("Comkp")
        assert tag.serialize() == "CoMKP"
        assert str(tag) == "CoMKP"

    def test_from_string_rejects_other_strings(self) -> None:
        with pytest.raises(ValueError, match="must be 'Comkp'"):
            ComkpTag.from_string("Bad")


class TestComupTag:
    def test_from_string_round_trip(self) -> None:
        tag = ComupTag.from_string("Comup")
        assert tag.serialize() == "CoMUP"
        assert str(tag) == "CoMUP"

    def test_from_string_rejects_other_strings(self) -> None:
        with pytest.raises(ValueError, match="must be 'Comup'"):
            ComupTag.from_string("Bad")


class TestIsotopeReplacement:
    def test_is_valid_true_for_known_isotope(self) -> None:
        ir = IsotopeReplacement(element=Element.C, isotope=13)
        assert ir.is_valid is True
        original, replaced = ir.get_isotope_replacements()
        assert original.mass_number is None
        assert replaced.mass_number == 13

    def test_is_valid_false_for_unknown_isotope(self) -> None:
        ir = IsotopeReplacement(element=Element.C, isotope=999)
        assert ir.is_valid is False
        assert ir.validate() is not None

    def test_get_mass_is_zero(self) -> None:
        assert IsotopeReplacement(element=Element.C, isotope=13).get_mass() == 0.0

    def test_get_composition_is_empty(self) -> None:
        assert IsotopeReplacement(element=Element.C, isotope=13).get_composition() == Counter()

    def test_get_charge_is_none(self) -> None:
        assert IsotopeReplacement(element=Element.C, isotope=13).get_charge() is None

    def test_from_string_serialize_round_trip(self) -> None:
        ir = IsotopeReplacement.from_string("13C")
        assert ir.serialize() == "13C"
        assert str(ir) == "13C"


class TestGlobalChargeCarrier:
    def test_charged_proton_factory(self) -> None:
        carrier = GlobalChargeCarrier.charged_proton(2)
        assert carrier.get_charge() == 2
        assert carrier.is_protonated is True

    def test_to_mz_paf(self) -> None:
        carrier = GlobalChargeCarrier.charged_proton(1)
        assert carrier.to_mz_paf() == "M+H"

    def test_to_mz_paf_includes_occurance_count(self) -> None:
        # mzPAF section 4.7's own example: "[M+2Na] denotes an adduct ion with
        # two sodium atoms."
        cf = ChargedFormula(formula=(FormulaElement(element=Element.Na, occurance=1),), charge=1)
        carrier = GlobalChargeCarrier(charged_formula=cf, occurance=2)
        assert carrier.to_mz_paf() == "M+2Na"

    def test_is_valid_true(self) -> None:
        carrier = GlobalChargeCarrier.from_string("Na:z+1")
        assert carrier.is_valid is True
        assert carrier.validate() is None

    def test_get_charge_raises_when_formula_has_no_charge(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.Na, occurance=1),), charge=None)
        carrier = GlobalChargeCarrier(charged_formula=cf, occurance=1)
        with pytest.raises(ValueError, match="no defined charge"):
            carrier.get_charge()

    def test_is_protonated_false_for_non_hydrogen(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.Na, occurance=1),), charge=1)
        carrier = GlobalChargeCarrier(charged_formula=cf, occurance=1)
        assert carrier.is_protonated is False

    def test_get_composition_scales_by_occurance(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.Na, occurance=1),), charge=1)
        carrier = GlobalChargeCarrier(charged_formula=cf, occurance=3)
        composition = carrier.get_composition()
        assert list(composition.values()) == [3]

    def test_get_mass_scales_by_occurance(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.Na, occurance=1),), charge=1)
        carrier = GlobalChargeCarrier(charged_formula=cf, occurance=2)
        assert carrier.get_mass() == pytest.approx(cf.get_mass() * 2)

    def test_serialize_and_str(self) -> None:
        carrier = GlobalChargeCarrier.from_string("Na:z+1")
        assert carrier.serialize() == "Na:z+1"
        assert str(carrier) == "Na:z+1"


class TestModificationTags:
    def test_is_valid_false_when_empty(self) -> None:
        tags = ModificationTags(tags=())
        assert tags.is_valid is False
        assert tags.validate() == "ModificationTags cannot be empty"

    def test_has_placement_tags_true(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD), LimitTag(limit=1)))
        assert tags.has_placement_tags is True

    def test_has_placement_tags_false(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        assert tags.has_placement_tags is False

    def test_placement_tags_collects_all_placement_types(self) -> None:
        position = PositionTag((PositionRule(terminal=Terminal.N_TERM),))
        limit = LimitTag(limit=2)
        comkp = ComkpTag()
        comup = ComupTag()
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD), position, limit, comkp, comup))
        assert tags.placement_tags == (position, limit, comkp, comup)

    def test_placement_tags_empty_when_none_present(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        assert tags.placement_tags == ()

    def test_position_tag_property(self) -> None:
        position = PositionTag((PositionRule(terminal=Terminal.N_TERM),))
        tags = ModificationTags(tags=(position,))
        assert tags.position_tag is position

    def test_position_tag_property_none_when_absent(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        assert tags.position_tag is None

    def test_limit_tag_property(self) -> None:
        limit = LimitTag(limit=2)
        tags = ModificationTags(tags=(limit,))
        assert tags.limit_tag is limit

    def test_limit_tag_property_none_when_absent(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        assert tags.limit_tag is None

    def test_comkp_tag_property(self) -> None:
        comkp = ComkpTag()
        tags = ModificationTags(tags=(comkp,))
        assert tags.comkp_tag is comkp

    def test_comkp_tag_property_none_when_absent(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        assert tags.comkp_tag is None

    def test_comup_tag_property(self) -> None:
        comup = ComupTag()
        tags = ModificationTags(tags=(comup,))
        assert tags.comup_tag is comup

    def test_comup_tag_property_none_when_absent(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        assert tags.comup_tag is None

    def test_validate_propagates_error_from_invalid_tag(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="UNKNOWN123", cv=CV.UNIMOD),))
        error = tags.validate()
        assert error is not None
        assert "UNKNOWN123" in error

    def test_get_charge_from_charged_formula_first_tag(self) -> None:
        cf = ChargedFormula(formula=(FormulaElement(element=Element.Na, occurance=1),), charge=1)
        tags = ModificationTags(tags=(cf,))
        assert tags.get_charge() == 1

    def test_get_charge_is_none_for_non_charged_formula_first_tag(self) -> None:
        tags = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        assert tags.get_charge() is None

    def test_getitem(self) -> None:
        ta = TagAccession(accession="35", cv=CV.UNIMOD)
        tags = ModificationTags(tags=(ta,))
        assert tags[0] is ta

    def test_from_string_serialize_round_trip(self) -> None:
        tags = ModificationTags.from_string("Oxidation")
        assert tags.serialize() == "Oxidation"
        assert str(tags) == "Oxidation"


class TestModificationAmbiguousPrimary:
    def test_raises_when_score_out_of_range(self) -> None:
        with pytest.raises(ValueError, match="Score must be between 0 and 1"):
            ModificationAmbiguousPrimary(
                label="1",
                tags=ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),)),
                score=1.5,
            )

    def test_raises_when_limit_not_positive(self) -> None:
        with pytest.raises(ValueError, match="Limit must be positive"):
            ModificationAmbiguousPrimary(
                label="1",
                tags=ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),)),
                limit=0,
            )

    def test_raises_when_tags_empty(self) -> None:
        with pytest.raises(ValueError, match="tags cannot be empty"):
            ModificationAmbiguousPrimary(label="1", tags=ModificationTags(tags=()))

    def test_len_matches_underlying_tags(self) -> None:
        mod = ModificationAmbiguousPrimary(
            label="1",
            tags=ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),)),
        )
        assert len(mod) == 1

    def test_from_string_serialize_round_trip(self) -> None:
        mod = ModificationAmbiguousPrimary.from_string("1#Oxidation")
        assert mod.serialize() == "1#Oxidation"
        assert str(mod) == "1#Oxidation"


class TestModificationAmbiguousSecondary:
    def test_raises_when_score_out_of_range(self) -> None:
        with pytest.raises(ValueError, match="Score must be between 0 and 1"):
            ModificationAmbiguousSecondary(label="1", score=2.0)

    def test_accepts_valid_score(self) -> None:
        mod = ModificationAmbiguousSecondary(label="1", score=0.5)
        assert mod.score == 0.5

    def test_from_string_serialize_round_trip(self) -> None:
        mod = ModificationAmbiguousSecondary.from_string("#1")
        assert mod.serialize() == "#1"
        assert str(mod) == "#1"


class TestModificationCrossLinker:
    def test_len_is_zero_when_tags_is_none(self) -> None:
        assert len(ModificationCrossLinker(label="xl1", tags=None)) == 0

    def test_len_matches_tags_when_present(self) -> None:
        linker = ModificationCrossLinker(
            label="xl1",
            tags=ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),)),
        )
        assert len(linker) == 1

    def test_from_string_serialize_round_trip(self) -> None:
        linker = ModificationCrossLinker.from_string("XLMOD:02001#XL1")
        assert linker.serialize() == "XLMOD:02001#XL1"
        assert str(linker) == "XLMOD:02001#XL1"


class TestFixedModification:
    def test_is_valid_true(self) -> None:
        fm = FixedModification(modifications=ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),)))
        assert fm.is_valid is True

    def test_get_mass_delegates_to_modifications(self) -> None:
        mods = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        fm = FixedModification(modifications=mods)
        assert fm.get_mass() == mods.get_mass()

    def test_get_composition_delegates_to_modifications(self) -> None:
        mods = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        fm = FixedModification(modifications=mods)
        assert fm.get_composition() == mods.get_composition()

    def test_get_charge_delegates_to_modifications(self) -> None:
        mods = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        fm = FixedModification(modifications=mods)
        assert fm.get_charge() == mods.get_charge()

    def test_len_matches_underlying_modifications(self) -> None:
        mods = ModificationTags(tags=(TagAccession(accession="35", cv=CV.UNIMOD),))
        fm = FixedModification(modifications=mods)
        assert len(fm) == 1

    def test_find_indexes_n_term(self) -> None:
        fm = FixedModification.from_string("[Oxidation]@N-term")
        assert fm.find_indexes("MPEPTIDEM") == [-1]

    def test_find_indexes_c_term(self) -> None:
        fm = FixedModification.from_string("[Oxidation]@C-term")
        assert fm.find_indexes("MPEPTIDEM") == [-2]

    def test_find_indexes_anywhere_matches_residue(self) -> None:
        fm = FixedModification.from_string("[Oxidation]@M")
        assert fm.find_indexes("MPEPTIDEM") == [0, 8]

    def test_find_indexes_combines_all_position_rules(self) -> None:
        fm = FixedModification.from_string("[Oxidation]@M,C-term,N-term")
        assert fm.find_indexes("MPEPTIDEM") == [0, 8, -2, -1]

    def test_from_string_serialize_round_trip(self) -> None:
        fm = FixedModification.from_string("[Oxidation]@M")
        assert fm.serialize() == "[Oxidation]@M"
        assert str(fm) == "[Oxidation]@M"


class TestSequenceElement:
    def test_get_mass_adds_amino_acid_and_modification_mass(self) -> None:
        se = SequenceElement(amino_acid=AminoAcid.A)
        aa = AA_LOOKUP.one_letter(AminoAcid.A)
        assert se.get_mass() == aa.monoisotopic_mass

    def test_get_mass_raises_for_amino_acid_with_no_defined_mass(self) -> None:
        # B (Asx) is an ambiguous placeholder amino acid with no defined mass in tacular.
        se = SequenceElement(amino_acid=AminoAcid.B)
        with pytest.raises(ValueError, match="Unknown mass for amino acid"):
            se.get_mass()

    def test_get_composition_returns_amino_acid_composition(self) -> None:
        se = SequenceElement(amino_acid=AminoAcid.A)
        aa = AA_LOOKUP.one_letter(AminoAcid.A)
        assert se.get_composition() == Counter(aa.composition)

    def test_get_composition_raises_for_amino_acid_with_no_defined_composition(self) -> None:
        se = SequenceElement(amino_acid=AminoAcid.B)
        with pytest.raises(ValueError, match="Unknown composition for amino acid"):
            se.get_composition()

    def test_get_composition_merges_modification_composition(self) -> None:
        se = SequenceElement.from_string("M[Oxidation]")
        unmodified = SequenceElement(amino_acid=AminoAcid.M)
        assert se.get_composition() != unmodified.get_composition()
        assert se.get_mass() > unmodified.get_mass()

    def test_from_string_serialize_round_trip(self) -> None:
        se = SequenceElement.from_string("M[Oxidation]")
        assert se.serialize() == "M[Oxidation]"
        assert str(se) == "M[Oxidation]"


class TestSequenceRegion:
    def test_get_mass_sums_sequence_and_modifications(self) -> None:
        region = SequenceRegion(
            sequence=(SequenceElement(amino_acid=AminoAcid.P), SequenceElement(amino_acid=AminoAcid.E)),
            modifications=(),
            ambiguous=False,
        )
        expected = sum(se.get_mass() for se in region.sequence)
        assert region.get_mass() == pytest.approx(expected)

    def test_get_composition_merges_sequence_and_modifications(self) -> None:
        region = SequenceRegion(
            sequence=(SequenceElement(amino_acid=AminoAcid.P), SequenceElement(amino_acid=AminoAcid.E)),
            modifications=(),
            ambiguous=False,
        )
        expected = region.sequence[0].get_composition() + region.sequence[1].get_composition()
        assert region.get_composition() == expected

    def test_from_string_serialize_round_trip(self) -> None:
        region = SequenceRegion.from_string("(PEPM[Oxidation]TIDE)")
        assert region.serialize() == "(PEPM[Oxidation]TIDE)"
        assert str(region) == "(PEPM[Oxidation]TIDE)"


class TestPeptidoform:
    """Peptidoform.from_string / serialize round-tripping.

    Note: parsing (`Peptidoform.from_string`) is not yet implemented upstream in
    `parse_peptidoform` -- it always raises NotImplementedError. Peptidoforms are
    otherwise constructed directly (e.g. via the annotation layer), so tests build
    them directly and only assert the documented NotImplementedError for parsing.
    """

    def _make(self) -> Peptidoform:
        return Peptidoform(sequence=(SequenceElement(amino_acid=AminoAcid.P), SequenceElement(amino_acid=AminoAcid.E)))

    def test_get_mass(self) -> None:
        pf = self._make()
        expected = sum(se.get_mass() for se in pf.sequence)
        assert pf.get_mass() == pytest.approx(expected)

    def test_get_composition(self) -> None:
        pf = self._make()
        expected = pf.sequence[0].get_composition() + pf.sequence[1].get_composition()
        assert pf.get_composition() == expected

    def test_from_string_not_implemented(self) -> None:
        with pytest.raises(NotImplementedError):
            Peptidoform.from_string("PEPTIDE")

    def test_serialize_and_str(self) -> None:
        pf = self._make()
        assert pf.serialize() == "PE"
        assert str(pf) == "PE"


class TestPeptidoformIon:
    def _make(self) -> PeptidoformIon:
        pf = Peptidoform(sequence=(SequenceElement(amino_acid=AminoAcid.P), SequenceElement(amino_acid=AminoAcid.E)))
        return PeptidoformIon(peptidoforms=(pf,))

    def test_get_mass_not_implemented(self) -> None:
        with pytest.raises(NotImplementedError):
            self._make().get_mass()

    def test_get_composition_not_implemented(self) -> None:
        with pytest.raises(NotImplementedError):
            self._make().get_composition()

    def test_from_string_not_implemented(self) -> None:
        with pytest.raises(NotImplementedError):
            PeptidoformIon.from_string("PEPTIDE/2")

    def test_serialize_and_str(self) -> None:
        pfi = self._make()
        assert pfi.serialize() == "PE"
        assert str(pfi) == "PE"


class TestCompoundPeptidoformIon:
    def _make(self) -> CompoundPeptidoformIon:
        pf = Peptidoform(sequence=(SequenceElement(amino_acid=AminoAcid.P), SequenceElement(amino_acid=AminoAcid.E)))
        pfi = PeptidoformIon(peptidoforms=(pf,))
        return CompoundPeptidoformIon(peptidoform_ions=(pfi,))

    def test_get_mass_not_implemented(self) -> None:
        with pytest.raises(NotImplementedError):
            self._make().get_mass()

    def test_get_composition_not_implemented(self) -> None:
        with pytest.raises(NotImplementedError):
            self._make().get_composition()

    def test_from_string_not_implemented(self) -> None:
        with pytest.raises(NotImplementedError):
            CompoundPeptidoformIon.from_string("PEPTIDE")

    def test_serialize_and_str(self) -> None:
        cpi = self._make()
        assert cpi.serialize() == "PE"
        assert str(cpi) == "PE"
