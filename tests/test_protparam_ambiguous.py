"""An ambiguous residue is resolved when every amino acid it stands for agrees, and refused when they do not.

The aliphatic index and the instability index are mole percents, so a residue counted in the
length but not in the sum dilutes the result. These tests pin which residues each ProtParam
number can resolve and which it refuses.
"""

import pytest

import peptacular as pt
from peptacular.diagnostics import PeptacularError
from peptacular.property.core import AMBIGUOUS_AMINO_ACID_MAP, aliphatic_index, instability_index, molar_extinction_coefficient

AMBIGUOUS = ["B", "J", "Z", "X"]
OUTSIDE_THE_TWENTY = ["U", "O"]


class TestAliphaticIndexResolvesWhatItCan:
    """J is leucine or isoleucine, both 3.9, so it is counted rather than diluted."""

    def test_j_counts_as_leucine_or_isoleucine(self):
        # In AJILAVIL the J sits where AVILAVIL has a valine, so the two sequences differ:
        # 1 + 3.9 + 3.9 + 3.9 + 1 + 2.9 + 3.9 + 3.9 = 24.4 over eight residues is 305.0,
        # against 292.5 for AVILAVIL. Counting the J as zero gave 256.25.
        assert aliphatic_index("AJILAVIL") == pytest.approx(305.0)
        assert pt.aliphatic_index("AJILAVIL") == pytest.approx(305.0)
        assert aliphatic_index("AJILAVIL") == pytest.approx(aliphatic_index("ALILAVIL"))
        assert aliphatic_index("AJILAVIL") == pytest.approx(aliphatic_index("AIILAVIL"))

    def test_j_alone_is_the_ile_leu_factor(self):
        assert aliphatic_index("J") == pytest.approx(390.0)
        assert aliphatic_index("JJJJ") == pytest.approx(390.0)

    @pytest.mark.parametrize("sequence", ["ABVIL", "AZVIL"])
    def test_b_and_z_are_zero_and_stay_in_the_length(self, sequence: str):
        # Asp, Asn, Glu and Gln are not aliphatic, so zero is exact and the residue belongs
        # in the denominator: 1 + 2.9 + 3.9 + 3.9 = 11.7 over five residues.
        assert aliphatic_index(sequence) == pytest.approx(234.0)

    def test_a_residue_outside_the_twenty_is_a_determinate_zero(self):
        # U and O are defined residues and neither is Ala, Val, Ile or Leu, so their
        # contribution is a fact rather than a guess.
        assert aliphatic_index("PEPUIDE") == pytest.approx(aliphatic_index("PEPTIDE"))
        assert aliphatic_index("PEPOIDE") == pytest.approx(aliphatic_index("PEPTIDE"))

    def test_x_is_refused_rather_than_diluted(self):
        with pytest.raises(PeptacularError, match="aliphatic index of the ambiguous residue 'X'"):
            aliphatic_index("AVILAVIX")
        with pytest.raises(PeptacularError, match="aliphatic index of the ambiguous residue 'X'"):
            pt.aliphatic_index("AVILAVIX")

    def test_the_refusal_explains_itself(self):
        with pytest.raises(PeptacularError) as info:
            aliphatic_index("X")
        message = str(info.value)
        assert "any amino acid" in message
        assert "differ" in message

    @pytest.mark.parametrize("sequence", ["X", "AX", "AVILAVIX", "XAVILAVIL"])
    def test_x_is_never_counted_as_zero(self, sequence: str):
        # The defect this replaces: X stayed in the length and contributed nothing, so
        # eight aliphatic residues plus one X read 243.75 instead of refusing.
        with pytest.raises(PeptacularError):
            aliphatic_index(sequence)

    def test_annotation_accessor_refuses_x_too(self):
        annot = pt.parse("AVILAVIX")
        with pytest.raises(PeptacularError):
            _ = annot.prop.aliphatic_index


class TestExtinctionCoefficientResolvesWhatItCan:
    """B, J and Z never stand for a chromophore, so they are exact zeros. X can be Trp or Tyr."""

    @pytest.mark.parametrize("code", ["B", "J", "Z"])
    def test_a_non_chromophore_ambiguity_code_is_an_exact_zero(self, code: str):
        assert molar_extinction_coefficient(f"WY{code}CC") == 6990
        assert molar_extinction_coefficient(f"WY{code}CC", cystines=True) == 7115

    def test_x_is_refused(self):
        with pytest.raises(PeptacularError, match="molar extinction coefficient of the ambiguous residue 'X'"):
            molar_extinction_coefficient("WYX")
        with pytest.raises(PeptacularError):
            pt.molar_extinction_coefficient("WYX")

    def test_the_refusal_names_the_chromophores_it_could_be(self):
        with pytest.raises(PeptacularError) as info:
            molar_extinction_coefficient("X")
        message = str(info.value)
        assert "tryptophan" in message
        assert "tyrosine" in message

    def test_a_residue_outside_the_twenty_does_not_absorb(self):
        assert molar_extinction_coefficient("WYU") == 6990
        assert molar_extinction_coefficient("WYO") == 6990

    def test_the_benchmark_values_are_unchanged(self):
        assert molar_extinction_coefficient("WYCCWY") == 13980
        assert molar_extinction_coefficient("WYCCWY", cystines=True) == 14105
        assert molar_extinction_coefficient("PEPTIDE") == 0

    def test_annotation_accessor_refuses_x_too(self):
        annot = pt.parse("WYX")
        with pytest.raises(PeptacularError):
            annot.prop.molar_extinction_coefficient()


class TestInstabilityIndexResolvesNothing:
    """No ambiguity code can be resolved here: the DIWV values of the two expansions differ."""

    @pytest.mark.parametrize("code", AMBIGUOUS + OUTSIDE_THE_TWENTY)
    def test_a_residue_without_a_diwv_entry_is_refused(self, code: str):
        with pytest.raises(PeptacularError, match="No DIWV entry for dipeptide"):
            instability_index(f"PEP{code}IDE")

    def test_the_refusal_names_the_dipeptide(self):
        with pytest.raises(PeptacularError, match="'PX'"):
            instability_index("PEPXIDE")

    def test_j_is_refused_although_the_aliphatic_index_resolves_it(self):
        # The two numbers disagree on J on purpose: both Leu and Ile occupy the same
        # aliphatic volume, while their DIWV rows differ, so one is determined and the
        # other is not.
        assert aliphatic_index("AJILAVIL") == pytest.approx(305.0)
        with pytest.raises(PeptacularError):
            instability_index("AJILAVIL")

    def test_the_benchmark_values_are_unchanged(self):
        assert instability_index("PEPTIDE") == pytest.approx(60.9143, abs=1e-4)
        assert instability_index("ACDEFGHIKLMNPQRSTVWY") == pytest.approx(84.74, abs=1e-4)


class TestAResolvedCodeAgreesWithEveryAminoAcidItStandsFor:
    """The rule itself: a code is resolved only when all of its expansions give one answer."""

    @pytest.mark.parametrize("code", ["B", "J", "Z"])
    def test_aliphatic_index_matches_each_expansion(self, code: str):
        for option in AMBIGUOUS_AMINO_ACID_MAP[code]:
            assert aliphatic_index(f"AVIL{code}AVIL") == pytest.approx(aliphatic_index(f"AVIL{option}AVIL"))

    @pytest.mark.parametrize("code", ["B", "J", "Z"])
    def test_extinction_coefficient_matches_each_expansion(self, code: str):
        for option in AMBIGUOUS_AMINO_ACID_MAP[code]:
            assert molar_extinction_coefficient(f"WY{code}CC", cystines=True) == molar_extinction_coefficient(f"WY{option}CC", cystines=True)

    def test_x_is_refused_by_all_three(self):
        for fn in (aliphatic_index, molar_extinction_coefficient, instability_index):
            with pytest.raises(PeptacularError):
                fn("AVILXAVIL")
