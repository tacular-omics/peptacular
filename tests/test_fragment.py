import unittest

import peptacular as pt


class TestFragment(unittest.TestCase):
    def test_calculate_mass1(self):
        sequence1 = "PE[+10]PTIDE[+10]"
        sequence2 = "<[+10]@E>PEPTIDE"

        # ensure the mass calculations are the same
        mass1 = pt.mass(sequence1, charge=0, ion_type="y", monoisotopic=True)
        mass2 = pt.mass(sequence2, charge=0, ion_type="y", monoisotopic=True)
        self.assertEqual(mass1, mass2)

    def test_calculate_mz_with_unmodified_peptide(self):
        sequence = "PEPTIDE"
        places = 2

        self.assertAlmostEqual(
            799.359964,
            pt.mz(sequence, charge=0, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            800.367241,
            pt.mz(sequence, charge=1, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            400.687259,
            pt.mz(sequence, charge=2, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            267.460598,
            pt.mz(sequence, charge=3, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            200.847268,
            pt.mz(sequence, charge=4, ion_type="y", monoisotopic=True),
            places,
        )

        # Average mass is off by 0.004 Da
        self.assertAlmostEqual(
            799.822520,
            pt.mz(sequence, charge=0, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            800.829796,
            pt.mz(sequence, charge=1, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            400.918536,
            pt.mz(sequence, charge=2, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            267.614783,
            pt.mz(sequence, charge=3, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            200.962906,
            pt.mz(sequence, charge=4, ion_type="y", monoisotopic=False),
            places,
        )

    def test_calculate_mz_with_modified_peptide(self):
        sequence = "[+15]-P[-10]EPTIDE[+100]"
        places = 2

        self.assertAlmostEqual(
            799.359964 + 105,
            pt.mz(sequence, charge=0, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            800.367241 + 105 / 1,
            pt.mz(sequence, charge=1, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            400.687259 + 105 / 2,
            pt.mz(sequence, charge=2, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            267.460598 + 105 / 3,
            pt.mz(sequence, charge=3, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            200.847268 + 105 / 4,
            pt.mz(sequence, charge=4, ion_type="y", monoisotopic=True),
            places,
        )

        self.assertAlmostEqual(
            799.822520 + 105,
            pt.mz(sequence, charge=0, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            800.829796 + 105 / 1,
            pt.mz(sequence, charge=1, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            400.918536 + 105 / 2,
            pt.mz(sequence, charge=2, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            267.614783 + 105 / 3,
            pt.mz(sequence, charge=3, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            200.962906 + 105 / 4,
            pt.mz(sequence, charge=4, ion_type="y", monoisotopic=False),
            places,
        )

    def test_calculate_mass_with_modified_peptide(self):
        sequence = "PEPTIDE"
        places = 2

        self.assertAlmostEqual(
            799.359964,
            pt.mass(sequence, charge=0, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            799.359964 + pt.PROTON_MASS * 1,
            pt.mass(sequence, charge=1, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            799.359964 + pt.PROTON_MASS * 2,
            pt.mass(sequence, charge=2, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            799.359964 + pt.PROTON_MASS * 3,
            pt.mass(sequence, charge=3, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            799.359964 + pt.PROTON_MASS * 4,
            pt.mass(sequence, charge=4, ion_type="y", monoisotopic=True),
            places,
        )

        self.assertAlmostEqual(
            799.822520,
            pt.mass(sequence, charge=0, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            799.822520 + pt.PROTON_MASS * 1,
            pt.mass(sequence, charge=1, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            799.822520 + pt.PROTON_MASS * 2,
            pt.mass(sequence, charge=2, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            799.822520 + pt.PROTON_MASS * 3,
            pt.mass(sequence, charge=3, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            799.822520 + pt.PROTON_MASS * 4,
            pt.mass(sequence, charge=4, ion_type="y", monoisotopic=False),
            places,
        )

    def test_calculate_mass_with_unmodified_peptide(self):
        sequence = "[+15]-P[-10]EPTIDE[+100]"
        places = 2

        self.assertAlmostEqual(
            799.359964 + 105,
            pt.mass(sequence, charge=0, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            799.359964 + pt.PROTON_MASS * 1 + 105,
            pt.mass(sequence, charge=1, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            799.359964 + pt.PROTON_MASS * 2 + 105,
            pt.mass(sequence, charge=2, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            799.359964 + pt.PROTON_MASS * 3 + 105,
            pt.mass(sequence, charge=3, ion_type="y", monoisotopic=True),
            places,
        )
        self.assertAlmostEqual(
            799.359964 + pt.PROTON_MASS * 4 + 105,
            pt.mass(sequence, charge=4, ion_type="y", monoisotopic=True),
            places,
        )

        self.assertAlmostEqual(
            799.822520 + 105,
            pt.mass(sequence, charge=0, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            799.822520 + pt.PROTON_MASS * 1 + 105,
            pt.mass(sequence, charge=1, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            799.822520 + pt.PROTON_MASS * 2 + 105,
            pt.mass(sequence, charge=2, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            799.822520 + pt.PROTON_MASS * 3 + 105,
            pt.mass(sequence, charge=3, ion_type="y", monoisotopic=False),
            places,
        )
        self.assertAlmostEqual(
            799.822520 + pt.PROTON_MASS * 4 + 105,
            pt.mass(sequence, charge=4, ion_type="y", monoisotopic=False),
            places,
        )

    def test_all_aa_and_ions_charge1(self):
        seq = "VWPSDCYTAIMHGQENFLKR"

        pyteomics_fragments = {
            "a": [
                2349.1267039886393,
                2193.0255929650393,
                2064.9306299510395,
                1951.8465659739097,
                1804.7781520609199,
                1690.73522461978,
                1561.69263153181,
                1433.63405402653,
                1376.61259030596,
                1239.55367844751,
                1108.51319353452,
                995.4291295573901,
                924.39201577268,
                823.34433730427,
                660.2810087717199,
                557.27182398701,
                442.24488096318,
                355.21285255890996,
                258.16008871005994,
                72.08077576019998,
            ],
            "b": [
                2377.1216186081992,
                2221.0205075845993,
                2092.9255445705994,
                1979.8414805934697,
                1832.7730666804798,
                1718.73013923934,
                1589.68754615137,
                1461.62896864609,
                1404.60750492552,
                1267.54859306707,
                1136.50810815408,
                1023.4240441769501,
                952.3869303922401,
                851.33925192383,
                688.27592339128,
                585.26673860657,
                470.23979558274,
                383.20776717846996,
                286.15500332961994,
                100.07569037975999,
            ],
            "c": [
                2394.1481677092092,
                2238.0470566856093,
                2109.9520936716094,
                1996.8680296944797,
                1849.7996157814898,
                1735.75668834035,
                1606.71409525238,
                1478.6555177471,
                1421.63405402653,
                1284.57514216808,
                1153.53465725509,
                1040.45059327796,
                969.4134794932501,
                868.36580102484,
                705.30247249229,
                602.29328770758,
                487.26634468375,
                400.23431627947997,
                303.18155243062995,
                117.10223948077,
            ],
            "x": [
                2421.111447847319,
                2322.0430339343293,
                2135.963720984469,
                2038.9109571356198,
                1951.8789287313498,
                1836.8519857075198,
                1733.8428009228096,
                1570.7794723902598,
                1469.73179392185,
                1398.69468013714,
                1285.61061616001,
                1154.57013124702,
                1017.5112193885698,
                960.4897556679997,
                832.4311781627198,
                703.3885850747499,
                589.3456576336098,
                442.27724372062,
                329.19317974348996,
                201.09821672949,
            ],
            "y": [
                2395.1321832918993,
                2296.0637693789095,
                2109.9844564290493,
                2012.9316925801998,
                1925.8996641759297,
                1810.8727211520998,
                1707.8635363673895,
                1544.8002078348397,
                1443.75252936643,
                1372.71541558172,
                1259.63135160459,
                1128.5908666916,
                991.5319548331498,
                934.5104911125798,
                806.4519136072998,
                677.40932051933,
                563.3663930781898,
                416.2979791652,
                303.21391518806996,
                175.11895217407,
            ],
            "z": [
                2378.1056341908893,
                2279.0372202778995,
                2092.9579073280393,
                1995.9051434791897,
                1908.8731150749197,
                1793.8461720510898,
                1690.8369872663795,
                1527.7736587338297,
                1426.72598026542,
                1355.68886648071,
                1242.60480250358,
                1111.56431759059,
                974.5054057321398,
                917.4839420115698,
                789.4253645062898,
                660.38277141832,
                546.3398439771798,
                399.27143006419,
                286.18736608705996,
                158.09240307305998,
            ],
        }

        places = 5

        annot = pt.ProFormaAnnotation(seq)
        for i in "abcxyz":
            frags = annot.fragment([i], [1])

            for n, (pt_frag, py_frag) in enumerate(zip(frags, pyteomics_fragments[i][::-1])):
                self.assertAlmostEqual(
                    pt_frag.mz,
                    py_frag,
                    places,
                    msg=f"Failed for {n}th {i} ion type, pt_frag: {pt_frag}, py_frag: {py_frag}",
                )

    def test_all_aa_and_ions_charge2(self):
        seq = "VWPSDCYTAIMHGQENFLKR"

        pyteomics_fragments = {
            "a": [
                1175.0669902277048,
                1097.0164347159048,
                1032.9689532089048,
                976.4269212203399,
                902.8927142638449,
                845.871250543275,
                781.34995399929,
                717.32066524665,
                688.8099333863651,
                620.28047745714,
                554.760235000645,
                498.21820301208004,
                462.699646119725,
                412.17580688552,
                330.644142619245,
                279.13955022689,
                221.626078714975,
                178.11006451283998,
                129.58368258841497,
                36.544026113484996,
            ],
            "b": [
                1189.0644475374847,
                1111.0138920256848,
                1046.9664105186848,
                990.4243785301198,
                916.8901715736249,
                859.868707853055,
                795.34741130907,
                731.31812255643,
                702.807390696145,
                634.27793476692,
                568.757692310425,
                512.21566032186,
                476.69710342950503,
                426.1732641953,
                344.641599929025,
                293.13700753667,
                235.623536024755,
                192.10752182262,
                143.58113989819498,
                50.541483423265,
            ],
            "c": [
                1197.5777220879897,
                1119.5271665761898,
                1055.4796850691898,
                998.9376530806248,
                925.4034461241299,
                868.38198240356,
                803.860685859575,
                739.831397106935,
                711.32066524665,
                642.791209317425,
                577.27096686093,
                520.728934872365,
                485.21037798001004,
                434.686538745805,
                353.15487447953,
                301.650282087175,
                244.13681057526,
                200.620796373125,
                152.09441444869998,
                59.05475797377,
            ],
            "x": [
                1211.0593621570447,
                1161.5251552005498,
                1068.4854987256197,
                1019.9591168011949,
                976.4431025990599,
                918.9296310871449,
                867.4250386947898,
                785.8933744285149,
                735.3695351943101,
                699.850978301955,
                643.30894631339,
                577.788703856895,
                509.2592479276699,
                480.74851606738486,
                416.7192273147449,
                352.19793077075997,
                295.1764670501899,
                221.642260093695,
                165.10022810513,
                101.05274659812999,
            ],
            "y": [
                1198.0697298793348,
                1148.5355229228398,
                1055.4958664479097,
                1006.9694845234849,
                963.4534703213499,
                905.9399988094349,
                854.4354064170798,
                772.9037421508049,
                722.3799029166,
                686.861346024245,
                630.31931403568,
                564.799071579185,
                496.2696156499599,
                467.7588837896749,
                403.7295950370349,
                339.20829849305,
                282.1868347724799,
                208.652627815985,
                152.11059582741998,
                88.06311432041998,
            ],
            "z": [
                1189.5564553288298,
                1140.0222483723348,
                1046.9825918974047,
                998.4562099729799,
                954.9401957708449,
                897.4267242589299,
                845.9221318665748,
                764.3904676002999,
                713.866628366095,
                678.34807147374,
                621.806039485175,
                556.28579702868,
                487.7563410994549,
                459.2456092391699,
                395.2163204865299,
                330.695023942545,
                273.6735602219749,
                200.13935326548,
                143.59732127691498,
                79.54983976991498,
            ],
        }

        places = 5
        annot = pt.ProFormaAnnotation.parse(seq)

        for i in "abcxyz":
            frags = annot.fragment(ion_types=[i], charges=[2])

            for n, (pt_frag, py_frag) in enumerate(zip(frags, pyteomics_fragments[i][::-1])):
                try:
                    self.assertAlmostEqual(
                        pt_frag.mz,
                        py_frag,
                        places,
                        msg=f"Failed for {n}th {i} ion type, pt_frag: {pt_frag}, py_frag: {py_frag}",
                    )
                except AssertionError as e:
                    self.fail(f"Assertion failed: {e}")

    def test_fragment_with_netural_delta(self):
        seq = "PEPTIDE"
        annot = pt.ProFormaAnnotation.parse(seq)

        frags = annot.fragment(ion_types=["y"], charges=[1], monoisotopic=False)

        expected_masses = [
            148.136536,
            263.223936,
            376.381576,
            477.485456,
            574.600636,
            703.714616,
            800.829796,
        ]

        for frag, expected_mass in zip(frags, expected_masses, strict=True):
            self.assertAlmostEqual(
                frag.mz,
                expected_mass,
                places=2,
                msg=f"Failed for fragment {frag}, expected mass: {expected_mass}",
            )

        frags = annot.fragment(ion_types=["y"], charges=[1], monoisotopic=False, neutral_deltas=["H2O"])

        expected_masses = [
            148.136536,
            148.136536 - 18.015286432429832,
            263.223936,
            263.223936 - 18.015286432429832,
            376.381576,
            376.381576 - 18.015286432429832,
            477.485456,
            477.485456 - 18.015286432429832,
            574.600636,
            574.600636 - 18.015286432429832,
            703.714616,
            703.714616 - 18.015286432429832,
            800.829796,
            800.829796 - 18.015286432429832,
        ]

        for frag, expected_mass in zip(frags, expected_masses, strict=True):
            self.assertAlmostEqual(
                frag.mz,
                expected_mass,
                places=2,
                msg=f"Failed for fragment {frag}, expected mass: {expected_mass}",
            )

        frags = annot.fragment(
            ion_types=["y"],
            charges=[1],
            monoisotopic=False,
            neutral_deltas=["H2O"],
            max_ndeltas=2,
        )

        expected_masses = [
            # E (Only one water to lose)
            148.136536,
            148.136536 - 18.015286432429832,
            # DE
            263.223936,
            263.223936 - 18.015286432429832,
            263.223936 - 18.015286432429832 * 2,
            # IDE
            376.381576,
            376.381576 - 18.015286432429832,
            376.381576 - 18.015286432429832 * 2,
            # TIDE
            477.485456,
            477.485456 - 18.015286432429832,
            477.485456 - 18.015286432429832 * 2,
            # PTIDE
            574.600636,
            574.600636 - 18.015286432429832,
            574.600636 - 18.015286432429832 * 2,
            # EPTIDE
            703.714616,
            703.714616 - 18.015286432429832,
            703.714616 - 18.015286432429832 * 2,
            # PEPTIDE
            800.829796,
            800.829796 - 18.015286432429832,
            800.829796 - 18.015286432429832 * 2,
        ]

        for frag, expected_mass in zip(frags, expected_masses, strict=True):
            self.assertAlmostEqual(
                frag.mz,
                expected_mass,
                places=2,
                msg=f"Failed for fragment {frag}, expected mass: {expected_mass}",
            )

    def test_fragmenty(self):
        annot = pt.parse("PEPTIDE/2")

        frags = annot.fragment(ion_types=["y"])
        self.assertEqual(len(frags), 7)
        self.assertEqual(frags[0].position, 1)
        self.assertEqual(frags[0].parent_sequence, "PEPTIDE/2")
        self.assertEqual(frags[0].charge_state, 2)
        self.assertEqual(frags[0].sequence, "E/2")
        self.assertAlmostEqual(frags[0].mz, pt.mz("E/2", ion_type="y"))

    def test_fragmentb(self):
        annot = pt.parse("PEPTIDE/2")

        frags = annot.fragment(ion_types=["b"])
        self.assertEqual(len(frags), 7)
        self.assertEqual(frags[0].position, 1)
        self.assertEqual(frags[0].parent_sequence, "PEPTIDE/2")
        self.assertEqual(frags[0].charge_state, 2)
        self.assertEqual(frags[0].sequence, "P/2")
        self.assertAlmostEqual(frags[0].mz, pt.mz("P/2", ion_type="b"))

    def test_fragment_immonium(self):
        annot = pt.parse("PEPTIDE/2")

        frags = annot.fragment(ion_types=["i"])
        self.assertEqual(len(frags), 7)
        self.assertEqual(frags[0].position, 1)
        self.assertEqual(frags[0].parent_sequence, "PEPTIDE/2")
        self.assertEqual(frags[0].charge_state, 2)
        self.assertEqual(frags[0].sequence, "P/2")
        self.assertAlmostEqual(frags[0].mz, pt.mz("P/2", ion_type="i"))

        self.assertEqual(frags[1].position, 2)
        self.assertEqual(frags[1].parent_sequence, "PEPTIDE/2")
        self.assertEqual(frags[1].charge_state, 2)
        self.assertEqual(frags[1].sequence, "E/2")
        self.assertAlmostEqual(frags[1].mz, pt.mz("E/2", ion_type="i"))

    def test_fragment_internal(self):
        annot = pt.parse("PEPT/2")

        frags = annot.fragment(ion_types=["by"])
        self.assertEqual(len(frags), 3)  # EP, E, P

        self.assertEqual(frags[0].position, (2, 2))
        self.assertEqual(frags[0].parent_sequence, "PEPT/2")
        self.assertEqual(frags[0].charge_state, 2)
        self.assertEqual(frags[0].sequence, "E/2")
        self.assertAlmostEqual(frags[0].mz, pt.mz("E/2", ion_type="by"))

        self.assertEqual(frags[1].position, (2, 3))
        self.assertEqual(frags[1].parent_sequence, "PEPT/2")
        self.assertEqual(frags[1].charge_state, 2)
        self.assertEqual(frags[1].sequence, "EP/2")
        self.assertAlmostEqual(frags[1].mz, pt.mz("EP/2", ion_type="by"))

        self.assertEqual(frags[2].position, (3, 3))
        self.assertEqual(frags[2].parent_sequence, "PEPT/2")
        self.assertEqual(frags[2].charge_state, 2)
        self.assertEqual(frags[2].sequence, "P/2")
        self.assertAlmostEqual(frags[2].mz, pt.mz("P/2", ion_type="by"))


if __name__ == "__main__":
    unittest.main()
