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
        # /2 precursor → default charges are [1] (range 1..charge-1)
        annot = pt.parse("PEPTIDE/2")

        frags = annot.fragment(ion_types=["y"])
        self.assertEqual(len(frags), 7)
        self.assertEqual(frags[0].position, 1)
        self.assertEqual(frags[0].parent_sequence, "PEPTIDE/1")
        self.assertEqual(frags[0].charge_state, 1)
        self.assertEqual(frags[0].sequence, "E/1")
        self.assertAlmostEqual(frags[0].mz, pt.mz("E/1", ion_type="y"))

    def test_fragmentb(self):
        # /2 precursor → default charges are [1] (range 1..charge-1)
        annot = pt.parse("PEPTIDE/2")

        frags = annot.fragment(ion_types=["b"])
        self.assertEqual(len(frags), 7)
        self.assertEqual(frags[0].position, 1)
        self.assertEqual(frags[0].parent_sequence, "PEPTIDE/1")
        self.assertEqual(frags[0].charge_state, 1)
        self.assertEqual(frags[0].sequence, "P/1")
        self.assertAlmostEqual(frags[0].mz, pt.mz("P/1", ion_type="b"))

    def test_fragment_immonium(self):
        # /2 precursor → default charges are [1] (range 1..charge-1)
        annot = pt.parse("PEPTIDE/2")

        frags = annot.fragment(ion_types=["i"])
        self.assertEqual(len(frags), 7)
        self.assertEqual(frags[0].position, 1)
        self.assertEqual(frags[0].parent_sequence, "PEPTIDE/1")
        self.assertEqual(frags[0].charge_state, 1)
        self.assertEqual(frags[0].sequence, "P/1")
        self.assertAlmostEqual(frags[0].mz, pt.mz("P/1", ion_type="i"))

        self.assertEqual(frags[1].position, 2)
        self.assertEqual(frags[1].parent_sequence, "PEPTIDE/1")
        self.assertEqual(frags[1].charge_state, 1)
        self.assertEqual(frags[1].sequence, "E/1")
        self.assertAlmostEqual(frags[1].mz, pt.mz("E/1", ion_type="i"))

    def test_fragment_internal(self):
        # /2 precursor → default charges are [1] (range 1..charge-1)
        annot = pt.parse("PEPT/2")

        frags = annot.fragment(ion_types=["by"])
        self.assertEqual(len(frags), 3)  # EP, E, P

        self.assertEqual(frags[0].position, (2, 2))
        self.assertEqual(frags[0].parent_sequence, "PEPT/1")
        self.assertEqual(frags[0].charge_state, 1)
        self.assertEqual(frags[0].sequence, "E/1")
        self.assertAlmostEqual(frags[0].mz, pt.mz("E/1", ion_type="by"))

        self.assertEqual(frags[1].position, (2, 3))
        self.assertEqual(frags[1].parent_sequence, "PEPT/1")
        self.assertEqual(frags[1].charge_state, 1)
        self.assertEqual(frags[1].sequence, "EP/1")
        self.assertAlmostEqual(frags[1].mz, pt.mz("EP/1", ion_type="by"))

        self.assertEqual(frags[2].position, (3, 3))
        self.assertEqual(frags[2].parent_sequence, "PEPT/1")
        self.assertEqual(frags[2].charge_state, 1)
        self.assertEqual(frags[2].sequence, "P/1")
        self.assertAlmostEqual(frags[2].mz, pt.mz("P/1", ion_type="by"))

    def test_fragment_default_charges_range(self):
        # /3 precursor → default charges [1, 2]
        frags3 = pt.parse("PEPTIDE/3").fragment(ion_types=["b"])
        charges3 = sorted({f.charge_state for f in frags3})
        self.assertEqual(charges3, [1, 2])

        # /2 precursor → default charges [1]
        frags2 = pt.parse("PEPTIDE/2").fragment(ion_types=["b"])
        self.assertTrue(all(f.charge_state == 1 for f in frags2))

        # no charge annotation → last-resort [1]
        frags0 = pt.parse("PEPTIDE").fragment(ion_types=["b"])
        self.assertTrue(all(f.charge_state == 1 for f in frags0))

        # user-supplied charges always win
        frags_exp = pt.parse("PEPTIDE/3").fragment(ion_types=["b"], charges=[2])
        self.assertTrue(all(f.charge_state == 2 for f in frags_exp))


class TestFragmentMzPAF(unittest.TestCase):
    """Test mzPAF serialization of Fragment objects."""

    def test_b_ion_full(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.B, charge=2)
        self.assertEqual(frag.to_mzpaf(), "b7{PEPTIDE}^2")

    def test_b_ion_position(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.B, charge=2, position=3)
        self.assertEqual(frag.to_mzpaf(), "b3{PEP}^2")

    def test_y_ion_position(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.Y, charge=2, position=3)
        self.assertEqual(frag.to_mzpaf(), "y3{IDE}^2")

    def test_y_ion_charge_1(self):
        frag = pt.parse("PEPTIDE/1").frag(ion_type=pt.IonType.Y, charge=1, position=3)
        self.assertEqual(frag.to_mzpaf(), "y3{IDE}")

    def test_immonium(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.IMMONIUM, charge=2, position=3)
        self.assertEqual(frag.to_mzpaf(), "IP^2")

    def test_immonium_with_mod(self):
        frag = pt.parse("PEP[+10]TIDE/2").frag(ion_type=pt.IonType.IMMONIUM, charge=2, position=3)
        self.assertEqual(frag.to_mzpaf(), "IP[+10]^2")

    def test_internal_by(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.BY, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}^2")

    def test_internal_ay(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.AY, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}-CO^2")

    def test_internal_ax(self):
        # Regression: tacular>=1.1.0 corrected every non-"by" internal ion offset;
        # peptacular's mzPAF label table must track those corrected values.
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.AX, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}-2H^2")

    def test_internal_az(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.AZ, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}-HCONH2^2")

    def test_internal_bx(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.BX, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}+CO-2H^2")

    def test_internal_bz(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.BZ, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}-NH3^2")

    def test_internal_cx(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.CX, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}+CHNO^2")

    def test_internal_cy(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.CY, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}+NH3^2")

    def test_internal_cz(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.CZ, charge=2, position=(3, 5))
        self.assertEqual(frag.to_mzpaf(), "m3:5{PTI}^2")

    def test_internal_ion_labels_match_computed_mass(self):
        # The mzPAF label's implied mass delta must equal the actual computed delta
        # from "by" for every internal ion type (catches a label/value drift like the
        # one tacular>=1.1.0's fix exposed).
        import re

        MONO = {"H": 1.00782503223, "C": 12.0, "N": 14.00307400443, "O": 15.99491461957}

        def label_mass(label: str) -> float:
            if label is None:
                return 0.0
            total = 0.0
            for tok in re.findall(r"[+-][0-9]*[A-Za-z0-9]+", label):
                sign = 1 if tok[0] == "+" else -1
                mult_m = re.match(r"(\d*)(.*)", tok[1:])
                mult = int(mult_m.group(1)) if mult_m.group(1) else 1
                comp: dict[str, int] = {}
                for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", mult_m.group(2)):
                    if el:
                        comp[el] = comp.get(el, 0) + (int(n) if n else 1)
                total += sign * mult * sum(MONO[e] * n for e, n in comp.items())
            return total

        annot = pt.parse("PEPTIDE/1")
        by_mass = annot.frag(ion_type=pt.IonType.BY, charge=1, position=(3, 5)).mass
        for ion_type in (
            pt.IonType.AX,
            pt.IonType.AY,
            pt.IonType.AZ,
            pt.IonType.BX,
            pt.IonType.BY,
            pt.IonType.BZ,
            pt.IonType.CX,
            pt.IonType.CY,
            pt.IonType.CZ,
        ):
            frag = annot.frag(ion_type=ion_type, charge=1, position=(3, 5))
            mzpaf = frag.to_mzpaf(include_sequence=False)
            label = mzpaf[2:].split("^")[0] or None  # strip leading "m3:5", trailing charge
            actual_diff = frag.mass - by_mass
            implied_diff = label_mass(label)
            self.assertAlmostEqual(actual_diff, implied_diff, places=4, msg=f"{ion_type}: label {label!r}")

    def test_precursor(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.PRECURSOR, charge=2)
        self.assertEqual(frag.to_mzpaf(), "p^2")

    def test_no_sequence(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.Y, charge=2, position=3)
        self.assertEqual(frag.to_mzpaf(include_sequence=False), "y3^2")

    def test_adduct(self):
        frags = pt.parse("PEPTIDE/2").fragment(ion_types=["y"], charges=["Na:z+1"])
        self.assertIn("[M+Na]", frags[0].to_mzpaf())

    def test_isotope_c13(self):
        frags = pt.parse("PEPTIDE/1").fragment(ion_types=["y"], charges=[1], isotopes=[1])
        self.assertIn("+i", frags[0].to_mzpaf())

    def test_isotope_c13_x2(self):
        frags = pt.parse("PEPTIDE/1").fragment(ion_types=["y"], charges=[1], isotopes=[2])
        self.assertIn("+2i", frags[0].to_mzpaf())

    def test_isotope_custom(self):
        frags = pt.parse("PEPTIDE/1").fragment(ion_types=["y"], charges=[1], isotopes=[{"17O": 2}])
        self.assertIn("+2i17O", frags[0].to_mzpaf())

    def test_neutral_loss(self):
        frags = pt.parse("PEPTIDE/1").fragment(ion_types=["y"], charges=[1], neutral_deltas=["H2O"])
        labels = [f.to_mzpaf() for f in frags]
        # mzPAF reuses ProForma's own atom-then-count formula notation; "-H2O" is
        # the spec's own canonical water-loss example (section 4.5).
        self.assertTrue(any("-H2O" in label for label in labels))

    def test_numeric_neutral_loss_rejected(self):
        # mzPAF's neutral_loss grammar only accepts a chemical formula or a bracketed
        # reference-group name after the sign (section 4.5); there is no
        # representation for an arbitrary unnamed mass delta.
        frags = pt.parse("PEPTIDE/1").fragment(ion_types=["y"], charges=[1], deltas=[15.9949])
        with self.assertRaises(ValueError):
            frags[0].to_mzpaf()

    def test_adduct_repeat_count(self):
        # mzPAF section 4.7's own example: "[M+2Na] denotes an adduct ion with two
        # sodium atoms."
        annot = pt.parse("PEPTIDE/[Na:z+1^2]")
        frag = annot.frag(ion_type=pt.IonType.Y, position=3, charge=None)
        self.assertIn("[M+2Na]", frag.to_mzpaf())

    def test_multiple_adducts_are_alphabetized(self):
        # mzPAF section 4.7: "If there are multiple types of atoms/molecules,
        # alphabetical order SHOULD be followed, e.g. [M+2H+Na] rather than
        # [M+Na+2H]."
        annot = pt.parse("PEPTIDE/[Na:z+1,H:z+1^2]")
        frag = annot.frag(ion_type=pt.IonType.Y, position=3, charge=None)
        self.assertIn("[M+2H+Na]", frag.to_mzpaf())

    def test_serialize_format_default(self):
        frag = pt.parse("PEPTIDE/1").frag(ion_type=pt.IonType.Y, charge=1, position=3)
        self.assertEqual(frag.serialize(format="default"), str(frag))

    def test_serialize_format_mzpaf(self):
        frag = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.Y, charge=2, position=3)
        self.assertEqual(frag.serialize(format="mzpaf"), "y3{IDE}^2")

    def test_serialize_invalid_format(self):
        frag = pt.parse("PEPTIDE/1").frag(ion_type=pt.IonType.Y, charge=1, position=3)
        with self.assertRaises(ValueError):
            frag.serialize(format="invalid")

    def test_modification_in_sequence(self):
        frag = pt.parse("PEPT[Phospho]IDE/2").frag(ion_type=pt.IonType.B, charge=2, position=4)
        label = frag.to_mzpaf()
        self.assertEqual(label, "b4{PEPT[Phospho]}^2")

    def test_negative_charge_z_minus_1(self):
        frag = pt.parse("PEPTIDE/-2").frag(ion_type=pt.IonType.B, charge=-1, position=3)
        self.assertEqual(frag.charge_state, -1)
        self.assertGreater(frag.mz, 0)
        self.assertEqual(frag.to_mzpaf(), "b3{PEP}^1")  # mzPAF: charge is a bare magnitude, no minus sign
        self.assertEqual(frag.serialize(format="mzpaf"), "b3{PEP}^1")

    def test_negative_charge_z_minus_2(self):
        frag = pt.parse("PEPTIDE/-3").frag(ion_type=pt.IonType.B, charge=-2, position=3)
        self.assertEqual(frag.charge_state, -2)
        self.assertGreater(frag.mz, 0)
        self.assertEqual(frag.to_mzpaf(), "b3{PEP}^2")  # mzPAF: charge is a bare magnitude, no minus sign

    def test_negative_charge_y_ion(self):
        frag = pt.parse("PEPTIDE/-2").frag(ion_type=pt.IonType.Y, charge=-1, position=3)
        self.assertEqual(frag.charge_state, -1)
        self.assertGreater(frag.mz, 0)
        self.assertEqual(frag.to_mzpaf(), "y3{IDE}^1")  # mzPAF: charge is a bare magnitude, no minus sign

    def test_negative_charge_mz_less_than_positive(self):
        # Negative-mode b3 loses a proton; positive-mode adds one — so neg mz < pos mz
        frag_pos = pt.parse("PEPTIDE/2").frag(ion_type=pt.IonType.B, charge=1, position=3)
        frag_neg = pt.parse("PEPTIDE/-2").frag(ion_type=pt.IonType.B, charge=-1, position=3)
        self.assertAlmostEqual(frag_pos.mz - frag_neg.mz, 2 * 1.007276, places=4)


class TestFragmentRepeatedAdducts(unittest.TestCase):
    """Repeated charge carriers must keep their occurrence count through the fragment."""

    def test_two_sodium_carriers_not_collapsed(self):
        # Two 'Na:z+1' list entries are one carrier repeated twice; the count must not be
        # dropped when rebuilding the fragment's adduct tuple (regression: adjust_mass_mz
        # / adjust_comp built the tuple from _mods.keys(), collapsing it to a single Na).
        frag = pt.parse("PEPTIDE").frag("p", charge=["Na:z+1", "Na:z+1"])
        self.assertEqual(frag.charge_adducts.serialize(), "[Na:z+1,Na:z+1]")
        # Only one Na would leave neutral_mass off by a full sodium mass (~822 vs ~799).
        self.assertAlmostEqual(frag.neutral_mass, 799.35996, places=4)

    def test_two_sodium_carriers_mzpaf(self):
        # mzPAF must fold the repeat count into the adduct prefix: [M+2Na], not [M+Na].
        frag = pt.parse("PEPTIDE").frag("p", charge=["Na:z+1", "Na:z+1"])
        self.assertEqual(frag.to_mzpaf(include_sequence=False), "p[M+2Na]^2")

    def test_repeated_list_matches_occurance_syntax(self):
        # ['Na:z+1','Na:z+1'] (count=2) and 'Na:z+1^2' (occurance=2) are the same species.
        frag_list = pt.parse("PEPTIDE").frag("p", charge=["Na:z+1", "Na:z+1"])
        frag_occ = pt.parse("PEPTIDE").frag("p", charge="Na:z+1^2")
        self.assertEqual(frag_list.to_mzpaf(include_sequence=False), frag_occ.to_mzpaf(include_sequence=False))
        self.assertAlmostEqual(frag_list.neutral_mass, frag_occ.neutral_mass, places=6)

    def test_mixed_adducts_sorted_with_counts(self):
        # Alphabetical ordering (mzPAF 4.7) with a repeated proton: [M+2H+Na].
        frag = pt.parse("PEPTIDE").frag("p", charge=["H:z+1", "H:z+1", "Na:z+1"])
        self.assertEqual(frag.to_mzpaf(include_sequence=False), "p[M+2H+Na]^3")


class TestFragmentSequenceInternalCharge(unittest.TestCase):
    """Fragment.sequence must emit external charge, not charge_state (which includes internal)."""

    def test_sequence_uses_external_charge(self):
        # b2 carries external_charge=1 (one proton) plus an internal +1 formula charge.
        # Fragment.sequence must serialize /1 (external), not /2 (regression: it used
        # charge_state and double-counted the internal formula charge).
        frags = pt.parse("PE[Formula:CH2:z+1]PTIDE").set_charge(1).fragment(ion_types="b", charges=[1])
        b2 = next(f for f in frags if f.position == 2)
        self.assertEqual(b2.external_charge, 1)
        self.assertEqual(b2.charge_state, 2)
        self.assertEqual(b2.sequence, "PE[Formula:CH2:z+1]/1")


class TestFragmentLazyCompositionIonType(unittest.TestCase):
    """Fragment.composition (lazy path) must apply the fragment's ion-type offset.

    Regression: the lazy path (calculate_composition=False) returned the sub-sequence's
    *precursor* composition, ignoring the ion-type offset, so a b-ion's composition was heavier
    than its own mass by H2O. y-ions coincidentally matched (y neutral == precursor).
    """

    def _elem_mass(self, comp):
        return sum(el.get_mass() * n for el, n in comp.items())

    def test_b_ion_lazy_composition_matches_mass(self):
        # Use charge 0 so the neutral mass and the element-sum align exactly (no electron term).
        b4 = next(f for f in pt.parse("EVTKLE").fragment(ion_types=["b"], charges=[0], calculate_composition=False) if f.position == 4)
        self.assertAlmostEqual(self._elem_mass(b4.composition), b4.mass, places=3)

    def test_lazy_matches_eager_across_ion_types(self):
        annot = pt.parse("PEM[Oxidation]TIDE")
        for ion in ("b", "y", "a", "c", "x", "z"):
            lazy = annot.fragment(ion_types=[ion], charges=[1], calculate_composition=False)
            eager = annot.fragment(ion_types=[ion], charges=[1], calculate_composition=True)
            for fl, fe in zip(lazy, eager, strict=True):
                self.assertEqual(dict(fl.composition), dict(fe.composition), f"{ion} pos={fl.position}")

    def test_intact_ions_support_lazy_composition_and_sequence(self):
        annot = pt.parse("PEM[Oxidation]TIDE")
        for ion in ("p", "n"):
            lazy = annot.fragment(ion_types=[ion], charges=[1], calculate_composition=False)[0]
            eager = annot.fragment(ion_types=[ion], charges=[1], calculate_composition=True)[0]

            self.assertIsNone(lazy.position)
            self.assertEqual(lazy.sequence, "PEM[Oxidation]TIDE/1")
            self.assertEqual(dict(lazy.composition), dict(eager.composition))


if __name__ == "__main__":
    unittest.main()
