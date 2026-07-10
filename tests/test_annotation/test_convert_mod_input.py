from peptacular.annotation.annotation import convert_moddict_input


def test_convert_mod_input_dict():
    assert convert_moddict_input({"Oxidation": 1}) == {"Oxidation": 1}
    assert convert_moddict_input({"Oxidation": 2, "Phospho": 1}) == {
        "Oxidation": 2,
        "Phospho": 1,
    }
    # Test non-string keys
    assert convert_moddict_input({1: 1}) == {"1": 1}


def test_convert_mod_input_str():
    assert convert_moddict_input("Oxidation") == {"Oxidation": 1}
    assert convert_moddict_input("Phospho") == {"Phospho": 1}
    assert convert_moddict_input("") == {"": 1}


def test_convert_mod_input_number():
    assert convert_moddict_input(1) == {"+1": 1}
    assert convert_moddict_input(-1) == {"-1": 1}
    assert convert_moddict_input(1.5) == {"+1.5": 1}
    assert convert_moddict_input(-1.5) == {"-1.5": 1}
    assert convert_moddict_input(0) == {"+0": 1}


def test_convert_mod_input_iterable():
    assert convert_moddict_input(["Oxidation", "Oxidation"]) == {"Oxidation": 2}
    assert convert_moddict_input(["Oxidation", "Phospho"]) == {
        "Oxidation": 1,
        "Phospho": 1,
    }
    assert convert_moddict_input((1, 2, 1)) == {"1": 2, "2": 1}


def test_convert_mod_input_empty():
    assert convert_moddict_input({}) == {}
    assert convert_moddict_input([]) == {}


def test_convert_mod_input_none():
    assert convert_moddict_input(None) == {}


def test_convert_mod_input_other():
    # Test with an object that doesn't match any criteria
    class SomeObj:
        pass

    assert convert_moddict_input(SomeObj()) == {}


def test_convert_mod_input_strips_whitespace_before_interning():
    # Differently-padded-but-equivalent inputs must collapse onto the same
    # interned string so peptides sharing a modification actually share memory
    # and the downstream parser's @lru_cache, instead of each paying for its
    # own copy and its own cache miss.
    padded = convert_moddict_input("  Oxidation  ")
    clean = convert_moddict_input("Oxidation")
    assert padded == clean == {"Oxidation": 1}
    assert next(iter(padded)) is next(iter(clean))

    padded_dict = convert_moddict_input({"  Phospho  ": 2})
    assert padded_dict == {"Phospho": 2}

    padded_iterable = convert_moddict_input(["  Oxidation  ", "Oxidation"])
    assert padded_iterable == {"Oxidation": 2}
