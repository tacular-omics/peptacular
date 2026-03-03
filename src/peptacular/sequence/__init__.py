"""
The sequence module contains the functional API for peptacular.

All sequence functions support a single sequence/annotation input or a list of sequences/annotations for batch processing.
When multiple sequences/annotations are provided, a multiprocessing pool is used to parallelize the computations.

Sequecnes are converted to a proforma annotation and then thier properties are computed. If the output is a sequence,
it is converted back to a string representation via serialize.

This is great for processing large datasets of sequences/annotations efficiently, especially if only a few operations are
needed per sequence/annotation. If many operations are needed, it is more efficient to create a ProFormaAnnotation object once and
call the methods directly on that object. This will avoid the overhead of converting back and forth between string and object representations.
"""

from .basic import (
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
from .combinatoric import (
    combinations,
    combinations_with_replacement,
    permutations,
    product,
)
from .converters import (
    convert_casanovo_sequence,
    convert_diann_sequence,
    convert_ip2_sequence,
)
from .digestion import (
    cleavage_sites,
    digest,
    left_semi_digest,
    nonspecific_digest,
    right_semi_digest,
    semi_digest,
    simple_cleavage_sites,
    simple_digest,
)
from .fragmentation import fragment, frag
from .isotope import isotopic_distribution
from .mass_funcs import comp, mass, mz
from .mod_builder import (
    append_mods,
    modify,
    condense_static_mods,
    condense_to_peptidoform,
    extend_mods,
    filter_mods,
    from_ms2_pip,
    get_mods,
    pop_mods,
    remove_mods,
    set_mods,
    strip_mods,
    to_ms2_pip,
)
from .parallel import (
    get_available_start_methods,
    get_start_method,
    set_start_method,
)
from .properties import (
    aa_property_percentage,
    alpha_helix_percent,
    aromaticity,
    average_buried_area,
    beta_sheet_percent,
    beta_turn_percent,
    bulkiness,
    calc_property,
    calc_window_property,
    charge_at_ph,
    codons,
    coil_percent,
    flexibility,
    hplc,
    hydrophilicity,
    hydrophobicity,
    mutability,
    pi,
    polarity,
    property_partitions,
    recognition_factors,
    refractivity,
    secondary_structure,
    surface_accessibility,
    transmembrane_tendency,
)
from .subseqs import (
    coverage,
    find_subsequence_indices,
    is_subsequence,
    modification_coverage,
    percent_coverage,
)
from .transformations import (
    join,
    reverse,
    shift,
    shuffle,
    sort,
    span_to_sequence,
    split,
)

__all__ = [
    # basic
    "parse",
    "serialize",
    "parse_chimeric",
    "serialize_chimeric",
    "sequence_length",
    "is_ambiguous",
    "is_modified",
    "count_residues",
    "percent_residues",
    "annotate_ambiguity",
    "validate",
    "generate_random",
    # combinatoric
    "permutations",
    "combinations",
    "combinations_with_replacement",
    "product",
    # subseqs
    "is_subsequence",
    "find_subsequence_indices",
    "coverage",
    "percent_coverage",
    "modification_coverage",
    # digestion
    "left_semi_digest",
    "right_semi_digest",
    "semi_digest",
    "nonspecific_digest",
    "digest",
    "simple_digest",
    "cleavage_sites",
    "simple_cleavage_sites",
    # fragmentation
    "fragment",
    # mass_funcs
    "mass",
    "mz",
    "comp",
    # mod_builder
    "modify",
    "get_mods",
    "set_mods",
    "remove_mods",
    "pop_mods",
    "append_mods",
    "extend_mods",
    "condense_static_mods",
    "strip_mods",
    "filter_mods",
    "to_ms2_pip",
    "from_ms2_pip",
    "condense_to_peptidoform",
    # properties
    "calc_property",
    "hydrophobicity",
    "flexibility",
    "hydrophilicity",
    "surface_accessibility",
    "polarity",
    "mutability",
    "codons",
    "bulkiness",
    "recognition_factors",
    "transmembrane_tendency",
    "average_buried_area",
    "hplc",
    "refractivity",
    "calc_window_property",
    "charge_at_ph",
    "pi",
    "aa_property_percentage",
    "aromaticity",
    "secondary_structure",
    "alpha_helix_percent",
    "beta_sheet_percent",
    "beta_turn_percent",
    "coil_percent",
    "property_partitions",
    # transformations
    "reverse",
    "shuffle",
    "shift",
    "split",
    "span_to_sequence",
    "sort",
    "join",
    # converters
    "convert_ip2_sequence",
    "convert_diann_sequence",
    "convert_casanovo_sequence",
    # isotope
    "isotopic_distribution",
    # parallel
    "get_available_start_methods",
    "get_start_method",
    "set_start_method",
    "frag",
]
