import re
import warnings
from collections.abc import Generator


def get_regex_match_indices(input_str: str, regex_str: str | re.Pattern[str], offset: int = 0) -> Generator[int, None, None]:
    """
    Identify the starting indexes of occurrences of a given regex pattern within a string.

    :param input_str: The sequence in which to search.
    :type input_str: str
    :param regex_str: The regex pattern to search for.
    :type regex_str: str
    :param offset: An optional offset to add to each identified index. Default is 0.
    :type offset: int
    :return: A list of starting indexes where the regex pattern is found in the sequence.
    :rtype: List[int]

    .. code-block:: python

        >>> list(get_regex_match_indices("PEPTIDE", "(?=P)"))
        [0, 2]

        >>> list(get_regex_match_indices("PEPTIDE", re.compile("(?=E)")))
        [1, 6]

        >>> list(get_regex_match_indices("PEPTIDE", '(?=E)'))
        [1, 6]

        # More complex regex
        >>> list(get_regex_match_indices("PEPTIDE", "(?=P[ST])"))
        [2]

        >>> list(get_regex_match_indices("PPPP", "(?=PP)"))
        [0, 1, 2]

        >>> list(get_regex_match_indices("PEPCTIDE", "(?=C)"))
        [3]

        >>> list(get_regex_match_indices("PEPCTIDE", "(?<=P)(?=E)"))
        [1]

        >>> list(get_regex_match_indices("PEPCTIDE", "(?<=C)"))
        [4]


    """

    if not isinstance(regex_str, re.Pattern):
        regex_pattern = re.compile(regex_str)
    else:
        regex_pattern = regex_str

    # Use overlapping search
    pos = 0
    while pos < len(input_str):
        match = regex_pattern.search(input_str, pos)

        if match is None:
            break

        if match.start() != match.end():
            warnings.warn(
                message="The regex pattern has a non-zero-length match. Using the match end index (the position after the matched motif) as the cleavage site.",
                category=UserWarning,
                stacklevel=2,
            )
            yield match.end() + offset
        else:
            yield match.start() + offset

        # Advance by 1 to find overlapping matches
        pos = match.start() + 1


def get_regex_match_range(input_str: str, regex_str: str | re.Pattern[str], offset: int = 0) -> list[tuple[int, int]]:
    """
    Identify the starting indexes of occurrences of a given regex pattern within a string.

    :param input_str: The sequence in which to search.
    :type input_str: str
    :param regex_str: The regex pattern to search for.
    :type regex_str: str
    :param offset: An optional offset to add to each identified index. Default is 0.
    :type offset: int
    :return: A list of starting indexes where the regex pattern is found in the sequence.
    :rtype: List[int]

    .. code-block:: python

        >>> get_regex_match_range("PEPTIDE", "P")
        [(0, 1), (2, 3)]

        >>> get_regex_match_range("PEPTIDE", re.compile("E"))
        [(1, 2), (6, 7)]

        # More complex regex
        >>> get_regex_match_range("PEPTIDE", "P[ST]")
        [(2, 4)]

        >>> get_regex_match_range("PPPP", "PP")
        [(0, 2), (1, 3), (2, 4)]

    """

    if not isinstance(regex_str, re.Pattern):
        regex_pattern = re.compile(regex_str)
    else:
        regex_pattern = regex_str

    # Use overlapping search
    matches: list[tuple[int, int]] = []
    pos = 0
    while pos < len(input_str):
        match = regex_pattern.search(input_str, pos)

        if match is None:
            break

        matches.append((match.start() + offset, match.end() + offset))

        # Advance by 1 to find overlapping matches
        pos = match.start() + 1

    return matches
