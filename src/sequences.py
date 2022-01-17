
import re


_COMPLEMENT_DICT = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
    'R': 'Y',
    'Y': 'R',
    'W': 'W',
    'S': 'S',
    'K': 'M',
    'M': 'K',
    'H': 'D',
    'V': 'B',
    'B': 'V',
    'D': 'H',
    'N': 'N',
}


_NON_SEQUENCE_PATTERN = r'[^{}]'.format(
    ''.join(_COMPLEMENT_DICT.keys())
)


def reverse_complement(seq):
    return ''.join(map(_get_complement_base, seq))[::-1]
# end def


def _get_complement_base(base):
    return _COMPLEMENT_DICT[base]
# end def


def verify_sequence(string):
    is_valid_sequence = re.search(_NON_SEQUENCE_PATTERN, string) is None
    return is_valid_sequence
# end def


def get_non_iupac_chars(string):
    return set(re.findall(_NON_SEQUENCE_PATTERN, string))
# end def
