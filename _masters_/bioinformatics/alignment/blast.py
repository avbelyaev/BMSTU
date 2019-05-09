#!/usr/bin/env python3
import itertools

from _masters_.bioinformatics.alignment.smith_waterman import SubstMatrix, BLOSUM62_FILENAME

QUERY_SEQUENCE = 'PQGEQG'

# just random 'words'
PROTEIN_DATABASE = [
    'DPPEGVVJGK', 'FGHKJNBFLKJG', 'QSGHJHIHKLHU', 'JKFGHSLAA', 'DPVQGZZFQJV',
    'DKHQJKKBF', 'GFDLWPQK', 'AJKHJEWBFJHWE', 'KOLBASA', 'EGDKASFJL', 'ASJKGHDSKQ'
]

SUBST_MATRIX = SubstMatrix(BLOSUM62_FILENAME)
K_WORD_LEN = 3
T_THRESHOLD_VALUE = 12


class HighscorePair:
    def __init__(self, origin_substring: str, matching_sbstring: str):
        self.origin = origin_substring
        self.match = matching_sbstring
        self.score = 0
        # hsp fields
        self.db_entry = None
        self.hsp_query = None
        self.hsp_db_entry = None
        self.accumulated_score = 0

    def calc_match_score(self):
        assert len(self.origin) == len(self.match)
        for i in range(len(self.origin)):
            self.score += SUBST_MATRIX.lookup(self.origin[i], self.match[i])

    def extend_to_highscore_segment_pair(self):
        q_from, db_from = self._extend_left()
        q_to, db_to = self._extend_right()
        self.hsp_query = QUERY_SEQUENCE[q_from: q_to]
        self.hsp_db_entry = self.db_entry[db_from: db_to]

    def _extend_left(self) -> (int, int, int):
        i = QUERY_SEQUENCE.index(self.origin)
        j = self.db_entry.index(self.match)

        score = 0
        # <---
        while score >= 0 and i >= 0 and j >= 0:
            # increase score
            self.accumulated_score += score
            score = SUBST_MATRIX.lookup(QUERY_SEQUENCE[i], self.db_entry[j])
            i -= 1
            j -= 1
        return i + 1, j + 1

    def _extend_right(self):
        i = QUERY_SEQUENCE.index(self.origin) + 1
        j = self.db_entry.index(self.match) + 1

        score = 0
        # <---
        while score >= 0 and i < len(QUERY_SEQUENCE) and j < len(self.db_entry):
            # increase score
            self.accumulated_score += score
            score = SUBST_MATRIX.lookup(QUERY_SEQUENCE[i], self.db_entry[j])
            i += 1
            j += 1
        return i - 1, j - 1

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'[{self.origin}]{QUERY_SEQUENCE} <-> [{self.match}]{self.db_entry}'


def make_word_list(s: str, word_len) -> list:
    words = []

    i = 0
    while i + word_len <= len(s):
        words.append(s[i: i + word_len])
        i += 1
    return words


def list_possible_matches(k_words: list) -> list:
    # make letter combinations (with length of k) of all proteins available
    all_combinations = list(itertools.combinations(SUBST_MATRIX.proteins, r=K_WORD_LEN))
    all_combinations = list(map(lambda tuple: ''.join(tuple), all_combinations))

    pairs = []
    for k_word in k_words:
        for combo in all_combinations:
            pairs.append(HighscorePair(k_word, combo))

    [pair.calc_match_score() for pair in pairs]
    above_threshold = list(filter(lambda pair: pair.score >= T_THRESHOLD_VALUE, pairs))
    print(f'reducing matching pairs: {len(pairs)} -> {len(above_threshold)}')
    return above_threshold


def find_database_matches(high_scoring_words: list) -> list:
    for highscore_word in high_scoring_words:
        for db_entry in PROTEIN_DATABASE:
            if highscore_word.match in db_entry:
                highscore_word.db_entry = db_entry

    return list(filter(lambda e: e.db_entry is not None, high_scoring_words))


# hsp == high-scoring segment pair
# hsp-extending step is from original BLAST, not BLAST2!
# extend input sequences's match to db-entry
def extend_matches_to_hsp(matches: list) -> list:
    return [m.extend_to_highscore_segment_pair() for m in matches]


def blast():
    print(f'query sequence: {QUERY_SEQUENCE}, subst matrix: BLOSUM62')

    k_words = make_word_list(QUERY_SEQUENCE, K_WORD_LEN)
    print(f'words of len {K_WORD_LEN}: {k_words}')

    high_scoring_words = list_possible_matches(k_words)
    print(f'high-scroing words: {high_scoring_words}')

    exact_matches = find_database_matches(high_scoring_words)
    print(f'exact database matches (word, db-entry): {exact_matches}')

    hsps = extend_matches_to_hsp(exact_matches)
    print(f'HSPs: {hsps}')


if __name__ == '__main__':
    blast()
