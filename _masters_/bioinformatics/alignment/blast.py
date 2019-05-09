#!/usr/bin/env python3
import itertools

from _masters_.bioinformatics.alignment.smith_waterman import SubstMatrix, BLOSUM62_FILENAME

K_WORD_LEN = 3
T_THRESHOLD_VALUE = 12


def make_word_list(s: str, word_len) -> list:
    words = []

    i = 0
    while i + word_len <= len(s):
        words.append(s[i: i + word_len])
        i += 1
    return words


def list_possible_matches(k_words: list, subst_matrix: SubstMatrix) -> list:
    def match(s1: str, s2: str) -> int:
        assert len(s1) == len(s2)
        score = 0
        for i in range(len(s1)):
            score += subst_matrix.lookup(s1[i], s2[i])
        return score

    # make letter combinations (with length of k) of all proteins available
    all_combinations = list(itertools.combinations(subst_matrix.proteins, r=K_WORD_LEN))
    all_combinations = list(map(lambda tuple: ''.join(tuple), all_combinations))

    matches = {}
    for k_word in k_words:
        for combo in all_combinations:
            matches[k_word, combo] = match(k_word, combo)

    above_threshold = {k: v for k, v in matches.items() if v >= T_THRESHOLD_VALUE}
    print(f'matched combinations: {len(matches)} -> {len(above_threshold)}')

    # keep only combinations. k-words are not needed
    combinations_above_threshold = []
    for k, v in above_threshold.items():
        k_word = k[0]
        combo = k[1]
        combinations_above_threshold.append(combo)

    return combinations_above_threshold


def main():
    seq = 'PQGEFG'
    matrix = SubstMatrix(BLOSUM62_FILENAME)

    k_words = make_word_list(seq, K_WORD_LEN)
    print(f'words of len {K_WORD_LEN}: {k_words}')

    high_scoring_words = list_possible_matches(k_words, matrix)


if __name__ == '__main__':
    main()
