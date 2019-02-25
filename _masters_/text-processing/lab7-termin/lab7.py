import re
import math
import pymorphy2
from collections import defaultdict


morph = pymorphy2.MorphAnalyzer()

# file_name = "c1.txt"
file_name = "cs-basics.txt"

TERMS_TO_OUTPUT = 50

TAG_ADJECTIVE = 'ADJF'
TAG_NOUN = 'NOUN'


def read_words(filename: str) -> list:
    contents = open(filename, 'r').read().casefold()
    clean = re.sub(r'[^а-яё]', ' ', contents)
    splitted = re.compile(r'\s+').split(clean)
    return splitted[1:-1] # remove empty words - first and last


def print_terms(terms: list, filename: str):
    with open(filename, 'w+') as out:
        for term in terms[:TERMS_TO_OUTPUT]:
            print(f'{term[0][0]} {term[0][1]} \\\\', file=out)


def make_bigrams(words: list) -> list:
    bigrams = []
    i = 1
    while i < len(words):
        bigrams.append((words[i - 1], words[i]))
        i += 1
    return bigrams


def normalize_words(words: list) -> list:
    tags_to_remove = ['NPRO', 'PRED', 'PREP', 'CONJ', 'PRCL', 'INTJ']
    normalized = []
    for w in words:
        parsed = morph.parse(w)[0]
        # leave only words that matter - nouns, verbs, adjectives, etc.
        if (parsed.tag.POS not in tags_to_remove) and 3 < len(parsed.normal_form):
            normalized.append(parsed.normal_form)

    return normalized


def find_terms(words: list) -> list:
    is_term = lambda bi: bi[0].tag.POS == TAG_ADJECTIVE \
                         and bi[1].tag.POS == TAG_NOUN

    parsed_words = list(map(lambda w: morph.parse(w)[0], words))
    parsed_bigrams = make_bigrams(parsed_words)

    # find pairs of ADJECTIVE-NOUN among bigrams
    potential_terms = list(filter(is_term, parsed_bigrams))
    # leave original words
    return list(map(lambda t: (t[0].word, t[1].word), potential_terms))


def count_mutual_info(term_count: dict, word_count: dict) -> list:
    N = sum(map(lambda w: word_count.get(w), word_count))
    term_measures = []
    for term in term_count:
        a = term[0]
        b = term[1]
        if a == 'берестяной':
            print(f'берестяной {term_count[term]}, {word_count[a]}, {word_count[b]}, N:{N}')
        MI = math.log2(term_count[term] * N / (word_count[a] * word_count[b]))
        term_measures.append((term, MI))

    # sort by Mutual Information measure DESC
    return sorted(term_measures, key=lambda term: term[1], reverse=True)


def main():
    words = read_words(file_name)

    normalized_words = normalize_words(words)
    word_count = defaultdict(int)
    for word in normalized_words:
        word_count[word] += 1


    normalized_terms = find_terms(normalized_words)
    term_count = defaultdict(int)
    for term in normalized_terms:
        term_count[term] += 1


    # get top terms by frequency
    terms_by_frequency = sorted(term_count.items(), key=lambda term: term[1], reverse=True)
    print_terms(terms_by_frequency, 'out_freq.txt')

    # get top terms by mutual information measure
    terms_by_mutual_info = count_mutual_info(term_count, word_count)
    print_terms(terms_by_mutual_info, 'out_mi.txt')


if __name__ == '__main__':
    main()
