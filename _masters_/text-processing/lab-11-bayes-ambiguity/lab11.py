from collections import defaultdict
import re
from math import log
from typing import Optional
import pymorphy2


morph = pymorphy2.MorphAnalyzer()

AMBIGUOUS_WORD = morph.parse('рак')[0].normal_form


class Sentence:
    def __init__(self, sentence: str, vocabulary: list, klass: str = None):
        self.sent = sentence
        self.words = string_to_words(sentence)
        self.norm = normalize_words(self.words)
        self.klass = klass
        self.vector = {}
        self.apply_pos_features()
        self.apply_lexical_features(vocabulary)

    # apply part-of-speech features
    def apply_pos_features(self):
        amb_word_index = self.norm.index(AMBIGUOUS_WORD)
        word_left_2 = item_at_index_or_none(self.norm, amb_word_index - 2)
        word_left_1 = item_at_index_or_none(self.norm, amb_word_index - 1)
        word_right_1 = item_at_index_or_none(self.norm, amb_word_index + 1)
        word_right_2 = item_at_index_or_none(self.norm, amb_word_index + 2)
        self.vector = {
            # word to the left is present
            'l_exists': ft_word_exists(word_left_1),
            # word to the right is present
            'r_exists': ft_word_exists(word_right_1),
            # pos features
            # left   -2
            'l_2_noun': ft_noun(word_left_2),
            'l_2_adjf': ft_adjf(word_left_2),
            'l_2_verb': ft_verb(word_left_2),
            'l_2_comp': ft_comp(word_left_2),
            'l_2_numr': ft_numr(word_left_2),
            'l_2_advb': ft_advb(word_left_2),
            'l_2_prtf': ft_prtf(word_left_2),
            # left   -1
            'l_1_noun': ft_noun(word_left_1),
            'l_1_adjf': ft_adjf(word_left_1),
            'l_1_verb': ft_verb(word_left_1),
            'l_1_comp': ft_comp(word_left_1),
            'l_1_numr': ft_numr(word_left_1),
            'l_1_advb': ft_advb(word_left_1),
            'l_1_prtf': ft_prtf(word_left_1),
            # right  +1
            'r_1_noun': ft_noun(word_right_1),
            'r_1_adjf': ft_adjf(word_right_1),
            'r_1_verb': ft_verb(word_right_1),
            'r_1_comp': ft_comp(word_right_1),
            'r_1_numr': ft_numr(word_right_1),
            'r_1_advb': ft_advb(word_right_1),
            'r_1_prtf': ft_prtf(word_right_1),
            # right  +2
            'r_2_noun': ft_noun(word_right_2),
            'r_2_adjf': ft_adjf(word_right_2),
            'r_2_verb': ft_verb(word_right_2),
            'r_2_comp': ft_comp(word_right_2),
            'r_2_numr': ft_numr(word_right_2),
            'r_2_advb': ft_advb(word_right_2),
            'r_2_prtf': ft_prtf(word_right_2)
        }

    # takes dictionary and turns every word 'x' into feature: 'sentence contains "x"'
    def apply_lexical_features(self, vocabulary: list):
        for word in set(vocabulary):
            self.vector[word] = word in self.norm

    def satisfies_feature(self, feature_name: str) -> bool:
        return True == self.vector[feature_name]

    # extract features that evaluate to true(1) for given sentence
    def extract_active_features(self) -> list:
        active_features = dict(filter(lambda ft: True == ft[1], self.vector.items()))
        return list(active_features.keys())

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.sent



# =========== Feature predicates ===========
ft_noun = lambda w: 'NOUN' == tag(w)
ft_adjf = lambda w: 'ADJF' == tag(w)
ft_verb = lambda w: 'VERB' == tag(w)
ft_comp = lambda w: 'COMP' == tag(w)
ft_numr = lambda w: 'NUMR' == tag(w)
ft_advb = lambda w: 'ADVB' == tag(w)
ft_prtf = lambda w: 'PRTF' == tag(w)
ft_word_exists = lambda w: w is not None
# ==========================================


def tag(w: str) -> Optional[str]:
    return morph.parse(w)[0].tag.POS if w is not None else None


def item_at_index_or_none(lst: list, index: int) -> Optional[str]:
    try:
        return lst[index]
    except IndexError:
        return None


def read_sentences(filename: str) -> list:
    contents = open(filename, 'r').read()
    return list(filter(lambda s: '' != s, contents.split('\n')))


def string_to_words(s: str) -> list:
    clean = re.sub(r'[^а-яё]', ' ', s.casefold())
    return re.compile(r'\s+').split(clean)


# normalizes words and removes functional words
def normalize_words(words: list) -> list:
    tags_to_remove = ['NPRO', 'PRED', 'PREP', 'CONJ', 'PRCL', 'INTJ']
    normalized = []
    for w in words:
        parsed = morph.parse(w)[0]
        if (parsed.tag.POS not in tags_to_remove) and 3 <= len(parsed.normal_form):
            normalized.append(parsed.normal_form)
    return normalized


def extract_unique_features(sentences: list) -> list:
    return sentences[0].vector.keys()


def count_sents_in_class(sentences: list, klass: str) -> int:
    return len(list(filter(lambda s: klass == s.klass, sentences)))


def count_sents_in_class_satisfying_feature(sentences: list, klass: str, feature: str) -> int:
    sent_of_klass = list(filter(lambda s: klass == s.klass, sentences))
    sent_satisfying_ft = list(filter(lambda s: s.satisfies_feature(feature), sent_of_klass))
    return len(sent_satisfying_ft)


def train_bernoulli(klasses: list, sents: list) -> (dict, dict, dict):
    V = extract_unique_features(sents)
    N = len(sents)
    cond_prod = defaultdict(dict)
    prior = {}

    for klass in klasses:
        Nc = count_sents_in_class(sents, klass)
        prior[klass] = Nc / N

        for feature in V:
            Nct = count_sents_in_class_satisfying_feature(sents, klass, feature)
            cond_prod[feature][klass] = (Nct + 1) / (Nc + 2)
    return V, prior, cond_prod


def apply_bernoulli(klasses: list, V: list, prior: dict, cond_prob: dict, new_sent: Sentence) -> dict:
    Vd = new_sent.extract_active_features()
    score = {}

    for klass in klasses:
        score[klass] = log(prior[klass])
        for feature in V:
            if feature in Vd:
                score[klass] += log(cond_prob[feature][klass])
            else:
                score[klass] += log(1 - cond_prob[feature][klass])
    return score


def main():
    # make dictionary of normalized words from all train files
    train_1 = read_sentences('c1.txt')
    train_2 = read_sentences('c2.txt')
    vocabulary = set()
    for s in train_1 + train_2:
        normalized = normalize_words(string_to_words(s))
        vocabulary.update(normalized)
    vocabulary = list(vocabulary)

    # make sentences (Docs) for class 1
    klass_1 = 'cancer'
    sent_1 = list(map(lambda s: Sentence(s, vocabulary, klass_1), train_1))

    # make sentences (Docs) for class 2
    klass_2 = 'crayfish'
    sent_2 = list(map(lambda s: Sentence(s, vocabulary, klass_2), train_2))

    V, prior, cond_prod = train_bernoulli([klass_1, klass_2], sent_1 + sent_2)

    test_sentences = list(map(lambda s: Sentence(s, vocabulary), read_sentences('test.txt')))
    for sent in test_sentences:
        scores = apply_bernoulli([klass_1, klass_2], V, prior, cond_prod, sent)

        res = klass_1 if max(scores.values()) == scores[klass_1] else klass_2
        print(f'{res}: {sent}: cancer:{scores[klass_1]:.3f}, crayfish:{scores[klass_2]:.3f}')
        # print(f'{sent} & {scores[klass_1]:.3f} & {scores[klass_2]:.3f} \\\\')


if __name__ == '__main__':
    main()


