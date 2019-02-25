import math
from collections import defaultdict

import nltk.data
import pymorphy2
import re


morph = pymorphy2.MorphAnalyzer()
nltk.download('punkt')   # required to split text into sentences

ARTICLES = [
    'art1-game.txt',
    'art2-siege.txt',
    'art3-beer.txt',
    'art4-jack.txt',
    'art5-jack-rule-63.txt'
]

FACTS = [
    'В рецензии на компьютерную игру критик пожаловался на то, что смерть заставляет начинать уровень заново.',
    'Под стенами осаждённой шведами русской крепости немцы побили шотландцев за пиво.',
    'Есть версия, что Джек Потрошитель был женщиной.'
]

VOCABULARY = set()

RESULTS_TO_OUTPUT = 10


# both document and query can be represented as vector
class Vectorizable:

    def __init__(self, sentence: str):
        self.sentence = sentence
        self.words = sentence.split(' ')
        self.vector = {}

    # weight of term in the document is its frequency
    def vectorize(self):
        for term in VOCABULARY:
            self.vector[term] = self.sentence.count(term)  # count == TermFrequency

    def norm(self) -> float:
        n = 0
        for word in self.vector.keys():
            n += self.vector[word] ** 2
        return math.sqrt(n)

    def __hash__(self):
        return hash(self.sentence)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.sentence



def normalize_sentence(sentence: str) -> str:
    tags_to_remove = ['NPRO', 'PRED', 'PREP', 'CONJ', 'PRCL', 'INTJ']
    normalized = []
    for w in sentence.split(' '):
        parsed = morph.parse(w)[0]
        if (parsed.tag.POS not in tags_to_remove) and 3 <= len(parsed.normal_form):
            normalized.append(parsed.normal_form)
    return ' '.join(normalized)


def clean_up_sentence(s: str) -> str:
    no_punct = re.sub(r'[^а-яё]', ' ', s.casefold())        # remove punctuation
    no_duplicate_spaces = re.sub(r'\s+', ' ', no_punct)     # remove duplicate spaces
    return no_duplicate_spaces.strip()                      # remove trailing whitespaces


def read_sentences(filename: str) -> list:
    text = open(filename, 'r').read()
    sentences = nltk.sent_tokenize(text)                                        # split into sentences
    clean_sentences = list(map(lambda s: clean_up_sentence(s), sentences))      # clean up
    normalized = list(map(lambda s: normalize_sentence(s), clean_sentences))    # normalize words
    return normalized


def scalar_product(v1: dict, v2: dict):
    # make sure vectors have same dimension
    assert len(v1) == len(v2)
    p = 0
    for w in v1.keys():
        p += v1[w] * v2[w]
    return p


def vector_space_model(query: Vectorizable, docs: list) -> list:
    cos_similarity = {}
    for doc in docs:
        cos_similarity[doc] = scalar_product(query.vector, doc.vector) / (query.norm() * doc.norm())

    sorted_by_weight = sorted(cos_similarity.items(), key=lambda doc: doc[1], reverse=True)
    return sorted_by_weight[:RESULTS_TO_OUTPUT]


def tf_idf(query: Vectorizable, docs: list):
    # TF == frequency of a term in a document or query (both are Vectorizable)
    def tf(term: str, vect: Vectorizable) -> int:
        return vect.words.count(term)

    # IDF == inverted frequency of a term among all documents in collection
    def idf(term: str) -> float:
        containing_term = list(filter(lambda doc: term in doc.words, docs))
        df = len(containing_term)
        return math.log(len(docs) / df)

    recalculated_doc_weights = {}
    for doc in docs:
        weight = defaultdict(int)
        for word in doc.words:
            weight[word] = tf(word, doc) * idf(word)
        recalculated_doc_weights[doc] = weight

    tip = list(recalculated_doc_weights)[:30]

    for word in query.words:
        query_weight[word] = tf(word, query) * idf(word)

    # now that we have constant weight of query, we can count weight of each document
    doc_weight = {}
    for doc in docs:
        doc_weight[doc]

    for doc in docs:

        query_weight = []
        for word in query.words:
            query_weight[word] = tf(word, doc)



def main():
    sentences = []
    for article_filename in ARTICLES:
        sentences.extend(read_sentences(article_filename))

    # remove sentences that are too short (less than 3 words)
    sentences = list(filter(lambda s: len(s.split(' ')) >= 3, sentences))

    # create vocabulary from words of all sentences
    for sentence in sentences:
        words = sentence.split(' ')
        VOCABULARY.update(words)

    # remove empty word that could appear by mistake :)
    if '' in VOCABULARY:
        VOCABULARY.remove('')

    docs = list(map(lambda s: Vectorizable(s), sentences))
    for d in docs:
        d.vectorize()

    clean_queries = list(map(lambda s: clean_up_sentence(s), FACTS))            # clean up
    norm_queries = list(map(lambda s: normalize_sentence(s), clean_queries))    # normalize
    queries = list(map(lambda q: Vectorizable(q), norm_queries))

    # find most relevant docs with vector-space model (no IDF)
    for q in queries:
        q.vectorize()
        matched_docs = vector_space_model(q, docs)
        print(q)
        for match in matched_docs:
            print(f'\t{match[1]:.3f} {match[0]}')

    # -//- with TF-IDF model
    for q in queries:
        matched_docs = tf_idf(q, docs)
        print(q)
        for match in matched_docs:
            print(f'\t{match[1]:.3f} {match[0]}')

    print('asdsa')


if __name__ == '__main__':
    main()

