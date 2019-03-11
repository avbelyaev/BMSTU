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

RESULTS_TO_OUTPUT = 5


# both document and query can be represented as vector
class Vectorizable:

    def __init__(self, sentence: str):
        self.sentence = sentence
        self.words = sentence.split(' ')
        self.vector = {}
        self.tfidf_vector = {}

    # weight of term in the document is its frequency
    def vectorize(self):
        for term in VOCABULARY:
            self.vector[term] = self.sentence.count(term)  # count == TermFrequency

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
    text = open(filename, 'r', encoding='UTF-8').read()
    preprocessed = re.sub(r'[,:;-]', ' ', text)
    sentences = nltk.sent_tokenize(preprocessed)                                # split into sentences
    clean_sents = list(map(lambda s: clean_up_sentence(s), sentences))          # clean up
    no_short_sents = list(filter(lambda s: len(s.split(' ')) >= 3, clean_sents))# remove short sentences
    normalized = list(map(lambda s: normalize_sentence(s), no_short_sents))     # normalize words
    return normalized


def scalar_product(v1: dict, v2: dict):
    p = 0
    for w in v1.keys():
        p += v1[w] * v2[w]
    return p


def norm(vect: dict) -> float:
    n = 0
    for word in vect.keys():
        n += vect[word] ** 2
    return math.sqrt(n)


def vector_space_model(query: Vectorizable, docs: list) -> list:
    cos_similarity = {}
    for doc in docs:
        multiplied_norms = norm(query.vector) * norm(doc.vector)
        cos_similarity[doc] = scalar_product(query.vector, doc.vector) / multiplied_norms

    sorted_by_weight = sorted(cos_similarity.items(), key=lambda doc: doc[1], reverse=True)
    return sorted_by_weight[:RESULTS_TO_OUTPUT]


def tf_idf(query: Vectorizable, docs: list) -> list:
    # TF == frequency of a term in a document or query (both are Vectorizable)
    def tf(term: str, vect: Vectorizable) -> int:
        return vect.words.count(term)

    # IDF == inverted frequency of a term among all documents in collection
    def idf(term: str) -> float:
        containing_term = list(filter(lambda doc: term in doc.words, docs))
        df = len(containing_term)
        if 0 == df: df = 1
        return math.log(len(docs) / df)

    # count TFIDF weight for each document
    for doc in docs:
        weight = defaultdict(int)
        for word in doc.words:
            weight[word] = tf(word, doc) * idf(word)
        doc.tfidf_vector = weight

    # count TFIDF for a query
    weight = defaultdict(int)
    for word in query.words:
        weight[word] = tf(word, query) * idf(word)
    query.tfidf_vector = weight

    # now just apply cosine similarity like for vector-space model
    cos_similarity = {}
    for doc in docs:
        multiplied_norms = norm(query.tfidf_vector) * norm(doc.tfidf_vector)
        cos_similarity[doc] = scalar_product(query.tfidf_vector, doc.tfidf_vector) / multiplied_norms

    sorted_by_weight = sorted(cos_similarity.items(), key=lambda doc: doc[1], reverse=True)
    return sorted_by_weight[:RESULTS_TO_OUTPUT]


def main():
    sentences = []
    for article_filename in ARTICLES:
        sentences.extend(read_sentences(article_filename))

    # create vocabulary from words of all sentences
    for sentence in sentences:
        VOCABULARY.update(sentence.split(' '))

    # remove empty word that could appear by mistake :)
    if '' in VOCABULARY:
        VOCABULARY.remove('')

    clean_queries = list(map(lambda s: clean_up_sentence(s), FACTS))            # clean up
    norm_queries = list(map(lambda s: normalize_sentence(s), clean_queries))    # normalize

    # add words from query to dictionary
    for q in norm_queries:
        VOCABULARY.update(q.split(' '))

    queries = list(map(lambda q: Vectorizable(q), norm_queries))
    [q.vectorize() for q in queries]

    docs = list(map(lambda s: Vectorizable(s), sentences))
    [d.vectorize() for d in docs]

    print('Vector space model')
    for q in queries:
        q.vectorize()
        matched_docs = vector_space_model(q, docs)
        print('\hline \n \\bold{' + q.sentence + '} \\\\ \n \hline')
        for match in matched_docs:
            print(f'\t{match[1]:.3f} {match[0]} \\\\')

    print('TF-IDF')
    for q in queries:
        matched_docs = tf_idf(q, docs)
        print('\hline \n \\bold{' + q.sentence + '} \\\\ \n \hline')
        for match in matched_docs:
            print(f'\t{match[1]:.3f} {match[0]} \\\\')


if __name__ == '__main__':
    main()

