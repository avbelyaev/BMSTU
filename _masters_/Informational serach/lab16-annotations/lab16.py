import re

import nltk
import pymorphy2

from collections import defaultdict

morph = pymorphy2.MorphAnalyzer()
nltk.download('punkt')  # required to split text into sentences

ARTICLES = ['cybersport.txt', 'fb.txt', 'telegram.txt']

SENTENCES_TO_OUTPUT = 4


class Sentence:
    def __init__(self, sentence: str):
        self.sent = sentence
        self.norm = normalize_sentence(sentence)
        self.words = self.norm.split(' ')

    def count_overall_frequency(self, wordcount: dict) -> int:
        freq = 0
        for w in self.words:
            freq += wordcount[w]
        return freq

    def count_avg_frequency(self, wordcount: dict) -> int:
        return self.count_overall_frequency(wordcount) // len(self.words)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.sent


class Article:
    def __init__(self, sentences: list):
        self.sents = sentences
        self.annt_overall = []
        self.annt_average = []
        self.wordcount = defaultdict(int)
        self._count_word_freq()

    def annotate(self):
        self._annotate_average()
        self._annotate_overall()

    def _annotate_overall(self):
        self.sents.sort(key=lambda sent: sent.count_overall_frequency(self.wordcount), reverse=True)
        self.annt_overall = self.sents[:SENTENCES_TO_OUTPUT]

    def _annotate_average(self):
        self.sents.sort(key=lambda sent: sent.count_avg_frequency(self.wordcount), reverse=True)
        self.annt_average = self.sents[:SENTENCES_TO_OUTPUT]

    def _count_word_freq(self):
        for s in self.sents:
            words = s.norm.split(' ')
            for w in words:
                self.wordcount[w] += 1


def normalize_sentence(sentence: str) -> str:
    tags_to_remove = ['NPRO', 'PRED', 'PREP', 'CONJ', 'PRCL', 'INTJ']
    normalized = []
    for w in sentence.split(' '):
        parsed = morph.parse(w)[0]
        if (parsed.tag.POS not in tags_to_remove) and 3 <= len(parsed.normal_form):
            normalized.append(parsed.normal_form)
    return ' '.join(normalized)


def clean_up_sentence(s: str) -> str:
    no_punct = re.sub(r'[^а-яё]', ' ', s.casefold())
    no_duplicate_spaces = re.sub(r'\s+', ' ', no_punct)
    return no_duplicate_spaces.strip()


def read_sentences(filename: str) -> list:
    text = open(filename, 'r', encoding='UTF-8').read()
    preprocessed = re.sub(r'[,:;-]', ' ', text)
    sentences = nltk.sent_tokenize(preprocessed)
    clean_sents = list(map(lambda s: clean_up_sentence(s), sentences))
    return list(filter(lambda s: len(s.split(' ')) >= 3, clean_sents))


def main():
    articles = []
    for article_filename in ARTICLES:
        sentences = read_sentences(article_filename)
        sentences = [Sentence(s) for s in sentences]
        articles.append(Article(sentences))

    [a.annotate() for a in articles]

    for a in articles:
        print('Overall')
        [print(f'\t{s}') for s in a.annt_overall]
        print('Average')
        [print(f'\t{s}') for s in a.annt_average]


if __name__ == '__main__':
    main()
