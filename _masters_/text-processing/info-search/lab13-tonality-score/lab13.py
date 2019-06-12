from xml.dom import minidom
from xml.dom.minidom import Element

import nltk
from nltk.corpus import stopwords

from lab16 import clean_up_sentence, normalize_sentence

TRAIN_DATA = 'news_eval_train.xml'
TRAIN_DATA_2 = 'news_eval_train_SHORT.xml'
TEST_DATA = 'news_eval_test.xml'

TONALITY = ['+', '-', '0']

nltk.download('stopwords')
STOP_WORDS = stopwords.words('russian')


class Cite:
    def __init__(self, s: Element):
        self.speech = s.getElementsByTagName('speech')[0].childNodes[0].nodeValue.strip()
        self.evaluation = s.getElementsByTagName('evaluation')[0].childNodes[0].nodeValue.strip()
        self.tokenized = tokenize(self.speech)

    def __str__(self):
        return f'{self.evaluation}: {self.speech}'

    def __repr__(self):
        return self.__str__()


# [[1, 2, 3], [4], [5, 6]] -> [1, 2, 3, 4, 5, 6]
def flatten(lst: list) -> list:
    flat_list = []
    for sublist in lst:
        for item in sublist:
            flat_list.append(item)
    return flat_list


def tokenize(text: str) -> list:
    def remove_stopwords(sentence: str) -> list:
        words = sentence.split(' ')
        return list(filter(lambda w: w not in STOP_WORDS, words))

    sentences = nltk.sent_tokenize(text)
    clean_sents = list(map(lambda s: clean_up_sentence(s), sentences))
    norm_sents = list(map(lambda s: normalize_sentence(s), clean_sents))
    no_stopwords_sents = list(map(lambda s: remove_stopwords(s), norm_sents))
    return flatten(no_stopwords_sents)


def parse(filename: str) -> list:
    xmldoc = minidom.parse(filename)
    document = xmldoc.getElementsByTagName('document')[0]
    sentences = document.getElementsByTagName('sentence')

    citations = [Cite(s) for s in sentences]
    return citations


def main():
    all_train = parse(TRAIN_DATA_2)
    train_cites = list(filter(lambda cite: cite.evaluation in TONALITY, all_train))
    print(f'train: {len(all_train)} -> {len(train_cites)}')
    #
    # all_test = parse(TEST_DATA)
    # test_cites = list(filter(lambda cite: cite.evaluation in TONALITY, all_test))
    # print(f'test : {len(all_test)} -> {len(test_cites)}')


if __name__ == '__main__':
    main()
