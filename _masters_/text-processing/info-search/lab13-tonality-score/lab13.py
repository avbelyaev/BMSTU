from datetime import datetime
from xml.dom import minidom
import nltk
from nltk.corpus import stopwords
from sklearn.ensemble import AdaBoostClassifier
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn import metrics
from lab16 import clean_up_sentence, normalize_sentence

TRAIN_DATA = 'news_eval_train.xml'
TRAIN_DATA_2 = 'news_eval_train_SHORT.xml'
TEST_DATA = 'news_eval_test.xml'

TONALITY = ['+', '-', '0']

nltk.download('stopwords')
STOP_WORDS = stopwords.words('russian')


class Cite:
    def __init__(self, speech_: str, eval_: str):
        self.speech = speech_.strip()
        self.evaluation = eval_.strip()
        self.tokenized = tokenize(self.speech)
        self.tokenized_str = ' '.join(self.tokenized)

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


def classify(vectorizer, train_cites: list, test_cites: list):
    tokenized_train = [c.tokenized_str for c in train_cites]
    x_train = vectorizer.fit_transform(tokenized_train).toarray()
    y_train = [c.evaluation for c in train_cites]

    tokenized_test = [c.tokenized_str for c in test_cites]
    x_test = vectorizer.transform(tokenized_test).toarray()
    y_test = [c.evaluation for c in test_cites]

    classifier = AdaBoostClassifier(n_estimators=10)
    classifier.fit(x_train, y_train)
    predicted = classifier.predict(x_test)

    print(metrics.classification_report(y_test, predicted))


def parse(filename: str) -> list:
    xmldoc = minidom.parse(filename)
    document = xmldoc.getElementsByTagName('document')[0]
    sentences = document.getElementsByTagName('sentence')

    citations = []
    for s in sentences:
        speech = s.getElementsByTagName('speech')[0].childNodes[0].nodeValue
        evaluation = s.getElementsByTagName('evaluation')[0].childNodes[0].nodeValue.strip()
        if evaluation in TONALITY:
            citations.append(Cite(speech, evaluation))

    return citations


def main():
    print(f'parse train: {datetime.now()}')
    train_cites = parse(TRAIN_DATA)

    print(f'parse test:  {datetime.now()}')
    test_cites = parse(TEST_DATA)

    print(f'classify:    {datetime.now()}')
    # vectorizer = TfidfVectorizer()
    vectorizer = CountVectorizer()
    classify(vectorizer, train_cites, test_cites)
    print(f'classified:  {datetime.now()}')


if __name__ == '__main__':
    main()
