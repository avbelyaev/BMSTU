from xml.dom import minidom


TRAIN_DATA = 'news_eval_train.xml'
TEST_DATA = 'news_eval_test.xml'


TONE_POS = '+'
TONE_NEG = '-'
TONE_NEUT = '0'
TONALITY = [TONE_POS, TONE_NEG, TONE_NEUT]


class Cite:
    def __init__(self, _id: str, speech: str, evaluation: str, url: str):
        self._id = _id.strip()
        self.speech = speech.strip()
        self.evaluation = evaluation.strip()
        self.url = url.strip()

    def __str__(self):
        return f'{self.evaluation}: {self.speech}'

    def __repr__(self):
        return self.__str__()


def parse(filename: str) -> list:
    xmldoc = minidom.parse(filename)
    document = xmldoc.getElementsByTagName('document')[0]
    sentences = document.getElementsByTagName('sentence')

    citations = []
    for s in sentences:
        sid = s.attributes['id'].value
        speech = s.getElementsByTagName('speech')[0].childNodes[0].nodeValue
        evaluation = s.getElementsByTagName('evaluation')[0].childNodes[0].nodeValue
        url = s.getElementsByTagName('url')[0].childNodes[0].nodeValue

        citations.append(Cite(sid, speech, evaluation, url))

    return citations


def main():
    all_train = parse(TRAIN_DATA)
    train_cites = list(filter(lambda cite: cite.evaluation in TONALITY, all_train))
    print(f'train: {len(all_train)} -> {len(train_cites)}')

    all_test = parse(TEST_DATA)
    test_cites = list(filter(lambda cite: cite.evaluation in TONALITY, all_test))
    print(f'test : {len(all_test)} -> {len(test_cites)}')


if __name__ == '__main__':
    main()
