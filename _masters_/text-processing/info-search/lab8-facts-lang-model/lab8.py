from lab3 import read_sentences, clean_up_sentence, normalize_sentence

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

COLLECTION = []

RESULTS_TO_OUTPUT = 5

LAMBDA = 0.9


def lang_model(query: str, docs: list) -> list:
    def collection_model (term: str) -> float:
        return COLLECTION.count(term) / len(COLLECTION)

    def document_model(term: str, document: list) -> float:
        if 0 == document.count(term):
            return COLLECTION.count(term) / len(COLLECTION)     # cf / cs
        else:
            return document.count(term) / len(document)         # tf(t,d) / len(d)

    similarity = {}
    for doc in docs:
        p = 1
        for term in query.split(' '):
            p *= (1 - LAMBDA) * collection_model(term) + LAMBDA * document_model(term, doc)
            assert 0 != p    # make sure smoothing works and we dont have 'p(term) == 0'
        similarity[doc] = p

    sorted_by_weight = sorted(similarity.items(), key=lambda doc: doc[1], reverse=True)
    return sorted_by_weight[:RESULTS_TO_OUTPUT]


def main():
    sentences = []
    for article_filename in ARTICLES:
        sentences.extend(read_sentences(article_filename))

    # create collection from words of all sentences
    for sentence in sentences:
        COLLECTION.extend(sentence.split(' '))

    # remove empty word that could appear by mistake :)
    if '' in COLLECTION:
        COLLECTION.remove('')

    clean_queries = [clean_up_sentence(s) for s in FACTS]
    queries = [normalize_sentence(s) for s in clean_queries]

    # smoothing: add words from query to collection to avoid 'p(term) == 0'
    for q in queries:
        COLLECTION.extend(q.split(' '))

    for q in queries:
        matched_docs = lang_model(q, sentences)
        print(q)
        for match in matched_docs:
            print(f'\t{match[1]} {match[0]}')


if __name__ == '__main__':
    main()
