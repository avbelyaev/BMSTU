import re
from collections import defaultdict

file_train = 'train.txt'
file_test = 'c1.txt'


WORD = 0
COUNT = 1


def read_words(filename: str) -> list:
    contents = open(filename, 'r').read().casefold()
    clean = re.sub(r'[^a-zа-яёЁ0-9]', ' ', contents)
    splitted = re.compile(r'\s+').split(clean)
    return splitted[:-1] # remove last empty word


def count_bigrams(word_list: list) -> dict:
    bigram_count = {}
    i = 1
    while i < len(word_list):
        bi = (word_list[i - 1], word_list[i])
        if bi in bigram_count:
            bigram_count[bi] += 1
        else:
            bigram_count[bi] = 1
        i += 1
    return bigram_count


def unigram_additive_model(train_words: list, test_words: list) -> (dict, float):
    # ----- train set -----
    word_count = defaultdict(int)
    for w in train_words:
        word_count[w] += 1

    N = len(train_words)
    V = len(word_count)
    print(f'N (num of words): {N}, V (num of types): {V}')

    word_probability = {}
    for w in word_count:
        word_probability[w] = (word_count[w] + 1) / (N + V)

    # ----- test set -----
    product = 1
    for w in test_words:
        if w in train_words:
            product *= word_probability[w]
        else:
            product *= 1 / (N + V)

    perplexity = product ** (-1 / len(test_words))

    return word_probability, perplexity



def bigram_witten_bell_model(train_words: list, test_words: list, unigram_model: dict) -> (dict, float):
    # ----- train set -----
    train_bigram_count = count_bigrams(train_words)
    word_count = defaultdict(int)
    for w in train_words:
        word_count[w] += 1

    print(f'num of bigrams: {len(train_bigram_count)}')

    # returns number of N-gram types which start with a given word 'prefix'
    N = lambda prefix: len(list(filter(lambda bi: bi[0] == prefix, list(train_bigram_count.keys()))))

    # Max likelihood: count(bigram) / count (bigram-prefix)
    Pml = lambda bigram: train_bigram_count[bigram] / word_count[bigram[0]]

    # Number of bigram histories (bigrams that end with given word)
    def E(word: str) -> int:
        end_with_word = list(filter(lambda bi: bi[0][1] == word, train_bigram_count.items()))
        return sum(w[1] for w in end_with_word)

    # 1 - lambda = Num-of-prefixes / (Num-of-prefixes + Total-occurrence)
    Lambda = lambda bigram: 1 - (N(bigram[0]) / (N(bigram[0]) + E(bigram[1])))

    bigram_probability = {}
    for bigram in train_bigram_count:
        bigram_probability[bigram] = Lambda(bigram) * Pml(bigram) + \
                                     (1 - Lambda(bigram)) * unigram_model[bigram[1]]

    # ----- test set -----
    V = len(word_count)
    test_bigrams = count_bigrams(test_words)
    product = 1
    for bi in test_bigrams:
        w1 = bi[0]
        w2 = bi[1]
        if bi in train_bigram_count:
            product *= bigram_probability[bi]

        elif w1 in train_words and w2 in train_words:
            product *= (1 - Lambda(bi)) * unigram_model[w1]

        elif w1 in train_words and w2 not in train_words:
            product *= (1 - Lambda(bi)) * unigram_model[w1]

        else:
            product *= 1 / (V ** 2)

    perplexity = product ** (-1 / len(test_words))

    return bigram_probability, perplexity


def main():
    train_set = read_words(file_train)
    test_set = read_words(file_test)

    # count models and perplexities
    uni_train_model, pp_unigram = unigram_additive_model(train_set, test_set)
    bi_train_model, pp_bigram = bigram_witten_bell_model(train_set, test_set, uni_train_model)
    print(f'Unigram perplexity: {pp_unigram}\n'
          f'Bigram  perplexity: {pp_bigram}')

    with open('out_unigrams.txt', 'w+') as out:
        for word, prob in uni_train_model.items():
            print(f'{word}\t\t&\t\t{prob} \\\\', file=out)

    with open('out_bigrams.txt', 'w+') as out:
        for word, prob in bi_train_model.items():
            print(f'{word[0]} {word[1]}\t\t&\t\t{prob} \\\\', file=out)


if __name__ == '__main__':
    main()
