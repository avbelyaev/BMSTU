import re
import pymorphy2
from collections import Counter


morph = pymorphy2.MorphAnalyzer()

file_name = "a_dance_with_dragons.txt"

N_WORDS_TO_PRINT = 50


def read_words(filename: str) -> list:
    contents = open(filename, 'r').read().casefold()
    clean = re.sub(r'[^a-zа-яёЁ0-9]', ' ', contents)
    splitted = re.compile(r'\s+').split(clean)
    return splitted[:-1] # remove last empty word


def main():
    word_list = read_words(file_name)

    normalized = list(map(lambda word: morph.parse(word)[0].normal_form, word_list))
    sorted_by_rank = Counter(normalized).most_common()

    with open('out.txt', 'w+') as out:
        for word, count in sorted_by_rank:
            print(f'{word}\t & \t{count} \\\\', file=out)


if __name__ == '__main__':
    main()
