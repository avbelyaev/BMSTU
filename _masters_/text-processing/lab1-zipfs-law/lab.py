import re
import pymorphy2
import matplotlib
import numpy as np
from collections import Counter

matplotlib.use('TkAgg') # macos fix

morph = pymorphy2.MorphAnalyzer()

file_name = "a_dance_with_dragons.txt"


def draw_plot(x_values, y_values):
    from matplotlib import pyplot as plt
    plt.plot(np.array(x_values), np.array(y_values))

    plt.xlabel('Rank')
    plt.ylabel('Frequency')
    plt.title("Zipf's law for George Martin's 'A dance with dragons'")
    plt.show()


def main():
    with open(file_name, encoding="utf-8", mode="r") as f:
        contents = f.read().lower()

        clean_text = re.sub("[^a-zа-я0-9]", " ", contents)
        splitted = re.compile("\s+").split(clean_text)
        normalized = list(map(lambda word: morph.parse(word)[0].normal_form, splitted))
        sorted_by_rank = Counter(normalized).most_common()

        with open('out.txt', 'w+') as out:
            for word, count in sorted_by_rank:
                print(f'{word}\t{count}', file=out)

        frequencies_list = list(map(lambda word_count_tuple: word_count_tuple[1], sorted_by_rank))
        draw_plot(list(range(1, len(frequencies_list) + 1)), frequencies_list)


def zipfs_number():
    out_file_name = 'out.txt'
    with open(file=out_file_name, mode='r') as f:
        contents = f.read()
        wordcounts = contents.split('\n')
        filtered_empty_words =  list(filter(lambda wc: '' != wc, wordcounts))
        counts = list(map(lambda wc: int(wc.split('\t')[1]), filtered_empty_words))
        enumerated_counts = list(enumerate(counts))
        indexed_from_1 = list(map(lambda index_count: (index_count[0] + 1, index_count[1]), enumerated_counts))

        zipfs_numbers = list(map(lambda rank_count: rank_count[0] * rank_count[1], indexed_from_1))
        with open('out_zipf.txt', 'w+') as out:
            for zipf in zipfs_numbers:
                print(f'{zipf}', file=out)


        from matplotlib import pyplot as plt
        x = list(range(1, len(zipfs_numbers) + 1))
        y = zipfs_numbers
        plt.plot(np.array(x), np.array(y))
        plt.show()


if __name__ == '__main__':
    main()
