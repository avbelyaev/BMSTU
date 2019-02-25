import pymorphy2
import pymorphy2_dicts_ru
import re
from collections import Counter
import nltk
import time

morph = pymorphy2.MorphAnalyzer()

import time

delay = 10  # время задержки в секундах;

with open('c1.txt', 'r') as f:
    ls = []
    for line in f:
        lst = line.split()

        words = []
        for word in lst:
            p = morph.parse(re.sub(r'[^\w\s]', '', word).lower())[0]  # делаем разбор
            words.append(p)

        ls.append(words)

# print(ls)

for x in ls:  # для 1 строки
    # print(x)
    for item in x:  # для элементов строки
        pos_tag = item.tag.POS
        is_noun = pos_tag == 'NOUN'
        print(f"{item.word} {pos_tag} -> {is_noun}")
    #  if 'NOUN' in item.tag:
    # print (item.tag.POS)

# time.sleep(delay)


