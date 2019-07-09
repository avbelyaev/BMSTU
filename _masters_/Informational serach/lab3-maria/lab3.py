import numpy
import pymorphy2
import io
import re
import math


def scalar(v1: list, v2: list):
    p = 0
    i = 0
    while i < len(v1):
        p += v1[i] * v2[i]
        i += 1
    return p


def cosine(v1, v2):
    return scalar(v1, v2) / math.sqrt(scalar(v1, v1) * scalar(v2, v2))


file = io.open('input.txt', encoding='utf-8')
fileLines = list(file)

morph = pymorphy2.MorphAnalyzer()

allSentences = []
for line in fileLines:
    allSentences += line.split(".")

for i in range(0, len(allSentences)):
    whitespaced = re.sub('[^а-яё]', ' ', allSentences[i].lower())
    splited = re.split('\s+', whitespaced)
    allSentences[i] = ' '.join(splited).strip()

for sent in allSentences:
    if '' == sent or sent.isspace() or len(sent.split(' ')) < 3:
        while sent in allSentences:
            allSentences.remove(sent)


dictionary = []
functors_pos = {'CONJ', 'ADV-PRO', 'CONJ', 'PART', 'PR', 'PREP'}


def should_be_added(word: str) -> bool:
    normalized = morph.parse(re.sub(r'[^\w\s]', '', word).lower())[0].normal_form
    word_pos = morph.parse(normalized)[0].tag.POS
    return (normalized != "") and (normalized != '\n') and (word_pos not in functors_pos)


for sent in allSentences:
    sent = sent.lower()
    for word in re.split('[^а-яё]', sent):
        normalized = morph.parse(re.sub(r'[^\w\s]', '', word).lower())[0].normal_form
        if should_be_added(word) and normalized not in dictionary:
            dictionary.append(normalized)


# print (dictionary)
wordCount = len(dictionary)
sentCount = len(allSentences)
count = {}
matrix = numpy.zeros((sentCount, wordCount))
for i in range(0, sentCount):
    current_sentence = allSentences[i].lower()

    splitted_sentence = re.split('[^а-яё]', current_sentence)

    for word in splitted_sentence:
        normalized = morph.parse(re.sub(r'[^\w\s]', '', word).lower())[0].normal_form

        if should_be_added(word) and (normalized in dictionary):
            matrix[i, dictionary.index(normalized)] += 1
            if normalized in count:
                count[normalized] += 1
            else:
                count[normalized] = 1

# print (count)

# requests = [
#     'Прежде существовали электрические суда (туеры), получавшие энергию от контактной сети, подобно троллейбусам',
#     'ЮАР стала первой в мире страной, добровольно отказавшейся от использования ядерного оружия',
#     'Христианские миссионеры использовали деревянных кукол для посрамления безнравственных женщин'
# ]

while 1:
    print("введите запрос или пустую строку чтобы закончить")
    request = input()
    if request == "":
        exit(0)

    reqVec = numpy.zeros((1, wordCount))

    for word in re.split('[^а-яё]', request.lower()):
        word = morph.parse(word)[0].normal_form
        word_pos = morph.parse(word)[0].tag.POS
        if word != "" and word != '\n' and word in dictionary and word_pos != 'PREP':
            reqVec[0, dictionary.index(word)] += 1
    # print(reqVec)
    cosDistance = []
    vecTfIdf = numpy.zeros((1, wordCount))
    for word in re.split('[^а-яё]', request.lower()):
        word = morph.parse(word)[0].normal_form
        word_pos = morph.parse(word)[0].tag.POS
        if word != "" and word != '\n' and word in dictionary and word_pos != 'PREP':
            vecTfIdf[0, dictionary.index(word)] = reqVec[0, dictionary.index(word)] * math.log10(
                sentCount / count[word])
            # vecTfIdf[0, dictionary.index(word)] = math.log(1+reqVec[0, dictionary.index(word)])*math.log10(sentCount/count[word])

    for i in range(0, sentCount):
        # res = 1 - distance.cosine(reqVec[0],matrix[i])
        # res = 1 - distance.cosine(vecTfIdf[0], matrix[i])

        request_vector = reqVec[0].tolist()
        request_only_nulls = all(v == 0 for v in request_vector)

        sentence_vector = matrix[i].tolist()
        sentence_only_nulls = all(v == 0 for v in sentence_vector)

        if (request_only_nulls or sentence_only_nulls):
            print('nulls here!')

        resCusom = cosine(request_vector, sentence_vector)
        cosDistance.append((resCusom, i))


    cosDistance.sort(reverse=True)

    resDistances = []
    resSentences = []

    for i in range(0, sentCount):
        # print(cosDistance[i][0], " - " , allSentences[cosDistance[i][1]])
        resDistances.append(cosDistance[i][0])
        resSentences.append(allSentences[cosDistance[i][1]])


    print(request)

    i = 0
    while i < 10:
        print(f'{resDistances[i]} \t{resSentences[i]}')
        i += 1

    # df = pd.DataFrame({'cosDistances': resDistances, 'sentences': resSentences})
    # df.to_excel('filename.xlsx', sheet_name="sheet1", index=False)
