from math import sqrt
import numpy as np

np.set_printoptions(precision=3)

DICT = ['w0', 'w1', 'w2']
DOCS = [['w0', 'w1', 'w1'],
        ['w0', 'w1', 'w2'],
        ['w0', 'w2', 'w2']]

ITERS = 5
THEMES = 3


def rand_init_matrix(cols: int, rows: int) -> list:
    matrix = np.random.rand(rows, cols)
    m_sum = sum(np.sum(matrix, 0).tolist())
    matrix /= m_sum

    matrix = matrix.transpose().tolist()
    for i in range(len(matrix)):
        s = sum(matrix[i])
        for j in range(len(matrix[i])):
            matrix[i][j] /= sqrt(s)
    return np.array(matrix).transpose().tolist()


def em_step(m_words: list, m_docs: list, docs: list):
    def n_dwt(thm_id: int, word_id: int, doc_id: int) -> float:
        p = m_words[word_id][thm_id] * m_docs[thm_id][doc_id]
        s = 0
        for i in range(THEMES):
            s += m_words[word_id][i] * m_docs[i][doc_id]
        return docs[doc_id].count(DICT[word_id]) * (p / s)

    def n_wt(thm_id: int, word_id: int) -> float:
        s = 0
        for i in range(len(docs)):
            s += n_dwt(thm_id, word_id, i)
        return s

    def n_td(thm_id: int, doc_id: int) -> float:
        s = 0
        doc = docs[doc_id]
        for i in range(len(doc)):
            s += n_dwt(thm_id, i, doc_id)
        return s

    def n_t(thm_id: int) -> float:
        s = 0
        for i in range(len(DICT)):
            s += n_wt(thm_id, i)
        return s

    def n_d(doc_id: int) -> float:
        s = 0
        for i in range(THEMES):
            s += n_td(i, doc_id)
        return s

    m_fi = np.zeros((len(docs), THEMES))
    for i in range(len(docs)):
        for j in range(THEMES):
            m_fi[i][j] = n_wt(j, i) / n_t(j)

    m_psi = np.zeros((THEMES, len(DICT)))
    for i in range(THEMES):
        for j in range(len(DICT)):
            m_psi[i][j] = n_td(i, j) / n_d(j)

    return m_fi, m_psi


def main():
    m_words = rand_init_matrix(THEMES, len(DICT))
    m_docs = rand_init_matrix(len(DOCS), THEMES)

    print(f'words:\n{np.array(m_words)}')
    print(f'docs:\n{np.array(m_docs)}')

    for i in range(ITERS):
        m_words, m_docs = em_step(m_words, m_docs, DOCS)
        print(f'{i}: words:\n{np.array(m_words)}')
        print(f'{i}: docs:\n{np.array(m_docs)}')


if __name__ == '__main__':
    main()
