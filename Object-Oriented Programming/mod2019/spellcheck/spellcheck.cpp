#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

vector<string> make_bigram(const string &str);

vector<string> SINGLETON_VEC;


class Typo {
public:
    string original;
    vector<string> bigram;

//    best candidate
    string curr_candidate;
    float curr_similarity = 0;
    int curr_frequency = 0;

    explicit Typo(const string &s);
    void try_update_candidate(const string &other, vector<string> &other_bigrams, int other_freq);

private:
    float count_similarity(vector<string> &w1, vector<string> &w2);
};

Typo::Typo(const string &s)
{
    this->original = s;
    this->bigram = make_bigram(s);
// sort original word's bigram only once
    sort(this->bigram.begin(), this->bigram.end());
}

void Typo::try_update_candidate(const string &other, vector<string> &other_bigrams, int other_freq)
{
    float new_similarity = this->count_similarity(this->bigram, other_bigrams);

//  check similarity coefficient
    if (new_similarity > this->curr_similarity) {
        this->curr_candidate = other;
        this->curr_similarity = new_similarity;

//  proper float comparison: abs(a - b) < 0.01
    } else if (new_similarity == this->curr_similarity) {
        if (other_freq > this->curr_frequency) {
            this->curr_candidate = other;
            this->curr_similarity = new_similarity;
            this->curr_frequency = other_freq;

//      frequencies are equal
        } else if (other_freq == this->curr_frequency) {
//          new candidate is lexicographically shorter
            if (other > this->curr_candidate) {
                this->curr_candidate = other;
                this->curr_similarity = new_similarity;
                this->curr_frequency = other_freq;
            }
        }
    }
}

float Typo::count_similarity(vector<string> &w1, vector<string> &w2)
{
    // both bigram vectors are already sorted

    vector<string> intersect_vect;
    set_intersection(w1.begin(),w1.end(), w2.begin(),w2.end(), back_inserter(intersect_vect));
    float intersect_len = intersect_vect.size();

    vector<string> union_vect(w1.size() + w2.size());
    auto it = set_union(w1.begin(), w1.end(), w2.begin(), w2.end(), union_vect.begin());
    union_vect.resize(it - union_vect.begin());
    int union_len = union_vect.size();

    return intersect_len / union_len;
}

vector<string> make_bigram(const string &str)
{
    if (1 == str.length()) {
        SINGLETON_VEC.clear();
        SINGLETON_VEC.push_back(str);
        return SINGLETON_VEC;
    }
    int len = str.size() - 1;

    vector<string> bigrams;
    bigrams.reserve(len);
    for(int i = 0; i < len; ++i) {
        bigrams.push_back(string() + str[i] + str[i + 1]);
    }

    // sort bigram
    sort(bigrams.begin(), bigrams.end());

    return bigrams;
}



int main() {
//    string input = "prepearing\ngoverment\ncomming\nquickle\njouvenile\n";

    string input( (istreambuf_iterator<char>(cin)),(istreambuf_iterator<char>()) );
    input += "\n";

    vector<Typo> typos;

    stringstream ss(input);
    string buff;
    while (getline(ss, buff, '\n')) {
        typos.push_back(Typo(buff));
    }

    ifstream input_file("count_big.txt");

    string correct_word;
    int freq;
    while (input_file >> correct_word >> freq) {
        vector<string> correct_word_bigram = make_bigram(correct_word);
        for (Typo &w: typos) {
            w.try_update_candidate(correct_word, correct_word_bigram, freq);
        }
    }
    input_file.close();

    for (Typo &w: typos) {
//        cout << w.original << " -> " << w.curr_candidate << endl;
        cout << w.curr_candidate << endl;
    }

    return 0;
}