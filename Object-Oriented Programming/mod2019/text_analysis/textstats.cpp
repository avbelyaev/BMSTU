
#include <string>
#include <vector>
#include <unordered_set>
#include <map>
#include <iostream>
#include <regex>
#include <sstream>
#include <locale>
#include <set>
#include "textstats.hpp"

using namespace std;

void get_tokens(const string &s, const unordered_set<char> &delimiters, vector<string> &tokens)
{
    string delimiter_str;
    for (const char &delim: delimiters) {
        delimiter_str += delim;
    }

    locale loc;
    string str_lower;
    for(const char &c: s) {
        str_lower += tolower(c, loc);
    }

    stringstream stringStream(str_lower);
    string line;
    while(getline(stringStream, line))
    {
        size_t prev = 0, pos;
        while ((pos = line.find_first_of(delimiter_str, prev)) != std::string::npos)
        {
            if (pos > prev) {
                tokens.push_back(line.substr(prev, pos-prev));
            }
            prev = pos+1;
        }
        if (prev < line.length())
            tokens.push_back(line.substr(prev, std::string::npos));
    }
}

// make dict
void get_type_freq(const vector<string> &tokens, map<string, int> &freqdi)
{
    for (const string &token: tokens) {
        freqdi[token]++;
    }
}

void get_types(const vector<string> &tokens, vector<string> &wtypes)
{
//    remove duplicates
    set<string> s( tokens.begin(), tokens.end() );
    vector<string> unique_tokens( s.begin(), s.end() );
    wtypes = unique_tokens;

//    sort lexicographically
//    sort(wtypes.begin(), wtypes.end());
}

void get_x_length_words(const vector<string> &wtypes, int x, vector<string> &words)
{
    for (const string &s: wtypes) {
        if (s.length() >= x) {
            words.push_back(s);
        }
    }
}

void get_x_freq_words(const map<string, int> &freqdi, int x, vector<string> &words)
{
    for (const pair<const string, int> &entry: freqdi) {
        if (entry.second >= x) {
            words.push_back(entry.first);
        }
    }
    sort(words.begin(), words.end());
}

void get_words_by_length_dict(const vector<string> &wtypes, map<int, vector<string> > &lengthdi)
{
    for (const string &s: wtypes) {
        int len = s.length();
        lengthdi[len].push_back(s);
    }
}