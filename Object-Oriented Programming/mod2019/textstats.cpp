//
// Created by anthony on 2019-05-29.
//

#include <string>
#include <vector>
#include <unordered_set>
#include <map>
#include <iostream>
#include <regex>
#include <sstream>
#include <set>

using namespace std;

void get_tokens(const string &s, const unordered_set<char> &delimiters, vector<string> &tokens)
{
//    remove delimiters
    regex rm_delims("[^\\s\\w]+");
    string no_delims = regex_replace(s, rm_delims, " ");

//    set 1-len spaces
    regex rm_spaces("\\s+");
    string no_spaces = regex_replace(no_delims, rm_spaces, " ");

//    lowercase
    std::transform(no_spaces.begin(), no_spaces.end(), no_spaces.begin(), ::tolower);

//    tokenize
    stringstream str_stream(no_spaces);
    string buffer;
    while (str_stream >> buffer) {
        tokens.push_back(buffer);
    }
}

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
    sort(wtypes.begin(), wtypes.end());
}