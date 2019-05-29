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

using namespace std;

void get_tokens(const string &s, const unordered_set<char> &delimiters, vector<string> &tokens)
{
//    remove delimiters
    regex rm_delims("[^\\s\\w]+");
    string clean_up = regex_replace(s, rm_delims, " ");

//    set 1-len spaces
    regex rm_spaces("\\s+");
    string no_spaces = regex_replace(clean_up, rm_spaces, " ");

//    lowercase
    std::transform(no_spaces.begin(), no_spaces.end(), no_spaces.begin(), ::tolower);

//    tokenize
    stringstream str_stream(no_spaces);
    string buffer;
    while (str_stream >> buffer) {
        tokens.push_back(buffer);
    }
}