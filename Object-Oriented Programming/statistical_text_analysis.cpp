#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <unordered_set>
#include <set>
#include <map>
#include <string.h>
#include <algorithm>
#include "textstats.hpp"

using namespace std;

const unordered_set<char> delimiters {
	'~', '!', '@', '#', '$', '%', '^', '&', '*', '(', ')', '_',
	'+', '-', '=', '`', '{', '}', '[', ']', '|', '\\', ':', ';',
	'\"', '\'', ',', '.', '/', '<', '>', '?', ' ', '\t', '\n'
};
// get_tokens разбивает текст, закодированный в ASCII, на слова,
// с переводом всех букв к нижнему регистру.
// Слово -- это последовательность символов, не являющихся разделителями.
void get_tokens(
    // [ВХОД] Текстовая строка.
    const string &s,
    // [ВХОД] Множество символов-разделителей.
    const unordered_set<char> &delimiters,
    // [ВЫХОД] Последовательность слов.
    vector<string> &tokens
)
{
    string buffer = s;
    int buff_size = buffer.size();
    string newstr;
    newstr.clear();

    for (int i = 0; i < s.size(); i++) if (buffer[i] >= 'A' && buffer[i] <= 'Z') buffer[i] = tolower(buffer[i]);
    if (delimiters.find(buffer[buff_size]) == delimiters.end()) buffer.append(".");
    for (int i = 0, mark = 0; i < buff_size; i++, mark = 0) {
        if (delimiters.find(buffer[i]) == delimiters.end()) newstr += buffer[i];
        if (delimiters.find(buffer[i+1]) != delimiters.end()) mark = 1;
        if (!newstr.empty() && mark == 1) {
            tokens.push_back(newstr);
            newstr.clear();
        }
    }
}


// get_type_freq составляет частотный словарь текста, в котором
// для каждого слова указано количество его вхождений в текст.
void get_type_freq(
    // [ВХОД] Последовательность слов.
    const vector<string> &tokens, 
    // [ВЫХОД] Частотный словарь
    // (ключи -- слова, значения -- количества вхождений).
    map<string, int> &freqdi
)
{
    for (int i = 0; i < tokens.size(); i++) freqdi[tokens[i]]++;
}


// get_types составляет список уникальных слов, встречающихся в
// тексте. Список должен быть отсортирован лексикографически.
void get_types(
    // [ВХОД] Последовательность слов.
    const vector<string> &tokens, 
    // [ВЫХОД] Список уникальных слов.
    vector<string> &wtypes
)
{
    set<string> newset;
    for (int i = 0; i < tokens.size(); i++) newset.insert(tokens[i]);

    //for (int i = 0; i < newset.size(); i++) {
      for (auto& f : newset){
        wtypes.push_back(f);
    }
    //sort (wtypes.begin(), wtypes.end());
}


// get_x_length_words формирует из списка уникальных слов новый
// список, в который попадают только те слова, длина которых
// не меньше x символов.
void get_x_length_words(
    // [ВХОД] Список уникальных слов.
    const vector<string> &wtypes, 
    // [ВХОД] Минимальная длина слова.
    int x, 
    // [ВЫХОД] Список слов, длина которых не меньше x.
    vector<string> &words
)
{
    for (int i = 0; i < wtypes.size(); i++) if (wtypes[i].size() >= x) words.push_back(wtypes[i]);
}


// get_x_freq_words формирует лексикографически отсортированный
// список слов, взятых из частотного словаря, которые встречаются
// в тексте не меньше x раз.
void get_x_freq_words(
    // [ВХОД] Частотный словарь
    const map<string, int> &freqdi, 
    // [ВХОД] Минимальное количество вхождений.
    int x,
    // [ВЫХОД] Список слов, встречающихся не меньше x раз.
    vector<string> &words
)
{
    for (auto& f: freqdi) if (f.second >= x) words.push_back(f.first);
    sort (words.begin(), words.end());
}


// get_words_by_length_dict составляет словарь, в котором каждый
// ключ -- это длина слова, а значение -- это список слов заданной длины.
void get_words_by_length_dict(
    // [ВХОД] Список уникальных слов.
    const vector<string> &wtypes, 
    // [ВЫХОД] Словарь распределения слов по длинам.
    map<int, vector<string> > &lengthdi
)
{
    int j, mark;
    vector<int>sizes;
    vector <string> strings;

    for (int i = 0; i < wtypes.size(); i++) {
        for (j = 0, mark = 0; j < i; j++) {
            if (wtypes[i].size() == wtypes[j].size()) {
                mark = 1;
                break;
            }
        }
        if (0 == mark) sizes.push_back(wtypes[i].size());
    }
    for (int i = 0; i < sizes.size(); i++) {
        for (j = 0, mark = 0; j < wtypes.size(); j++) {
            if (wtypes[j].size() == sizes[i]) {
                strings.push_back(wtypes[j]);
                mark = 1;
            }
        }
        if (mark == 1) {
            lengthdi[sizes[i]] = strings;
            strings.clear();
        }
    }
}


int main()
{
	// Ввод текста из стандартного потока ввода
	// (при вводе текста из консоли в конце нужно нажать Ctrl-D).
        //string text( (istreambuf_iterator<char>(cin)), (istreambuf_iterator<char>()) );
	
//string text = "In Key MaP the key values ARE generally used TO idenTify the elements thE mAp valuES store the content associated to key\n";

	// Итератор для вывода слов через пробел.
	ostream_iterator<string> owords(cout, " ");

	// Разбиение текста на слова.
	vector<string> tokens;
	get_tokens(text, delimiters, tokens);
	copy(tokens.begin(), tokens.end(), owords);
	cout << endl;

	// Составление частотного словаря.
	map<string, int> freqdi;
	get_type_freq(tokens, freqdi);
	for (auto p : freqdi) {
		cout << "(" << p.first << " => " << p.second << ") ";
	}
	cout << endl;

	// Составление списка уникальных слов.
	vector<string> wtypes;
	get_types(tokens, wtypes);
	copy(wtypes.begin(), wtypes.end(), owords);
	cout << endl;

	// Вычисление средней длины уникальных слов.
	int sum = 0;
	for (auto w : wtypes) sum += w.length();
	int med = tokens.size() > 0 ? (sum + wtypes.size() - 1)/wtypes.size() : 0;

	// Формирование списка слов, длина которых не ниже средней.
	vector<string> long_wtypes;
	get_x_length_words(wtypes, med, long_wtypes);
	copy(long_wtypes.begin(), long_wtypes.end(), owords);
	cout << endl;

	// Формирование списка слов, встречающихся больше одного раза.
	vector<string> multi_wtypes;
	get_x_freq_words(freqdi, 2, multi_wtypes);
	copy(multi_wtypes.begin(), multi_wtypes.end(), owords);
	cout << endl;

	// Составление словаря распределения слов по длинам.
	map<int, vector<string> > lengthdi;
	get_words_by_length_dict(wtypes, lengthdi);
	for (auto p : lengthdi) {
		cout << p.first << " => ";
		copy(p.second.begin(), p.second.end(), owords);
		cout << endl;
	}

	return 0;
}

