#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <string.h>
#include <fstream>
#include <set>
#include <map>

using namespace std;

int get_tokens(const string s, vector<string> &tokens)
{
    string newstr;
    newstr.clear();
    for (int i = 0, mark = 0; i < s.size(); i++, mark = 0) {
        if (s[i] >= 'a' && s[i] <= 'z') newstr += s[i];
        if (s[i] == '\n') mark = 1;
        if (!newstr.empty() && mark == 1) {
            tokens.push_back(newstr);
            newstr.clear();
        }
    }
    return tokens.size();
}

float dice_coefficient(set<string> a, set<string> b)
{
    int intersection = 0, join;
    for(set<string>::iterator IT = b.begin(); IT != b.end(); IT++) intersection += a.count((*IT));
    join = a.size() + b.size() - intersection;
    float dice = (float)intersection / (float)join;
    return dice;
}

set<string> make_bi(string a)	//Назовём биграммой слова его подстроку длиной в две буквы.
{
    set<string> retset;
    for (int k = 0; k < a.length() - 1; k++) retset.insert(a.substr(k, 2));
    return retset;
}

int main()
{
//Tests:
    //string intext = "prepearing\ngoverment\ncomming\nquickle\njouvenile\n";               //1
    //string intext = "beatiful\ntogegether\nenvolving\nilness\nepidemia\n";              //2
    //string intext = "anual\nsincerly\ncluching\nmentaly\nballons\ngirle\nbrethless\n";    //4
    //string intext = "finaly\nexausted\ngrabed\npubliclly\nexcelent\ncontageous\nbegining\nnobady\nhappenin\ninnecessary\n"; //5
    //string intext = "epidemy\nbycicle\ndamadged\nstollen\ndeliceaus\npreventions\nhollidays\nfamilly\nradios\nsympthoms\n"; //6
    //string intext = "affraid\nmeasurment\nappologized\nsimptoms\natribute\npannic\nsincerly\nhuredly\nstoped\nacused\n";      //7


    string intext( (istreambuf_iterator<char>(cin)),(istreambuf_iterator<char>()) );
    intext += "\n";
    map    < string, string > in_out;		// map  < incorrect_word, correct_word >
    vector < set   < string > > check_word_bi;	// vector of sets of incorrects_word's bigramms
    vector < string > check_word;		// incorrect_words
    vector < int > freq_vect;			// correct_word's frequencies
    vector < float > dice_vect;			// dice_coefficients
    string word;
    int i, word_num = get_tokens(intext, check_word), freq;
    float newdice;

    for (i = 0; i < word_num; i++) {
        check_word_bi.push_back(make_bi(check_word[i]));   //makes a bigramm(set) for check_word and pushes it into vector
        dice_vect.push_back(0);
        freq_vect.push_back(0);
    }

    ifstream infile;
    infile.open("count_big.txt"); 		//text format: correct_word(word) word_frequency(number) \n

    while (infile >> word >> freq) {
        set <string> word_bi = make_bi(word);
        for (i = 0; i < word_num; i++) {
            newdice = dice_coefficient(word_bi, check_word_bi[i]);
            if (newdice > dice_vect[i]) {
                dice_vect[i] = newdice;
                in_out[check_word[i]] = word;
                freq_vect[i] = freq;
            } else {
                if (newdice == dice_vect[i] && freq > freq_vect[i]) {
                    in_out[check_word[i]] = word;
                    freq_vect[i] = freq;
                }
            }
        }
    }

    for (i = 0; i < word_num; i++) cout << in_out[check_word[i]] << endl;
    infile.close();
    return 0;
}

