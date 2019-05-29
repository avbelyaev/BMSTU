#include <iostream>
#include <string>
#include <fstream>

using namespace std;

const string FILENAME = "count_big.txt";

int main() {
    string input = "prepearing\ngoverment\ncomming\nquickle\njouvenile\n";

//    string input( (istreambuf_iterator<char>(cin)),(istreambuf_iterator<char>()) );
//    input += "\n";
    cout << "starting" << endl;

    ifstream infilestream(FILENAME);

    string word;
    int freq;
    while (infilestream >> word >> freq) {
        cout << word << " " << freq << endl;
    }
    infilestream.close();
}