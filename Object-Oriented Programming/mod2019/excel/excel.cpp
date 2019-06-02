#include <iostream>

using namespace std;


template <class T>
class Cell {


public:
    T data;

    Cell() {
        // empty
    }

    Cell<T>& operator= (const Cell<T>& other_cell) {
        data = other_cell.data;
        return *this;
    }

    Cell<T>& operator= (const T other_data) {
        data = other_data;
        return *this;
    }

};





template <class T>
class SuperCalc {
    int rows, cols;
    Cell<T> **data;

public:
    SuperCalc<T>(int rows_, int cols_) {
        rows = rows_;
        cols = cols_;

        data = new Cell<T> *[rows];
        for(int i = 0; i < rows; ++i){

            data[i] = new Cell<T> [cols];
            for(int j = 0; j < cols; ++j){
                data[i][j] = *(new Cell<T>());
            }
        }
    };

    Cell<T>& operator() (int row, int col) {
        return data[row][col];
    };
};


int main() {

    SuperCalc<int> s(1, 10);
    Cell<int> c1 = s(0, 10);

    cout << c1.data << endl;

    s(0, 10) = 5;
    c1 = s(0, 10);
    cout << c1.data << endl;


    return 0;
}