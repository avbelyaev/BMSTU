#include <iostream>
#include <map>


using namespace std;


//template <class T>
class Formula {
public:
    int arg1;
    int arg2;
    string op;

    Formula(int arg1_, int arg2_, const string& op_) {
        arg1 = arg1_;
        arg2 = arg2_;
        op = op_;
    }
};


//template <class T>
class Cell {


public:
    int data;
    bool isData;
    Cell* c1;
    Cell* c2;
    string op;

    Cell(int value_) {
        data = value_;
        isData = true;
    }

    Cell(Cell& c1_, Cell& c2_, string op_) {
        c1 = &c1_;
        c2 = &c2_;
        op = op_;
        isData = false;
    };

    Cell(Cell& c1_, int value, string op_) {
        c1 = &c1_;
        op = op_;
        isData = false;
    };

    int calc() {
        if (isData) {
            return data;
        }
        if ("+" == op) {
            return c1->calc() + c2->calc();
        }
        return 0;
    }

    Cell& operator= (const Cell& other_cell) {
        data = other_cell.data;
        return *this;
    }

    Cell& operator= (int other_data) {
        data = other_data;
        return *this;
    }

    friend Cell& operator+ (Cell& arg1, Cell& arg2) {
        Cell* c3 = new Cell(arg1, arg2, "+");
//        Cell c3(arg1, arg2, "+");
        return *c3;
    }

//    friend Cell& operator/ (const Cell& c, int value) {
//        Cell* c3 = new Cell(c, value, "/");
//        return *c3;
//    }
};



//template <class T>
class SuperCalc {
    map<int, Cell> data;

public:
    SuperCalc(int rows_, int cols_) {
//        map<int, int> m;
    };

//    Cell& operator() (int row, int col) {
//        return data[new pair(row, col)];
//    };
};


int main() {

//    SuperCalc s(1, 10);
//    Cell c1 = s(0, 10);
//
//    cout << c1.data << endl;
//
//    s(0, 10) = 5;
//    c1 = s(0, 10);
//    cout << c1.data << endl;

    Cell c1(1);
    Cell c2(2);

    Cell c3 = c1 + c2;
    cout << c3.calc() << endl;

    c1 = 5;

    cout << c3.calc() << endl;

    return 0;
}