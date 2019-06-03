#include <iostream>
#include <map>


using namespace std;


enum Op { EQ, ADD, MUL };


//template <class T>
class Cell {


public:
    bool isUsed = false;
    int data;
    bool isData = false;
    Cell* l;
    Cell* r;
    Op op;

    Cell() {
        isUsed = false;
        l = nullptr;
        r = nullptr;
    }

    Cell(int value_) {
        data = value_;
        isData = true;
        isUsed = true;
        l = nullptr;
        r = nullptr;
    }

    Cell(Cell& c1_, Cell& c2_, Op op_) {
        l = &c1_;
        r = &c2_;
        op = op_;
        isData = false;
        isUsed = true;
    };

    Cell(Cell& c1_, int value, Op op_) {
        l = &c1_;
        r = new Cell(value);
        op = op_;
        isData = false;
        isUsed = true;
    };

    Cell(Cell& c1_, Op op_) {
        l = &c1_;
        r = nullptr;
        op = op_;
        isData = false;
        isUsed = true;
    }

    int calc() {
        if (isUsed) {
            if (isData) {
                return data;
            }
            switch (op) {
                case EQ:
                    return l->calc();
                case ADD:
                    return l->calc() + r->calc();
                case MUL:
                    return l->calc() * r->calc();

                default:
                    cout << "ERR CALC";
                    return 0;
            }
        }
    }

    Cell& operator= (Cell& other_cell) = default;

    Cell& operator= (int other_data) {
        data = other_data;
        isUsed = true;
        isData = true;
        return *this;
    }

    friend Cell& operator+ (Cell& a, Cell& b) {
        Cell* c3 = new Cell(a, b, ADD);
        return *c3;
    }

    friend Cell& operator* (Cell& a, int b) {
        Cell* c3 = new Cell(a, b, MUL);
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
//    Cell l = s(0, 10);
//
//    cout << l.data << endl;
//
//    s(0, 10) = 5;
//    l = s(0, 10);
//    cout << l.data << endl;

    Cell c1;
    Cell c2;

    Cell c3 = c1 + c2;
    c1 = 2;
    c2 = 3;
    cout << c3.calc() << endl;

    c1 = 5;
    cout << c3.calc() << endl;

    c1 = c2 * 2;
    cout << c3.calc() << endl;

    return 0;
}