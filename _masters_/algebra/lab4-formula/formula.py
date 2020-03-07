from abc import ABC
from typing import List

counter = 0


def nextIndex() -> int:
    global counter
    counter += 1
    return counter


class Symbol(ABC):
    pass


class Connective(Symbol):
    def __init__(self, neutral: bool):
        super().__init__()
        self.neutral = neutral

    @staticmethod
    def conj() -> 'Connective':
        """ конъюнкция *"""
        return Connective(True)

    @staticmethod
    def disj() -> 'Connective':
        """ дизъюнкция +"""
        return Connective(False)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '*' if self.neutral else '+'


class Impl(Symbol):
    @staticmethod
    def imply() -> 'Impl':
        """ импликация >> => """
        return Impl()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '=>'


class Eq(Symbol):
    @staticmethod
    def equivalent() -> 'Eq':
        """ эквивалентность % =="""
        return Eq()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '=='


class Neg(Symbol):
    @staticmethod
    def negate() -> 'Neg':
        """ отрицание -"""
        return Neg()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '-'


class Formula:
    def __init__(self, sym: Symbol, args: 'List[Formula]'):
        """ use factory methods instead! """
        self.symbol = sym
        self.args = args
        self.isLiteral = True if 0 == len(args) else False

    @staticmethod
    def new(sym: Symbol, args: 'List[Formula]') -> 'Formula':
        return Formula(sym, args)

    @staticmethod
    def literal(sym: Symbol) -> 'Formula':
        return Formula(sym, [])

    def __mul__(self, other):
        return Formula.new(Connective.conj(), [self, other])

    def __add__(self, other):
        return Formula.new(Connective.disj(), [self, other])

    def __neg__(self):
        return Formula.new(Neg.negate(), [self])

    def __rshift__(self, other):
        return Formula.new(Impl.imply(), [self, other])

    def __mod__(self, other):
        return Formula.new(Eq.equivalent(), [self, other])

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.isLiteral:
            return f'{self.symbol}'

        if isinstance(self.symbol, Neg):
            return f'{self.symbol}{self.args[0]}'

        left = f'({self.args[0]})'
        if self.args[0].isLiteral:
            left = f'{self.args[0]}'

        right = f'({self.args[1]})'
        if self.args[1].isLiteral:
            right = f'{self.args[1]}'

        return f'{left} {self.symbol} {right}'


class Var(Symbol):
    def __init__(self, val: str):
        """ use factory methods instead! """
        self.val = val
        self.index = nextIndex()  # уникальный номер переменной

    @staticmethod
    def new(val: str) -> 'Var':
        return Var(val)

    def __mul__(self, other):
        return Formula.new(Connective.conj(), [Formula.literal(self), Formula.literal(other)])

    def __add__(self, other):
        return Formula.new(Connective.disj(), [Formula.literal(self), Formula.literal(other)])

    def __neg__(self):
        return Formula.new(Neg.negate(), [Formula.literal(self)])

    def __rshift__(self, other):
        return Formula.new(Impl.imply(), [Formula.literal(self), Formula.literal(other)])

    def __mod__(self, other):
        return Formula.new(Eq.equivalent(), [Formula.literal(self), Formula.literal(other)])

    def __repr__(self):
        return f'{self.val}'


def main():
    a = Var.new('a')
    b = Var.new('b')
    c = Var.new('c')
    d = Var.new('d')
    e = Var.new('e')

    f = ((a + b) * -c) % (d >> e)
    print(f)


if __name__ == '__main__':
    main()
