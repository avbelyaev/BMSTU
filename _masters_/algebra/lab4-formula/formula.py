import copy
from abc import ABC
from typing import List, Any

# @formatter:off
counter = 0
def nextIndex() -> int:
    global counter
    counter += 1
    return counter


class Symbol(ABC):
    pass


class Connective(Symbol):
    def __init__(self, neutral: bool):
        self.neutral = neutral

    @staticmethod
    def conjAnd() -> 'Connective': return Connective(True)

    @staticmethod
    def disjOr() -> 'Connective': return Connective(False)

    def __repr__(self): return self.__str__()
    def __str__(self): return '*' if self.neutral else '+'


class Impl(Symbol):
    @staticmethod
    def imply() -> 'Impl': return Impl()

    def __repr__(self): return self.__str__()
    def __str__(self): return '=>'


class Eq(Symbol):
    @staticmethod
    def equivalent() -> 'Eq': return Eq()

    def __repr__(self): return self.__str__()
    def __str__(self): return '<=>'


class Neg(Symbol):
    @staticmethod
    def negate() -> 'Neg': return Neg()

    def __repr__(self): return self.__str__()
    def __str__(self): return '!'


class Formula:
    def __init__(self, sym: Symbol, args: 'List[Any]'):
        """ use factory methods instead! """
        self._symbol = sym
        self._args = args
        self.isLiteral = True if 0 == len(args) else False
        self.canonize()

    def canonize(self):
        if isinstance(self._symbol, Eq):                     # инвариант-2: a <=> b === (a => b) AND (b => a)
            self._symbol = Connective.conjAnd()
            a = copy.deepcopy(self._args[0])
            b = copy.deepcopy(self._args[1])
            self._args[0] = Formula.new(Impl.imply(), [a, b])
            self._args[1] = Formula.new(Impl.imply(), [b, a])

        if isinstance(self._symbol, Impl):                   # инвариант-2: a => b === !a OR b
            self._symbol = Connective.disjOr()
            self._args[0] = Formula.new(Neg.negate(), [self._args[0]])

    @property
    def isConjunctionAnd(self) -> bool: return isinstance(self._symbol, Connective) and self._symbol.neutral

    @property
    def isDisjunctionOr(self) -> bool: return isinstance(self._symbol, Connective) and not self._symbol.neutral

    def addArgument(self, var: 'Var'):
        self._args.append(Formula.literal(var))
        self._args.sort(key=lambda arg: arg._symbol.value if arg.isLiteral else 'dummy')

    def mergeWithFormula(self, other: 'Formula'):
        self._args.extend(other._args)
        self._args.sort(key=lambda arg: arg._symbol.value if arg.isLiteral else 'dummy')

    @staticmethod
    def new(sym: Symbol, args: 'List[Any]') -> 'Formula':
        return Formula(sym, args)

    @staticmethod
    def literal(sym: Symbol) -> 'Formula':
        return Formula(sym, [])

    def __mul__(self, other):
        # (a * b) * c === *(a, b, c)
        if self.isConjunctionAnd and isinstance(other, Var):
            self.addArgument(other)
            return self

        # (a * b) * (c * d) === *(a, b, c, d)
        if self.isConjunctionAnd and isinstance(other, Formula) and other.isConjunctionAnd:
            self.mergeWithFormula(other)
            return self

        return Formula.new(Connective.conjAnd(), [self, other])

    def __add__(self, other):
        # (a + b) + c === +(a, b, c)
        if self.isDisjunctionOr and isinstance(other, Var):
            self.addArgument(other)
            return self

        # (a + b) + (c + d) === +(a, b, c, d)
        if self.isDisjunctionOr and isinstance(other, Formula) and other.isDisjunctionOr:
            self.mergeWithFormula(other)
            return self

        return Formula.new(Connective.disjOr(), [self, other])

    def __neg__(self): return Formula.new(Neg.negate(), [self])
    def __rshift__(self, other): return Formula.new(Impl.imply(), [self, other])
    def __mod__(self, other): return Formula.new(Eq.equivalent(), [self, other])

    def __repr__(self): return self.__str__()

    def __str__(self):
        if self.isLiteral:
            return f'{self._symbol}'

        left = f'({self._args[0]})'
        if self._args[0].isLiteral:
            left = f'{self._args[0]}'

        if isinstance(self._symbol, Neg):
            return f'{self._symbol}{left}'

        s = f'{self._symbol}('
        for arg in self._args:
            s += f'{arg} '
        s = s[0: len(s) - 1]
        s += ')'
        return s


class Var(Symbol):
    def __init__(self, val: str):
        """ use factory methods instead! """
        self.value = val
        self.index = nextIndex()  # уникальный номер переменной

    @staticmethod
    def new(val: str) -> 'Var':
        return Var(val)

    def __mul__(self, other):
        # a * (b * c) === *(a, b, c)
        if isinstance(other, Formula) and other.isConjunctionAnd:
            other.addArgument(self)
            return other
        return Formula.new(Connective.conjAnd(), [Formula.literal(self), Formula.literal(other)])

    def __add__(self, other):
        # a + (b + c) === +(a, b, c)
        if isinstance(other, Formula) and other.isDisjunctionOr:
            other.addArgument(self)
            return other
        return Formula.new(Connective.disjOr(), [Formula.literal(self), Formula.literal(other)])

    def __neg__(self): return Formula.new(Neg.negate(), [Formula.literal(self)])
    def __rshift__(self, other): return Formula.new(Impl.imply(), [Formula.literal(self), Formula.literal(other)])
    def __mod__(self, other): return Formula.new(Eq.equivalent(), [Formula.literal(self), Formula.literal(other)])

    def __repr__(self): return f'{self.value}'



class Normalizer:
    @staticmethod
    def toNegationNormalForm(formula: Formula, b: bool) -> Formula:
        f = copy.deepcopy(formula)

        return f


def main():
    a = Var.new('a')
    b = Var.new('b')
    c = Var.new('c')
    d = Var.new('d')
    e = Var.new('e')
    f = Var.new('f')

    ff = ((a + d + b) * -c) % (d >> e >> f)
    print(ff)

    ff = a * (b * (c * f * (e * d)))
    print(ff)

    ff = a >> b >> c
    print(ff)

    ff = a % b
    print(ff)


if __name__ == '__main__':
    main()

# @formatter:on
