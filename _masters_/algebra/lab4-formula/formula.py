import copy
from abc import ABC
from typing import List, Any, Optional

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
    def newConjunction(args: 'List[Any]') -> 'Formula': return Formula.new(Connective.conjAnd(), args)

    @staticmethod
    def newDisjunction(args: 'List[Any]') -> 'Formula': return Formula.new(Connective.disjOr(), args)

    @staticmethod
    def conjAnd() -> 'Connective': return Connective(True)

    @staticmethod
    def disjOr() -> 'Connective': return Connective(False)

    def __repr__(self): return self.__str__()
    def __str__(self): return '*' if self.neutral else '+'


class Impl(Symbol):
    @staticmethod
    def newImplication(args: List[Any]) -> 'Formula': return Formula.new(Impl(), args)

    def __repr__(self): return self.__str__()
    def __str__(self): return '=>'


class Eq(Symbol):
    @staticmethod
    def newEquivalence(args: List[Any]) -> 'Formula': return Formula.new(Eq(), args)

    def __repr__(self): return self.__str__()
    def __str__(self): return '<=>'


class Neg(Symbol):
    @staticmethod
    def newNegative(arg: Any) -> 'Formula': return Formula.new(Neg(), [arg])

    def __repr__(self): return self.__str__()
    def __str__(self): return '!'


class Var(Symbol):
    def __init__(self, value: str):
        """ use factory methods instead! """
        self.value = value
        self.index = nextIndex()  # уникальный номер переменной

    @staticmethod
    def new(value: str) -> 'Var':
        return Var(value)

    def __repr__(self): return f'{self.value}'

    # пришлось перегружать все операторы для Var, чтобы 2 рядом стоящие (a+b) и (b+c)
    # состояли из 3х Var'ов, но обе 'b' были независмы
    def __add__(self, other):
        right = other
        if isinstance(other, Var):
            right = Formula.literal(other)
        left = Formula.literal(self)
        return left + right

    def __mul__(self, other):
        right = other
        if isinstance(other, Var):
            right = Formula.literal(other)
        left = Formula.literal(self)
        return left * right

    def __rshift__(self, other):
        right = other
        if isinstance(other, Var):
            right = Formula.literal(other)
        left = Formula.literal(self)
        return left >> right

    def __mod__(self, other):
        right = other
        if isinstance(other, Var):
            right = Formula.literal(other)
        left = Formula.literal(self)
        return left % right

    def __neg__(self):
        left = Formula.literal(self)
        return -left


class Formula:
    def __init__(self, sym: Optional[Symbol], args: 'List[Any]'):
        """ use factory methods instead! """
        self._symbol = sym
        self._args = args
        self.canonize()

    @staticmethod
    def new(sym: Symbol, args: 'List[Formula]') -> 'Formula':
        return Formula(sym, args)

    @staticmethod
    def literal(var: 'Var') -> 'Formula':
        return Formula(None, [var])

    def canonize(self):
        if isinstance(self._symbol, Eq):                     # инвариант-2: a <=> b === (a => b) AND (b => a)
            self._symbol = Connective.conjAnd()
            a = copy.deepcopy(self._args[0])
            b = copy.deepcopy(self._args[1])
            self._args[0] = Impl.newImplication([a, b])
            self._args[1] = Impl.newImplication([b, a])

        if isinstance(self._symbol, Impl):                   # инвариант-2: a => b === !a OR b
            self._symbol = Connective.disjOr()
            self._args[0] = -self._args[0]

    @property
    def isConjunction(self) -> bool: return isinstance(self._symbol, Connective) and self._symbol.neutral

    @property
    def isDisjunction(self) -> bool: return isinstance(self._symbol, Connective) and not self._symbol.neutral

    @property
    def isNegated(self) -> bool: return isinstance(self._symbol, Neg)

    @property
    def isLiteral(self) -> bool: return isinstance(self._args[0], Var)

    def addArguments(self, args: 'List[Formula]'):
        self._args.extend(args)
        self._args.sort(key=lambda arg: arg._args[0].value if arg.isLiteral else 'dummy')

    def __mul__(self, other):
        # a * (b * c) === *(a, b, c)
        if self.isLiteral and other.isConjunction:
            other.addArguments([self])
            return other

        # (a * b) * c === *(a, b, c)
        if self.isConjunction and other.isLiteral:
            self.addArguments([other])
            return self

        # (a * b) * (c * d) === *(a, b, c, d)
        if self.isConjunction and other.isConjunction:
            self.addArguments(other._args)
            return self

        return Connective.newConjunction([self, other])

    def __add__(self, other):
        # a + (b + c) === +(a, b, c)
        if self.isLiteral and other.isDisjunction:
            other.addArguments([self])
            return other

        # (a + b) + c === +(a, b, c)
        if self.isDisjunction and other.isLiteral:
            self.addArguments([other])
            return self

        # (a + b) + (c + d) === +(a, b, c, d)
        if self.isDisjunction and other.isDisjunction:
            self.addArguments(other._args)
            return self

        return Connective.newDisjunction([self, other])

    def __neg__(self):
        if self.isLiteral:
            # need to point a new formula to the same var
            if self.isNegated:
                self._symbol = None
            else:
                self._symbol = Neg()
            return self
        if self.isNegated:
            return self._args[0]
        return Neg.newNegative(self)

    def __rshift__(self, other): return Impl.newImplication([self, other])
    def __mod__(self, other): return Eq.newEquivalence([self, other])

    def __repr__(self): return self.__str__()

    def __str__(self):
        symbol = f'{self._symbol}' if self._symbol is not None else ''
        if self.isLiteral:
            return f'{symbol}{self._args[0]}'

        s = f'{symbol}('
        for arg in self._args:
            s += f'{arg} '
        s = s[0: len(s) - 1]
        s += ')'
        return s


class Normalizer:
    @staticmethod
    def toNegationNormalForm(formula: Formula, b: bool = False) -> Formula:
        f = copy.deepcopy(formula)

        if f.isLiteral:
            if b:
                f = -f
            return f

        # Отрицание переносится вниз по синтаксическому дереву
        if f.isNegated:
            f = -f
            b = not b

        # Применение законов де Моргана
        for i, arg in enumerate(f._args):
            argCopy = copy.deepcopy(arg)
            f._args[i] = Normalizer.toNegationNormalForm(argCopy, b)

        if f.isConjunction and b:
            f._symbol = Connective.disjOr()
        elif f.isConjunction and not b:
            f._symbol = Connective.conjAnd()
        elif f.isDisjunction and b:
            f._symbol = Connective.conjAnd()
        elif f.isDisjunction and not b:
            f._symbol = Connective.disjOr()
        return f


def main():
    a = Var.new('a')
    b = Var.new('b')
    c = Var.new('c')
    d = Var.new('d')
    e = Var.new('e')
    f = Var.new('f')
    x = Var.new('x')
    y = Var.new('y')
    z = Var.new('z')

    # уплощение
    testSort = b * (-a * (c * f * (e * d)))
    assert str(testSort) == '*(!a b c d e f)'
    print(testSort)

    # де морган
    testDeMorgan = -(a * b)
    assert str(testDeMorgan) == '!(*(a b))'
    testDeMorgan = Normalizer.toNegationNormalForm(testDeMorgan)
    print(testDeMorgan)

    # пример из википедии
    testNNF = -((a >> b) + -(b >> c))
    assert str(testNNF) == '!(+(+(!a b) !(+(!b c))))'
    testNNF = Normalizer.toNegationNormalForm(testNNF)
    print(testNNF)


if __name__ == '__main__':
    main()

# @formatter:on
