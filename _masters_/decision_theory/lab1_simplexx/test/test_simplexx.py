#! /usr/bin/env python3.6
import unittest

# relative to venv ? venv now at /BMSTU/venv
from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx


class TestSimplexxMethods(unittest.TestCase):
    def test_change_basis(self):
        simplex = Simplexx()
        print(simplex)

        # TODO unittest change_basis()


if __name__ == '__main__':
    unittest.main()
