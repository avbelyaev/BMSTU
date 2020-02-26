from typing import List
from functools import partial

print = partial(print, end=' ')


# @formatter:off
def prefix(arr: List[int], len_arr: int) -> List[int]:
    pi = [0] * len_arr
    i = 1
    while i < len_arr:                      # condLoopOverArray
        print('A')
        t = pi[i - 1]                       # statementA
        while t > 0 and arr[t] != arr[i]:   # condLoopOverPrefix
            print('B')
            t = pi[t - 1]                   # statementB

        if arr[t] == arr[i]:                # condIf
            print('C')
            t += 1                          # statementC

        print('D')
        pi[i], i = t, i + 1                 # statementD
    return pi


if __name__ == '__main__':
    print("\nModel 1\n")
    print(prefix([5, 6], 2))

    print("\nModel 2\n")
    print(prefix([7, 7], 2))

    print("\nModel 3\n")
    print(prefix([3, 4, 3], 3))

    print("\nModel 4\n")
    print(prefix([1, 2, 3, 1], 4))
# @formatter:on
