isPrime :: Int -> Bool
 
isPrime num = primeTest num 2
    where
      primeTest :: Int -> Int -> Bool
      primeTest num x
            | x == num      = True
            | num `mod` x == 0 = False
            | otherwise    = primeTest num (x + 1)

