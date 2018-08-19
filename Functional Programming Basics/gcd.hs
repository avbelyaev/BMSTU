mygcd :: Int -> Int -> Int
mygcd n m
 | m == 0 = n
 | otherwise = mygcd m (mod n m)
 
mylcm :: Int -> Int -> Int
mylcm n m = div (n * m) (mygcd n m)

