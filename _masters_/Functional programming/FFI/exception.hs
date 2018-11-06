module Exc where

import Control.Exception


data IllegalArgumentException = IllegalArgumentException 
    { message :: String
    } deriving (Show)


instance Exception IllegalArgumentException


sublist :: [a] -> Int -> [a]
sublist list amount 
    | (length list) >= amount   = take amount list
    | otherwise                 = throw (IllegalArgumentException "fuck u")


main = do
    print "calling normal"
    print $ sublist [1, 2, 3, 4] 3

    print "expecting exception"
    print $ sublist [1, 2, 3] 4
    print "ok"

