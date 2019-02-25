module Exc where

import Control.Exception


data IllegalArgumentException = IllegalArgumentException 
    { message :: String
    } deriving (Show)


instance Exception IllegalArgumentException


sublist :: [a] -> Int -> IO [a]
sublist list amount 
    | (length list) >= amount   = return $ take amount list
    | otherwise                 = throw (IllegalArgumentException 
                                        "amount cannot be grater than list size")


main = do
    print "calling normal"
    print =<< sublist [1, 2, 3, 4] 3

    print "catching exception"
    --let handler = (\_ -> print "Error") :: IllegalArgumentException -> IO ()
    --catch (sublist [1, 2, 3] 4) handler
    result <- try (sublist [1, 2, 3] 4) :: IO (Either IllegalArgumentException [Integer])
    -- :: IO (Either IllegalArgumentException [Integer])
    -- print result

    print "not catching"
    print =<< sublist [1, 2, 1,2,1,2,1,2] 3
  