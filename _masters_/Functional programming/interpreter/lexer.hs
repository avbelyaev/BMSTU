module Interp where

--import Text.Read
import Text.Regex.Posix

data Token = Token
    { v :: String
    -- , len :: Int
    -- , offs :: Int
    } deriving (Show, Read)


isInteger s = case reads s :: [(Integer, String)] of
    [(_, "")] -> True
    _         -> False


-- from=1 to=3 [a, b, c, d, e] -> [b, c, d]
slice :: Int -> Int -> [a] -> [a]
slice from to xs = take (to - from + 1) (drop from xs)


-- item=c [a, b, c, d] -> 2
-- item=c [a, b, d] -> err
indexOfItem :: String -> [String] -> Int
indexOfItem item xs =
    if null xs 
    then error ("could not find " ++ show item ++ " among " ++ show xs)
    else if item == (head xs)
         then 0
         else 1 + indexOfItem item (tail xs)


function :: String -> Token
function arg
  | 0 == fst (res arg "define")       = Token { v = "KEYWORD" }
  | otherwise                         = Token { v = "NONE" }
  where res :: String -> String -> (Int,Int)
        res text re = text =~ re :: (MatchOffset,MatchLength)


function2 :: String -> [Int]
function2 arg
  | (x,y) <- arg =~ "define" :: (MatchOffset,MatchLength) = [x, y]
  | otherwise = [123, 456]


main = do
    print "start interpreting"

    let t7 = "                                      \
             \ define =0? dup 0 = end               \
             \ define gcd                           \
             \       =0? if drop exit endif         \
             \       swap over mod                  \
             \       gcd                            \
             \ end                                  \
             \ 90 99 gcd                            \
             \ 234 8100 gcd                         "

    -- let tokens = lexer t7 pattern
    -- print tokens
    print $ function "define"

