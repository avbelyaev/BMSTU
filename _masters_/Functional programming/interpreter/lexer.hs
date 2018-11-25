module Interp where

--import Text.Read
import Text.Regex.Posix
import Data.Char
import Debug.Trace


data Token = Token
    { v :: String
    --, len :: Int
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


{-
lexer :: String -> Token
lexer arg = case matcher of
     (_,len) <- arg =~ "de(fi)ne" :: (MatchOffset,MatchLength)
    (_,len):_         -> Token { v = "KEYWORD",   len = len }
    (offs,len):_      -> Token { v = "KEYWORD",   len = len }
    --_               -> Token { v = "NONE",      len = -1 }

    --where regexMatcher :: String -> String -> (Int,Int)
    --      regexMatcher text regex = text =~ regex :: (MatchOffset,MatchLength)


function2 :: String -> [Int]
function2 arg
  | (x,y) <- arg =~ "define" :: (MatchOffset,MatchLength) = [x, y]
  | otherwise = [123, 456]
-}


isEOP :: String -> (Bool, String)
isEOP follow = (matches, tail follow)
    where matches = '\0' == head follow 


isWhitespace follow = (matches, tail follow)
    where matches = head follow `elem` [' ', '\n', '\r', '\t']


isKeyword follow = checkPrefix "" follow
    where isKey word = word `elem` ["define", "end"]

          checkPrefix prefix str =
            if null str
            then (isKey prefix, str)

            else if isLetter $ head str
                then checkPrefix (prefix ++ [head str]) (tail str)
                else (isKey prefix, str)


recognizeIdent :: String -> String -> (String, String, Bool)
recognizeIdent ident follow 
    | (isDigit curr || isLetter curr) 
        && not (isDigit next || isLetter next)  = (ident ++ [curr], tail follow, False)

    | (isDigit curr || isLetter curr)           = recognizeIdent (ident ++ [curr]) (tail follow)
    | otherwise                                 = (ident, follow, True)

    where curr = head follow 
          next = head $ tail follow
    

isIdent follow = 
    if isLetter $ head follow
    then let (ident, rest, error) = recognizeIdent "" follow
         in case error of
                True   -> (False, follow)
                False  -> (True, rest)
    else (False, follow)





scan :: String -> [Token]
scan follow
    | (True, _)    <- isEOP follow          = Token { v = "EOP" } : []
    | (True, rest) <- isWhitespace follow   = scan rest
    | (True, rest) <- isKeyword follow      = Token { v = "KEY" } : scan rest  
    | (True, rest) <- isIdent follow        = Token { v = "IDT" } : scan rest
    
    | otherwise                             = Token { v = "-" } : (scan $ tail follow)
    

{-| isDigit $ head symbols        = Token { v = "DIGIT", len = len } : (lexx $ tail symbols)

    | (True, rest) <- isNumber follow       = Token { v = "NUM" } : scan rest
-}

main = do
    print "start scanning"

    let t7 = "                                      \
             \ define =0? dup 0 = end               \
             \ define gcd                           \
             \       =0? if drop exit endif         \
             \       swap over mod                  \
             \       gcd                            \
             \ end                                  \
             \ 90 99 gcd                            \
             \ 234 8100 gcd                         "

    --let idt = "id1337  <rest of programm>"
    --print $ recognizeIdent "" idt
    --print $ isIdent idt

    let prog = " define dend \n\r\t  dfdg   \0"
    print $ scan prog
    print "done"




