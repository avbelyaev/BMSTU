module Interp where

--import Text.Read
import Text.Regex.Posix
import Data.Char

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
isEOP symbols = (matches, rest)
    where matches = '\0' == head symbols
          rest = tail symbols


isWhitespace symbols = (matches, rest)
    where matches = head symbols `elem` [' ', '\n', '\r', '\t']
          rest = tail symbols


isKeyword symbols = checkPrefix "" symbols
    where isKey word = word `elem` ["define", "end"]

          checkPrefix prefix str =
            if null str
            then (isKey prefix, str)

            else if isLetter $ head str
                then checkPrefix (prefix ++ [head str]) (tail str)
                else (isKey prefix, str)


recognizeIdent :: String -> String -> (String, String, Bool)
recognizeIdent ident symbols 
    | (isDigit curr || isLetter curr) 
        && not (isDigit next || isLetter next)  = (ident ++ [curr], tail symbols, False)

    | (isDigit curr || isLetter curr)           = recognizeIdent (ident ++ [curr]) (tail symbols)
    | otherwise                                 = (ident, symbols, True)

    where curr = head symbols 
          next = head $ tail symbols
    

isIdent symbols = 
    let curr = head symbols
    in  if (isLetter curr || isDigit curr)
        then let (ident, follow, error) = recognizeIdent "" symbols
                in (True, symbols) -- TODO remove line
             --in case (follow, error) of
               --      (word, True) -> (False, follow)
                 --    (word, False) -> (True, follow)
        else (False, symbols)


nextToken :: String -> [Token]
nextToken symbols
    | (True, _) <- isEOP symbols            = [ Token { v = "EOP" } ]
    | (True, rest) <- isWhitespace symbols  = nextToken rest
    | (True, rest) <- isKeyword symbols     = Token { v = "KEY" } : nextToken rest
    | (True, rest) <- isIdent symbols       = Token { v = "IDT" } : nextToken rest
    | otherwise                             = Token { v = "-" } : (nextToken $ tail symbols)
    

{-| isDigit $ head symbols        = Token { v = "DIGIT", len = len } : (lexx $ tail symbols)
    | isVar $ head symbols          = Token { v = "VAR", len = len } : (lexx $ tail symbols)
    | isKeyword $ head symbols      = Token { v = "KEYWORD", len = len } : (lexx $ tail symbols)
    | isEOF $ head symbols          = Token { v = "EOF", len = len } : (lexx $ tail symbols)
    | otherwise                     = error("undefined symbol")
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

    let idt = "id1337  <rest of programm>"
    print $ recognizeIdent "" idt
    --let prog = " defineend \n\r\t  dfdg   \0"
    --print $ nextToken prog




