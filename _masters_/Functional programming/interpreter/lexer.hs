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

isLetterOrDigit c = isLetter c || isDigit c


isEOP :: String -> (Bool, String)
isEOP text = (matches, tail text)
    where matches = '\0' == head text 


isWhitespace text = (matches, tail text)
    where matches = head text `elem` [' ', '\n', '\r', '\t']


isKeyword text = checkPrefix "" text
    where isKey word = word `elem` ["define", "end"]

          checkPrefix prefix str =
            if null str
            then (isKey prefix, str)

            else if isLetterOrDigit $ head str
                then checkPrefix (prefix ++ [head str]) (tail str)
                else (isKey prefix, str)


recognizeIdent :: String -> String -> (String, String, Bool)
recognizeIdent ident text 
    | (isDigit curr || isLetter curr) 
        && not (isLetterOrDigit next)  = (ident ++ [curr], tail text, False)
    | (isDigit curr || isLetter curr)           = recognizeIdent (ident ++ [curr]) (tail text)
    | otherwise                                 = ("ident-err", text, True)
    where curr = head text 
          next = head $ tail text
    

isIdent text = 
    if isLetter $ head text
    then let (ident, follow, error) = recognizeIdent "" text
         in case error of
                True    -> (False, text)
                False   -> (True, follow)
    else (False, text)


recognizeDigit :: String -> String -> (String, String, Bool)
recognizeDigit digit text 
    | (isDigit curr) 
        && not (isLetterOrDigit next)  = (digit ++ [curr], tail text, False)
    | isDigit curr                              = recognizeDigit (digit ++ [curr]) (tail text)
    | otherwise                                 = ("digit-err", text, True)
    where curr = head text 
          next = head $ tail text


isNumeric text = 
    if isDigit $ head text
    then let (digit, follow, error) = recognizeDigit "" text
         in case error of 
                True    -> (False, text)
                False   -> (True, follow)
    else (False, text)


scan :: String -> [Token]
scan text
    | (True, _)    <- isEOP text            = Token { v = "EOP" } : []
    | (True, follow) <- isWhitespace text   = scan follow
    | (True, follow) <- isKeyword text      = Token { v = "KEY" } : scan follow  
    | (True, follow) <- isIdent text        = Token { v = "IDT" } : scan follow
    | (True, follow) <- isNumeric text      = Token { v = "NUM" } : scan follow
    | otherwise                             = error "Unexpected character!" 
    

{-| isDigit $ head symbols        = Token { v = "DIGIT", len = len } : (lexx $ tail symbols)

    
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

    let prog = " define gcd2 dend \n\r\t 123 dfdg  end \0"
    print $ scan prog
    print "done"




