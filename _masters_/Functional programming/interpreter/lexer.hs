module Lexer where
-- module's name should be capitalized filename

import Data.Char
import Data.Maybe


data TokType = EopToken
             | KeywordToken
             | BracketToken
             | OperatorToken
             | IdentToken
             | NumberToken
             deriving (Show, Eq, Enum)

data Token = Token
    { t :: TokType
    , v :: String   -- token value
    --, len :: Int
    --, offs :: Int
    } deriving (Show, Eq)


isLetterOrDigit c = isLetter c || isDigit c

isWhitespace c = isSeparator c || isControl c


-- (распознанное значение в случае успеха, остаток программы справа)
type RecognitionResult = (Maybe String, String)

recognizeIdent :: String -> String -> RecognitionResult
recognizeIdent ident text 
    | (isLetterOrDigit curr) 
        && not (isLetterOrDigit next)   = (Just $ ident ++ [curr], tail text)
    | (isDigit curr || isLetter curr)   = recognizeIdent (ident ++ [curr]) (tail text)
    | otherwise                         = (Nothing, text)
    where curr = head text 
          next = head $ tail text


recognizeDigit :: String -> String -> RecognitionResult
recognizeDigit digit text 
    | (isDigit curr) 
        && not (isLetterOrDigit next)   = (Just $ digit ++ [curr], tail text)
    | isDigit curr                      = recognizeDigit (digit ++ [curr]) (tail text)
    | otherwise                         = (Nothing, text)
    where curr = head text 
          next = head $ tail text


-- checks if prefix of given 'string' matches given 'predicate'
-- accumulator holds prefix of 'string': [] [fuck], [f] [uck], [fu] [ck], ...
recognizePrefix :: String -> String -> (String -> Bool) -> RecognitionResult
recognizePrefix accumulator string predicate =
    if not $ isWhitespace $ head string
    then recognizePrefix (accumulator ++ [head string]) (tail string) predicate
    else if predicate accumulator
         then (Just accumulator, string)
         else (Nothing, string)



-- (заматчился ли токен, значение токена если заматчился, оставшаяся справа часть программы)
type MatchResult = (Bool, Maybe String, String)


isEOP :: String -> MatchResult
isEOP text = (matches, Just "EOP", tail text)
    where matches = null text || '\0' == head text 


isWhSpace :: String -> MatchResult
isWhSpace text = (matches, Nothing, tail text) -- ws is not needed and so ignored
    where matches = isWhitespace $ head text


isKeyword :: String -> MatchResult
isKeyword text = 
    let isKeywordPredicate word = word `elem` ["define", "end", "var", "val", "return"]
        (maybeKeyword, follow)  = recognizePrefix "" text isKeywordPredicate
    in case maybeKeyword of
            Just keyword    -> (True, Just keyword, follow)
            Nothing         -> (False, Nothing, text)
          

isBracket :: String -> MatchResult
isBracket text = case (head text) of
    '(' -> (True, Just "(", tail text)
    ')' -> (True, Just ")", tail text)
    _   -> (False, Nothing, text)


isIdent :: String -> MatchResult
isIdent text = 
    if isLetter $ head text
    then let (maybeIdent, follow) = recognizeIdent "" text
         in case maybeIdent of
                Just ident  -> (True, Just ident, follow)
                Nothing     -> (False, Nothing, text)
    else (False, Nothing, text)


isNumeric :: String -> MatchResult
isNumeric text = 
    if isDigit $ head text
    then let (maybeDigit, follow) = recognizeDigit "" text
         in case maybeDigit of 
                Just digit  -> (True, Just digit, follow)
                Nothing     -> (False, Nothing, text)
    else (False, Nothing, text)


isOpertor :: String -> MatchResult
isOpertor text = 
    let isOperatorPredicate word = word `elem` ["+", "-", "*", "/", "="]
        (maybeOerator, follow) = recognizePrefix "" text isOperatorPredicate
    in case maybeOerator of 
        Just operator       -> (True, Just operator, follow)
        Nothing             -> (False, Nothing, text)



scan :: String -> [Token]
scan text
    | (True, v, _)      <- isEOP text       = Token { t = EopToken, v = fromJust v } : []
    | (True, _, follow) <- isWhSpace text   = scan follow
    | (True, v, follow) <- isKeyword text   = Token { t = KeywordToken,  v = fromJust v } : scan follow 
    | (True, v, follow) <- isBracket text   = Token { t = BracketToken,  v = fromJust v } : scan follow
    | (True, v, follow) <- isOpertor text   = Token { t = OperatorToken, v = fromJust v } : scan follow
    | (True, v, follow) <- isIdent text     = Token { t = IdentToken,    v = fromJust v } : scan follow
    | (True, v, follow) <- isNumeric text   = Token { t = NumberToken,   v = fromJust v } : scan follow
    | otherwise                             = error "Unexpected character!" 
    

-- run Lexer only: uncomment following and do `runhaskell lexer.hs`
{-
main = do
    print "start scanning"

    let t1 = " define sum (a b)         \
             \      var result = a + b  \
             \      return result       \
             \      end                 \
             \ sum(3 5)                 "

    let tokens = scan t1
    mapM_ print tokens
-}
