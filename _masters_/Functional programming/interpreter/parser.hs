module Parser where
-- module's name should be capitalized filename

import Text.Read
import Debug.Trace
import Data.Maybe

import Lexer


data Env = Env
    { tokens :: [Token]
    , stack  :: [Int]
    } deriving (Show)


isInteger s = case reads s :: [(Integer, String)] of
    [(_, "")] -> True
    _         -> False


toInt :: Token -> Maybe Int
toInt tok = case (isInteger $ v tok) of
    True    -> Just (read $ v tok :: Int)
    False   -> Nothing


binOp :: (Int -> Int -> Int) -> Env -> Env
binOp operation e = e { stack = binOpResult : modStack }
    where n1 = stack e !! 0
          n2 = stack e !! 1 
          binOpResult = operation n2 n1 
          modStack = tail $ tail (stack e)

addOp e = binOp (+) e
subOp e = binOp (-) e
mulOp e = binOp (*) e
divOp e = binOp (div) e

pushOp e val = e { stack = val : (stack e) }


nextToken e = modEnv
    where nextTok = head $ tokens e
          modEnv = e { tokens = tail $ tokens e }

-- <E> ::= <T> <E’>.
-- <T> ::= <F> <T’>.
-- <E’> ::= + <T> <E’> | - <T> <E’> | .
-- <T’> ::= * <F> <T’> | / <F> <T’> | .
-- <F> ::= <number> | <var> | ( <E> ) | - <F> .

parse :: [Token] -> Env
parse tokenList = parse_E e
    where e = Env { tokens = tokenList, stack = [] }


-- <E> ::= <T> <E’>.
parse_E :: Env -> Env
parse_E e = 
    parse_E1 $ trace " . E" (parse_T e)


-- <T> ::= <F> <T’>.
parse_T :: Env -> Env
parse_T e =
    parse_T1 $ trace " . . T" (parse_F e)


-- <E’> ::= + <T> <E’> | - <T> <E’> | .
parse_E1 :: Env -> Env
parse_E1 e = 
    let token = head $ tokens e
        op = v $ token
    in case op of
        "+" -> parse_E1 $ addOp $ parse_T $ trace " . . . E1 +" (nextToken e)
        "-" -> parse_E1 $ subOp $ parse_T $ trace " . . . E1 -" (nextToken e)
        _   -> trace " . . . E1 e" e


-- <T’> ::= * <F> <T’> | / <F> <T’> | .
parse_T1 :: Env -> Env
parse_T1 e = 
    let token = head $ tokens e
        op = v $ token
    in case op of
        "*" -> parse_T1 $ mulOp $ parse_F $ trace " . . . T1 *" (nextToken e)
        "/" -> parse_T1 $ divOp $ parse_F $ trace " . . . T1 /" (nextToken e)
        _   -> trace " . . . T1 e" e


-- <F> ::= <number> | <var> | ( <E> ) | - <F>.
parse_F :: Env -> Env
parse_F e =
    let token = head $ tokens e
        tokType = t $ token
        maybeIntToken = toInt token
    in case tokType of
        NumberToken     -> nextToken $ trace (" . . . F  " ++ v token) (pushOp e (fromJust maybeIntToken))
        BracketToken    -> nextToken $ parse_E $ nextToken e
        _               -> error "fuck"



-- run Parser after Lexer: 
-- - make sure Lexer's filename is lexer.hs,
-- - comment Lexer's `main`
-- - do `runhaskell parser.hs`
main = do
    --let p2 = "(5 - 1) * (2 + 1) "
    -- let p2 = "2 * 5 + 3 * 4 + 4 * 6 "
    let p2 = "1 + 2 + 3 - (4 * 5) + 6 "
    let tokens = scan p2
    mapM_ print tokens

    print "start parsing"
    let res = parse tokens
    print $ stack res
    
