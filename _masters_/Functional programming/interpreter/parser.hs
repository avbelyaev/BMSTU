module Parser where
-- module's name should be capitalized filename

import Text.Read
import Debug.Trace

import Lexer


data Env = Env
    { tokens :: [Token]
    , stack  :: [Int]
    } deriving (Show)


toInt s = read s :: Int


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


parse :: [Token] -> Env
parse tokenList = parse_E e
    where e = Env { tokens = tokenList, stack = [] }


-- <E> ::= <T> <E’>.
parse_E :: Env -> Env
parse_E e = 
    trace " E -> T E1" (parse_E1 $ parse_T e)


-- <T> ::= <F> <T’>.
parse_T :: Env -> Env
parse_T e =
    trace " T -> F T1" (parse_T1 $ parse_F e)


-- <E’> ::= + <T> <E’> | - <T> <E’> | .
parse_E1 :: Env -> Env
parse_E1 e = 
    let token = head $ tokens e
        op = trace ("E1:" ++ show token) (v $ token)
    in case op of
        "+" -> parse_E1 $ addOp $ parse_T $ nextToken e
        "-" -> parse_E1 $ subOp $ parse_T $ nextToken e
        _   -> e


-- <T’> ::= * <F> <T’> | / <F> <T’> | .
parse_T1 :: Env -> Env
parse_T1 e = 
    let token = head $ tokens e
        op = trace ("T1:" ++ show token) (v $ token)
    in case op of
        "*" -> parse_T1 $ mulOp $ parse_F $ nextToken e
        "/" -> parse_T1 $ divOp $ parse_F $ nextToken e
        _   -> trace "T1 skip" e


-- <F> ::= <number> | <var> | ( <E> ) | - <F>.
parse_F :: Env -> Env
parse_F e =
    let token = head $ tokens e
        tokType = trace (" F:" ++ show token) (t $ token)
    in case tokType of
        NumberToken     -> nextToken $ pushOp e (toInt $ v token)
        BracketToken    -> nextToken $ parse_E $ nextToken e
        _               -> e



-- run Parser after Lexer: 
-- - make sure Lexer's filename is lexer.hs,
-- - comment Lexer's `main`
-- - do `runhaskell parser.hs`
main = do
    let p2 = "(5 - 1) * 2 + 1 "
    let tokens = scan p2
    mapM_ print tokens

    print "start parsing"
    let res = parse tokens
    print $ stack res
    
