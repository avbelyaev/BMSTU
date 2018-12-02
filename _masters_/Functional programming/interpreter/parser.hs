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


parse :: Env -> Env
parse e = trace "S" (parse_E e)

-- <E> ::= <T> <E’>.
parse_E :: Env -> Env
parse_E e = 
    trace "E" (parse_E1 $ parse_T e)


-- <T> ::= <F> <T’>.
parse_T :: Env -> Env
parse_T e =
    trace "T" (parse_T1 $ parse_F e)


-- <E’> ::= + <T> <E’> | - <T> <E’> | .
parse_E1 :: Env -> Env
parse_E1 e = 
    let token = trace "E1" (head $ tokens e)
        op = trace (">>>E1:" ++ show token) (v $ token)
    in case op of
        "+" -> parse_E1 $ addOp $ parse_T $ nextToken e
        "-" -> nextToken $ parse_E1 $ subOp $ parse_T e
        _   -> e


-- <T’> ::= * <F> <T’> | / <F> <T’> | .
parse_T1 :: Env -> Env
parse_T1 e = 
    let token = trace "T1" (head $ tokens e)
        op = trace (">>>T1:" ++ show token) (v $ token)
    in case op of
        "*" -> nextToken $ parse_T1 $ mulOp $ parse_F e
        "/" -> nextToken $ parse_T1 $ divOp $ parse_F e
        _   -> e


-- <F> ::= <number> | <var> | ( <E> ) | - <F>.
parse_F :: Env -> Env
parse_F e =
    let token = trace "F" (head $ tokens e)
        tokType = trace (">>> F:" ++ show token) (t $ token)
    in case tokType of
        NumberToken -> nextToken $ pushOp e (toInt $ v token)
        _           -> e
 



-- run Parser after Lexer: 
-- - make sure Lexer's filename is lexer.hs,
-- - comment Lexer's `main`
-- - do `runhaskell parser.hs`
main = do
    let p2 = "1 + 2 + 51 "
    let tokens = scan p2
    mapM_ print tokens

    print "start parsing"
    let e = Env { tokens = tokens, stack = [] }
    let res = parse e
    print res
    
