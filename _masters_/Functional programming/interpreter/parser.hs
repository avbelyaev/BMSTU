module Parser where

import Text.Read
import Debug.Trace
import Data.Maybe

import Lexer


data Node = Node 
    { val       :: String
    , childs    :: [Node]
    } deriving (Show)

data Env = Env
    { tokens    :: [Token]
    , stack     :: [Int]
    , nodes     :: [Node]
    } deriving (Show)


toInt :: String -> Int
toInt val = read $ val :: Int


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

pushOp e token = e { stack = intVal : (stack e) }
    where intVal = toInt $ v token

nextToken e = modEnv
    where nextTok = head $ tokens e
          modEnv = e { tokens = tail $ tokens e }


-- pretty printer funcs

makeGenericNode val = Node { val = val, childs = [] }

makeEmptyNode = makeGenericNode "e"

makeNodeE e = e { nodes = newNode : modNodes }
    where modNodes = tail $ tail $ nodes e
          nodeT  = nodes e !! 0
          nodeE1 = nodes e !! 1
          newNode = Node { val = "E", childs = [nodeT, nodeE1]}

makeNodeT e = e { nodes = newNode : modNodes }
    where modNodes = tail $ tail $ nodes e
          nodeF  = nodes e !! 0
          nodeT1 = nodes e !! 1
          newNode = Node { val = "T", childs = [nodeF, nodeT1]}

makeNodeE1 e = e { nodes = newNode : modNodes }
    where modNodes = tail $ tail $ nodes e
          nodeT  = nodes e !! 0
          nodeE1 = nodes e !! 1
          nodeOp = makeGenericNode "(+/-)"
          newNode = Node { val = "E1", childs = [nodeOp, nodeT, nodeE1] }

makeNodeE1Empty e = e { nodes = newNode : (nodes e) }
    where newNode = Node { val = "E1", childs = [makeEmptyNode] }

makeNodeT1 e = e { nodes = newNode : modNodes }
    where modNodes = tail $ tail $ nodes e
          nodeF  = nodes e !! 0
          nodeT1 = nodes e !! 1
          nodeOp = makeGenericNode "(*/:)"
          newNode = Node { val = "T1", childs = [nodeOp, nodeF, nodeT1] }

makeNodeT1Empty e = e { nodes = newNode : (nodes e) }
    where newNode = Node { val = "T1", childs = [makeEmptyNode] }

makeNodeF token e = e { nodes = newNode : (nodes e) }
    where nodeVal = makeGenericNode $ v token
          newNode = Node { val = "F", childs = [nodeVal] }

makeNodeFBracket e = e { nodes = newNode : modNodes }
    where modNodes = tail $ nodes e
          nodeE = nodes e !! 0
          nodeLeftParen = makeGenericNode "("
          nodeRightParen = makeGenericNode ")"
          newNode = Node { val = "F", childs = [nodeLeftParen, nodeE, nodeRightParen] }

-- <E> ::= <T> <E’>
-- <T> ::= <F> <T’>
-- <E’> ::= + <T> <E’> | - <T> <E’> | e
-- <T’> ::= * <F> <T’> | / <F> <T’> | e
-- <F> ::= <number> | <var> | ( <E> ) | - <F>


-- <E> ::= <T> <E’>
parse_E :: Env -> Env
parse_E e = 
    makeNodeE $ parse_E1 $ parse_T e

-- <T> ::= <F> <T’>.
parse_T :: Env -> Env
parse_T e =
    makeNodeT $ parse_T1 $ parse_F e

-- <E’> ::= + <T> <E’> | - <T> <E’> | e
parse_E1 :: Env -> Env
parse_E1 e = 
    let token = head $ tokens e
        op = v $ token
    in case op of
        "+" -> makeNodeE1 $ parse_E1 $ addOp $ parse_T $ nextToken e
        "-" -> makeNodeE1 $ parse_E1 $ subOp $ parse_T $ nextToken e
        _   -> makeNodeE1Empty e

-- <T’> ::= * <F> <T’> | / <F> <T’> | e
parse_T1 :: Env -> Env
parse_T1 e = 
    let token = head $ tokens e
        op = v $ token
    in case op of
        "*" -> makeNodeT1 $ parse_T1 $ mulOp $ parse_F $ nextToken e
        "/" -> makeNodeT1 $ parse_T1 $ divOp $ parse_F $ nextToken e
        _   -> makeNodeT1Empty e

-- <F> ::= <number> | <var> | ( <E> ) | - <F>
parse_F :: Env -> Env
parse_F e =
    let token = head $ tokens e
        tokType = t $ token
    in case tokType of
        NumberToken     -> makeNodeF token (nextToken $ pushOp e token)
        BracketToken    -> makeNodeFBracket $ nextToken $ parse_E $ nextToken e
        _               -> error "unexpected token"



parse :: [Token] -> Env
parse tokenList = parse_E e
    where e = Env { tokens = tokenList
                  , stack = []
                  , nodes = [] }


-- run Parser after Lexer: 
-- - make sure Lexer's filename is lexer.hs,
-- - comment Lexer's `main`
-- - do `runhaskell parser.hs`
main = do
    --let p2 = "(5 - 1) * (2 + 1) "
    -- let p2 = "2 * 5 + 3 * 4 + 4 * 6 "
    let p2 = "1 + 2 "
    let tokens = scan p2
    mapM_ print tokens

    print "start parsing"
    let res = parse tokens
    print $ stack res

    print $ nodes res
    

    
