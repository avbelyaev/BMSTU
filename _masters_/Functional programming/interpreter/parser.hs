module Parser where

import Text.Read
import Debug.Trace
import System.IO.Unsafe
import Data.Unique

import Lexer


data Env = Env
    { tokens    :: [Token]
    , stack     :: [Int]    -- стек вычислений
    , nodes     :: [Node]   -- стек вершин, где вершина = терминал/нетерминал
    } deriving (Show)


-- ================================================
-- =================  Stack ops  ==================
-- ================================================
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

-- убрать верхний токен
nextToken e = e { tokens = tail $ tokens e }


-- ================================================
-- ================  AST builder  =================
-- ================================================
data Node = Node 
    { val       :: String
    , childs    :: [Node]
    , nodeId    :: Int    -- "UUID" ноды
    } deriving (Show, Eq)

-- deriving (Eq) для 2х вершин с одинаковыми val и childs должен выдавать разные хеши.
-- для этого добавлен nodeId -- хранящий хеш ноды
-- при этом вызывать надо для каждого объекта, иначе GHC инлайнит функцию
-- и она начинает возвращать одно и то же значение
makeNode val children = Node { val = val
                             , childs = children
                             , nodeId = hashUnique $ unsafePerformIO newUnique }

-- создание листовой вершины
makeLeafNode val = makeNode val []

-- создание вершины с 2 дочерними вершинами
makeBinaryNode val e = e { nodes = newNode : modNodes }
    where modNodes = tail $ tail $ nodes e
          n1 = nodes e !! 1
          n2 = nodes e !! 0
          newNode = makeNode val [n1, n2]


-- -//- с 3 вершинами
makeTernaryNode val extraNodeVal e = e { nodes = newNode : modNodes }
    where modNodes = tail $ tail $ nodes e
          n1 = nodes e !! 1
          n2 = nodes e !! 0
          nExtra = makeLeafNode extraNodeVal
          newNode = makeNode val [nExtra, n1, n2]


-- создание вершины с удинственной дочерней вершиной
makeNodeWithChild val childVal e = e { nodes = newNode : (nodes e) }
    where child = makeLeafNode childVal
          newNode = makeNode val [child]


-- создание вершины для правила F ::= ( <E> )
makeNodeFBracket e = e { nodes = newNode : modNodes}
    where modNodes = tail $ nodes e
          nE = nodes e !! 0
          nLeftParen = makeLeafNode "("
          nRightParen = makeLeafNode ")"
          newNode = makeNode "F" [nLeftParen, nE, nRightParen]


-- ================================================
-- ===============  LL(1) Parser   ================
-- ================================================
-- <E>  ::= <T> <E’>
-- <T>  ::= <F> <T’>
-- <E’> ::= + <T> <E’> | - <T> <E’> | e
-- <T’> ::= * <F> <T’> | / <F> <T’> | e
-- <F>  ::= <number> | <var> | ( <E> ) | - <F>


-- <E> ::= <T> <E’>
parse_E :: Env -> Env
parse_E e = 
    makeBinaryNode "E" (parse_E1 $ parse_T e)

-- <T> ::= <F> <T’>.
parse_T :: Env -> Env
parse_T e =
    makeBinaryNode "T" (parse_T1 $ parse_F e)

-- <E’> ::= + <T> <E’> | - <T> <E’> | e
parse_E1 :: Env -> Env
parse_E1 e = 
    let token = head $ tokens e
        op = v $ token
    in case op of
        "+" -> makeTernaryNode "E1" "+" (parse_E1 $ addOp $ parse_T $ nextToken e)
        "-" -> makeTernaryNode "E1" "-" (parse_E1 $ subOp $ parse_T $ nextToken e)
        _   -> makeNodeWithChild "E1" "e" e

-- <T’> ::= * <F> <T’> | / <F> <T’> | e
parse_T1 :: Env -> Env
parse_T1 e = 
    let token = head $ tokens e
        op = v $ token
    in case op of
        "*" -> makeTernaryNode "T1" "*" (parse_T1 $ mulOp $ parse_F $ nextToken e)
        "/" -> makeTernaryNode "T1" "/" (parse_T1 $ divOp $ parse_F $ nextToken e)
        _   -> makeNodeWithChild "T1" "e" e

-- <F> ::= <number> | ( <E> ) | - <F>
-- -<F> not implemented :)
parse_F :: Env -> Env
parse_F e =
    let token = head $ tokens e
        tokType = t $ token
        value = v $ token
    in case tokType of
        NumberToken     -> makeNodeWithChild "F" value (nextToken $ pushOp e token)
        BracketToken    -> makeNodeFBracket $ nextToken $ parse_E $ nextToken e
        _               -> error "unexpected token"



parse :: [Token] -> Env
parse tokenList = parse_E e
    where e = Env { tokens = tokenList, stack = [], nodes = [] }



-- ================================================
-- =============   Pretty printer   ===============
-- ================================================
printDigraph :: [Node] -> String
printDigraph ast = 
    "digraph { dummy [label = \"\", shape = none]\n" ++
        (printNodeList ast) ++ 
        (printEdgeList ast) ++ 
    "}"


dfs :: [Node] -> [Node] -> [Node]
dfs visited nodes =
    if null $ nodes
    then reverse visited
    else let curr = head nodes
             adjacent = childs curr
         in if curr `elem` visited
            then dfs visited (tail nodes)
            else dfs (curr:visited) (nodes ++ adjacent)


-- shoud produce line: "56 [label = "E1" shape = circle id = "56"]"
nodeToDigraphVertex node =
    "  " ++ show(nodeId node) ++ 
    " [label = \"" ++ (val node) ++ "\"" ++
    " shape = circle" ++
    " id = \"" ++ show(nodeId node) ++ "\"]"

-- печатает список вершин для AST
printNodeList :: [Node] -> String
printNodeList ast = unlines names
    where names = map nodeToDigraphVertex dfsNodes
          dfsNodes = dfs [] ast 


-- given: node{ v = 21 }, node{ v = 13 }
-- should produce line: "21 -> 13"
nodeToDigraphEdge :: Node -> Node -> String
nodeToDigraphEdge nodeFrom nodeTo =
    "  " ++ show(nodeId nodeFrom) ++ " -> " ++ show(nodeId nodeTo)

-- печатает все ребра для данной вершины
-- given: node{ v = 21, childNodes = [13, 18, 20] }
-- should produce line: "21 -> 13, 21 -> 18, 21 -> 20"
getOutgoingEdges :: Node -> String
getOutgoingEdges node = unwords edges
    where edges = edgeHelper node (childs node)
          edgeHelper node childNodes =
                if null childNodes
                then []
                else let edge = nodeToDigraphEdge node (head childNodes)
                     in [edge] ++ edgeHelper node (tail childNodes)

-- печатает лист ребер для данного AST
printEdgeList :: [Node] -> String
printEdgeList ast = unlines filteredEmptyLines
    where filteredEmptyLines = filter (\line -> not $ null line) edges
          edges = map getOutgoingEdges dfsNodes
          dfsNodes = dfs [] ast


-- run Parser after Lexer: 
-- - make sure Lexer's filename is lexer.hs,
-- - comment Lexer's `main`
-- - do `runhaskell parser.hs`
main = do
    putStrLn ">>> scanning (line should end with space!)"
    let p2 = "1 + 2 "
    -- let p2 = "2 * 5 + 3 * 4 + 4 * 6 "
    -- let p2 = "1 + ((2 + 3) * 1337 * (4 - 5)) * 4 - (1488 / 10 + 5) + 999 "
    let tokens = scan p2
    mapM_ print tokens

    putStrLn "\n>>> parsing"
    let res = parse tokens

    putStrLn "\n>>> result stack"
    print $ stack res

    putStrLn "\n>>> AST"
    let ast = nodes res
    print ast

    putStrLn "\n>>> Graphviz (go webgraphviz.com)"
    putStrLn $ printDigraph ast


    -- print ">>> DFS"
    -- let vis = dfs [] ast
    -- mapM_ print vis

    -- putStrLn "\n>>> edges to root vertex childs:"
    -- let root = head ast
    -- let edges = getOutgoingEdges root
    -- print edges


    
