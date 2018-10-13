module Interp where

import Text.Read
import Debug.Trace


data Env = Env
    { ws :: [String]            -- program
    , counter :: Int            -- pointer to current word
    , stack :: [Int]            -- data stack
    , retAddrs :: [Int]         -- return stack
    , funcs :: [(String, Int)]  -- dictionary (funcName, funcAddr)
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



-- basic binary ops
binOp :: (Int -> Int -> Int) -> Env -> Env
binOp operation e = e { stack = binOpResult : modStack, counter = counter e + 1 }
    where n1 = stack e !! 0 -- stack grows to the left. stack[0] is a top
          n2 = stack e !! 1 
          binOpResult = operation n2 n1 
          modStack = tail $ tail (stack e)

negOp e = e { stack = modStack, counter = counter e + 1 }
    where modStack = -1 * head (stack e) : tail (stack e)



-- stack ops
dropOp e = e { stack = tail (stack e), counter = counter e + 1 }

swapOp e = e { stack = n2 : (n1 : restOfStack), counter = counter e + 1 }
    where n1 = stack e !! 0
          n2 = stack e !! 1 
          restOfStack = tail $ tail (stack e) -- remove 2 values from top of stack

dupOp e = e { stack = head (stack e) : (stack e), counter = counter e + 1 }

overOp e = e { stack = modStack, counter = counter e + 1 }
    where secondFromLast = head $ tail $ stack e
          modStack = secondFromLast : (stack e)



-- comparison
customTrue = -1 -- why, forth? just why?
customFalse = 0

boolOp operation e = e { stack = comparisonResult : restOfStack, counter = counter e + 1 }
    where n2 = stack e !! 0 -- n2 is at the top of stack
          n1 = stack e !! 1
          comparisonResult = if operation n1 n2 then customTrue else customFalse
          restOfStack = tail $ tail (stack e) -- remove 2 values from top of stack



-- control-flow operations
-- define abs {...} end
defineCustomFunc e = 
    let funcName = ws e !! (counter e + 1) -- func name goes after 'define'
        funcBodyStartIndex = counter e + 2 -- skip 'define', skip <name>
        modFuncs = funcs e ++ [(funcName, funcBodyStartIndex)] -- add (name, index) to known already funcs

        restOfProgram = slice (counter e) (length $ ws e) (ws e) -- remove everything to the left
        funcBodyEndIndex = indexOfItem "end" restOfProgram -- index of function's 'end' from current pos
        modCounter = counter e + funcBodyEndIndex + 1 -- ctr = current + funcBody 
    in trace ("func: " ++ show funcName ++ 
                "[" ++ show funcBodyStartIndex ++ ".." ++ show funcBodyEndIndex ++ "]")
        $ e { counter = modCounter, funcs = modFuncs }


callCustomFunc e =
    let funcNameToCall = ws e !! (counter e)
        funcNames = map (\fn -> fst fn) (funcs e) -- get names of all known functions
        funcToCall = funcs e !! (indexOfItem funcNameToCall funcNames) -- find function by name
        -- get funcAddr from (funcName, funcAddr)
        funcCallAddress = snd funcToCall
        
        returnAddr = counter e + 1 -- memorize address where we should return after invocation
    in e { counter = funcCallAddress, retAddrs = returnAddr : retAddrs e }


endCustomFunc e = e { counter = addrToReturn, retAddrs = modRetAddrs }
    -- get retAddr from top of retAddrs stack
    where addrToReturn = head (retAddrs e)
          !modRetAddrs = tail (retAddrs e) -- remove it


-- determine, where to place counter - to 'if-block' or to 'else-block'
startBranch e = 
    let goToElseBlock = customFalse == (head $ stack e)
    in 
        if goToElseBlock
        then let restOfProgram = slice (counter e) (length $ ws e) (ws e) -- remove everything to the left
                 indexOfNearestEndif = indexOfItem "endif" restOfProgram
             in e { counter = counter e + indexOfNearestEndif, stack = tail (stack e) }
        else e { counter = counter e + 1, stack = tail (stack e) }


skip e = e { counter = counter e + 1 }




interp :: Env -> Env
interp e = 
    let cs = --"ws:" ++ show (ws e) ++
             " len=" ++ show (length (ws e)) ++
             " ctr:" ++ show (counter e) ++
             " funcs:" ++ show (funcs e) ++
             " ret:" ++ show (retAddrs e) ++
             " stack:" ++ show (stack e)
        progIsEmptyOrHasEnded = (null $ ws e) || (counter e == length (ws e))
    in 
        if progIsEmptyOrHasEnded
        then e
        else let !current = ws e !! counter e -- DEBUG trace cs $ ws e !! counter e
            in if isInteger current
                then let !newElem = (read current :: Int)       -- force this fucker to evaluate
                    in interp e { stack = newElem : (stack e), counter = counter e + 1 }
                else let !modEnv = case current of
                                -- basic ops
                                "+"         -> binOp (+) e
                                "-"         -> binOp (-) e
                                "*"         -> binOp (*) e
                                "/"         -> binOp (div) e
                                "mod"       -> binOp (mod) e
                                "neg"       -> negOp e

                                -- comparison
                                "="         -> boolOp (==) e
                                ">"         -> boolOp (>) e
                                "<"         -> boolOp (<) e

                                -- stack ops
                                "drop"      -> dropOp e
                                "swap"      -> swapOp e
                                "dup"       -> dupOp e
                                "over"      -> overOp e

                                -- control-flow objects
                                "define"    -> defineCustomFunc e
                                "end"       -> endCustomFunc e
                                "exit"      -> endCustomFunc e
                                "if"        -> startBranch e
                                "endif"     -> skip e
                                _           -> callCustomFunc e
                    in interp modEnv


-- initialize interpreter and evalate program
eva :: String -> [Int] -> Env
eva program initialStack = interp env
    where env = Env { ws = words program
                    , counter = 0
                    , stack = initialStack
                    , retAddrs = []
                    , funcs = [] }
    

assert prog expected = (expected == actualRes) && retStackIsEmpty
    where resultEnv = eva prog []
          actualRes = stack resultEnv
          retStackIsEmpty = null (retAddrs resultEnv)



main = do
    print "start interpreting"

    let t1 = assert "1 1 + 1 - 4 * 2 / neg 744 neg *" [1488]
    print t1

    print "-----"
    let t2 = assert " define inc            \
                    \       1 1 + *         \
                    \ end                   \
                    \ 3 inc                 " [6]
    print t2

    print "-----"
    let t3 = assert " define abs            \
                    \      dup 0 <          \
                    \      if neg endif     \
                    \ end                   \
                    \  9 abs                \
                    \ -9 abs                " [9, 9]
    print t3

    print "-----"
    let t4 = assert " define =0? dup 0 = end            \
                    \ define <0? dup 0 < end            \
                    \ define signum                     \
                    \       =0? if exit endif           \
                    \       <0? if drop -1 exit endif   \
                    \       drop 1                      \
                    \ end                               \
                    \  0 signum                         \
                    \ -5 signum                         \
                    \ 10 signum                         " [1, -1, 0]
    print t4


    print "-----"
    let t5 = assert "  define -- 1 - end                \
                     \ define =0? dup 0 = end           \
                     \ define =1? dup 1 = end           \
                     \ define factorial                 \
                     \      =0? if drop 1 exit endif    \
                     \      =1? if drop 1 exit endif    \
                     \      dup --                      \
                     \      factorial *                 \
                     \ end                              \
                     \ 0 factorial                      \
                     \ 1 factorial                      \
                     \ 2 factorial                      \
                     \ 3 factorial                      \
                     \ 4 factorial                      " [24, 6, 2, 1, 1]
    print t5


    print "-----"
    let t6 = assert " define =0? dup 0 = end                \
                    \ define =1? dup 1 = end                \
                    \ define -- 1 - end                     \
                    \ define fib                            \
                    \       =0? if drop 0 exit endif        \
                    \           =1? if drop 1 exit endif    \
                    \           -- dup                      \
                    \           -- fib                      \
                    \           swap fib +                  \
                    \ end                                   \
                    \ define make-fib                       \
                    \       dup 0 <                         \
                    \       if drop exit endif              \
                    \           dup fib                     \
                    \       swap --                         \
                    \       make-fib                        \
                    \ end                                   \
                    \ 10 make-fib                           " [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
    print t6


    print "-----"
    let t7 = assert " define =0? dup 0 = end        \
                    \ define gcd                    \
                    \       =0? if drop exit endif  \
                    \       swap over mod           \
                    \       gcd                     \
                    \ end                           \
                    \ 90 99 gcd                     \
                    \ 234 8100 gcd                  " [18, 9]
    print t7

