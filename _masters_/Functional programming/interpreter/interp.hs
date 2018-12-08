module Interp where

import Text.Read
import Debug.Trace


data Env = Env
    { ws :: [String]            -- program
    , counter :: Int            -- pointer to current word
    , stack :: [Int]            -- data stack
    , callStack :: [Int]         -- return stack
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


-- 'define abc 1 neg if exit'. counter points to '1' -> 'neg if exit'
cutOffLeftPart :: Env -> [String]
cutOffLeftPart e = slice (counter e) (length $ ws e) (ws e)




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
          restOfStack = tail $ tail (stack e)                 -- remove 2 values from top of stack

dupOp e = e { stack = head (stack e) : (stack e), counter = counter e + 1 }

overOp e = e { stack = modStack, counter = counter e + 1 }
    where secondFromLast = head $ tail $ stack e
          modStack = secondFromLast : (stack e)



-- comparison
customTrue = -1 -- why, forth? just why?
customFalse = 0

boolOp operation e = e { stack = comparisonResult : restOfStack, counter = counter e + 1 }
    where n2 = stack e !! 0                                     -- n2 is at the top of stack
          n1 = stack e !! 1
          comparisonResult = if operation n1 n2 then customTrue else customFalse
          restOfStack = tail $ tail (stack e)                   -- remove 2 values from top of stack



-- control-flow operations
-- define abs {...} end
defineCustomFunc e = 
    let funcName = ws e !! (counter e + 1)                      -- func name goes after 'define'
        funcBodyStartIndex = counter e + 2                      -- skip 'define', skip <name>
        modFuncs = (funcName, funcBodyStartIndex) : funcs e     -- add (name, addr) to known already funcs

        restOfProgram = cutOffLeftPart e                        -- remove everything to the left
        funcBodyEndIndex = indexOfItem "end" restOfProgram      -- index of function's 'end' from current pos
        modCounter = counter e + funcBodyEndIndex + 1           -- ctr = current + funcBodyLength
    
    in trace ("f: " ++ show funcName ++ 
              " [" ++ show funcBodyStartIndex ++ ".." ++ show funcBodyEndIndex ++ "]")
        $ e { counter = modCounter, funcs = modFuncs }


-- lookup for function among known ones, and call if present
callCustomFunc e =
    let funcNameToCall = ws e !! (counter e)                            -- name of function to be called
        funcNames = map (\fn -> fst fn) (funcs e)                       -- get names of all known functions
        funcToCall = funcs e !! (indexOfItem funcNameToCall funcNames)  -- find function among function names
        funcAddress = snd funcToCall                                    -- get funcAddr from (funcName, funcAddr)
        
        returnAddr = counter e + 1                                      -- remember address where we should return 
    in e { counter = funcAddress, callStack = returnAddr : callStack e }


-- find function's 'end', revert callStack, and return to caller
endCustomFunc e = e { counter = addrToReturn, callStack = modcallStack }
    where addrToReturn = head (callStack e)                              -- get address from top of callStack stack
          modcallStack = tail (callStack e)                              -- remove it


-- determine, where to place counter - to 'if-block' or to 'else-block'
startBranch e = 
    let goToElseBlock = customFalse == (head $ stack e)                 -- should be go straight to 'Else'?
    in 
        if goToElseBlock
        then let restOfProgram = cutOffLeftPart e
                 indexOfNearestEndif = indexOfItem "endif" restOfProgram    -- find where is the nearest 'endif'
             in e { counter = counter e + indexOfNearestEndif, stack = tail (stack e) } -- remove flag from stack
        else e { counter = counter e + 1, stack = tail (stack e) }          -- remove flag from stack


-- just inc the counter to skip 'endif' word
skip e = e { counter = counter e + 1 }



-- takes initialized environment and evaluates prog word-by-word
interp :: Env -> Env
interp e = 
    if counter e == length (ws e)
    then e                                          
    else let current = ws e !! counter e
         in if isInteger current
            then interp e { stack = (read current :: Int) : (stack e), counter = counter e + 1 }
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
                            -- control-flow
                            "define"    -> defineCustomFunc e
                            "end"       -> endCustomFunc e
                            "exit"      -> endCustomFunc e
                            "if"        -> startBranch e
                            "endif"     -> skip e
                            _           -> callCustomFunc e
                in interp modEnv


-- initializes interpreter environment. returns interpreter's stack
eva :: String -> [Int] -> [Int]
eva program initialStack = stack evaluated 
    where initialEnv = Env { ws = words program
                           , counter = 0
                           , stack = initialStack
                           , callStack = []
                           , funcs = [] }
          evaluated = interp initialEnv


-- launch: $ runhaskell interp.hs
main = do
    print "start interpreting"

    let t1 = eva "1 1 + 1 - 4 * 2 / neg 744 neg *" [1337]
    print $ t1 == [1488, 1337]


    print "-----"
    let t2 = eva " define inc            \
                 \       1 1 + *         \
                 \ end                   \
                 \ 3 inc                 " []
    print $ t2 == [6]


    print "-----"
    let t3 = eva " define abs            \
                 \      dup 0 <          \
                 \      if neg endif     \
                 \ end                   \
                 \  9 abs                \
                 \ -9 abs                " []

    print $ t3 == [9, 9]


    print "-----"
    let t4 = eva " define =0? dup 0 = end               \
                 \ define <0? dup 0 < end               \
                 \ define signum                        \
                 \       =0? if exit endif              \
                 \       <0? if drop -1 exit endif      \
                 \       drop 1                         \
                 \ end                                  \
                 \  0 signum                            \
                 \ -5 signum                            \
                 \ 10 signum                            " []

    print $ t4 == [1, -1, 0]


    print "-----"
    let t5 = eva "  define -- 1 - end                   \
                 \ define =0? dup 0 = end               \
                 \ define =1? dup 1 = end               \
                 \ define factorial                     \
                 \      =0? if drop 1 exit endif        \
                 \      =1? if drop 1 exit endif        \
                 \      dup --                          \
                 \      factorial *                     \
                 \ end                                  \
                 \ 0 factorial                          \
                 \ 1 factorial                          \
                 \ 2 factorial                          \
                 \ 3 factorial                          \
                 \ 4 factorial                          " []

    print $ t5 == [24, 6, 2, 1, 1]


    print "-----"
    let t6 = eva " define =0? dup 0 = end                \
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
                 \ 10 make-fib                           " []

    print $ t6 == [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55]


    print "-----"
    let t7 = eva " define =0? dup 0 = end               \
                 \ define gcd                           \
                 \       =0? if drop exit endif         \
                 \       swap over mod                  \
                 \       gcd                            \
                 \ end                                  \
                 \ 90 99 gcd                            \
                 \ 234 8100 gcd                         " []

    print $ t7 == [18, 9]
    

