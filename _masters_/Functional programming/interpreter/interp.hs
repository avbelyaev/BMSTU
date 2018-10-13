module Interp where

import Text.Read
import Debug.Trace


data Env = Env
    { ws :: [String]        -- program
    , counter :: Int            -- word counter
    , stack :: [Int]            -- data stack. осн вычисления здесь
    , retAddrs :: [Int]           -- return stack. для процедур
    , funcs :: [(String, Int)]         -- dictionary (процедура, номер слова)
    } deriving (Show, Read)


isInteger s = case reads s :: [(Integer, String)] of
    [(_, "")] -> True
    _         -> False


slice :: Int -> Int -> [a] -> [a]
slice from to xs = take (to - from + 1) (drop from xs)


indexOfItem :: String -> [String] -> Int
indexOfItem item xs =
    let eqq = trace ("eq(" ++ show item ++ "," ++ show (head xs) ++ "):" ++ show ((head xs) == item)) 
                $ item == (head xs)
    in 
        if null xs 
        then error ("could not find " ++ show item ++ " among " ++ show xs)
        else if eqq
             then 0
             else 1 + indexOfItem item (tail xs)



-- basic ops
doOp :: (Int -> Int -> Int) -> Env -> Env
doOp operation e = e { stack = res : modStack, counter = counter e + 1 }
    where res = operation (stack e !! 1) (stack e !! 0) -- stack grows to the left
          modStack = tail $ tail (stack e)

negOp e = e { stack = modStack, counter = counter e + 1 }
    where modStack = -1 * head (stack e) : tail (stack e)


-- stack ops
dropOp e = e { stack = tail (stack e), counter = counter e + 1 }

swapOp e = e { stack = n2 : (n1 : restOfStack), counter = counter e + 1 }
    where n2 = stack e !! 0 -- n2 is at the top of stack
          n1 = stack e !! 1
          restOfStack = tail $ tail (stack e) -- remove 2 values from top of stack

dupOp e = e { stack = head (stack e) : (stack e), counter = counter e + 1 }


-- comparison
customTrue = -1 -- why, forth? just why?
customFalse = 0

boolOp operation e = e { stack = comparisonResult : restOfStack, counter = counter e + 1 }
    where n2 = stack e !! 0 -- n2 is at the top of stack
          n1 = stack e !! 1
          comparisonResult = if operation n1 n2 then customTrue else customFalse
          restOfStack = tail $ tail (stack e) -- remove 2 values from top of stack



-- control-flow operations
defineCustomFunc e = 
    let funcName = ws e !! (counter e + 1) -- func name goes after 'define'
        funcBodyStartIndex = counter e + 2 -- skip 'define', skip <name>
        modFuncs = funcs e ++ [(funcName, funcBodyStartIndex)] -- add (name, index) to known funcs

        restOfProgram = slice (counter e) (length $ ws e) (ws e) -- remove everything to the left
        funcBodyEndIndex = indexOfItem "end" restOfProgram -- index of function's 'end' from current pos
        modCounter = counter e + funcBodyEndIndex + 1 -- ctr = current + funcBody 
    in trace ("func: " ++ show funcName ++ 
                "[" ++ show funcBodyStartIndex ++ ":" ++ show funcBodyEndIndex ++ "]")
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
    let goToElseBlock = trace ("goToElse=" ++ show (customFalse == (head $ stack e)))
            $ customFalse == (head $ stack e)
    in 
        if goToElseBlock
        then let restOfProgram = slice (counter e) (length $ ws e) (ws e) -- remove everything to the left
                 indexOfNearestEndif = indexOfItem "endif" restOfProgram
             in e { counter = counter e + indexOfNearestEndif, stack = tail (stack e) }
        else e { counter = counter e + 1, stack = tail (stack e) }


skip e = e { counter = counter e + 1 }




interp :: Env -> Env
interp e = 
    let cs = "ws:" ++ show (ws e) ++
             " len=" ++ show (length (ws e)) ++
             " ctr:" ++ show (counter e) ++
             " funcs:" ++ show (funcs e) ++
             " ret:" ++ show (retAddrs e) ++
             " stack:" ++ show (stack e) 
    in 
        if (null $ ws e) || (counter e == length (ws e))
        then e
        else let !current = trace cs $ ws e !! counter e -- force this fucker to evaluate
            in if isInteger (trace ("curr:" ++ show current) current)
                then let !newElem = (read current :: Int)
                    in interp e { stack = newElem : (stack e), counter = counter e + 1 }
                else let !modEnv = case current of
                                -- basic ops
                                "+"     -> doOp (+) e
                                "-"     -> doOp (-) e
                                "*"     -> doOp (*) e
                                "/"     -> doOp (div) e
                                "mod"   -> doOp (mod) e
                                "neg"   -> negOp e

                                -- comparison
                                "="     -> boolOp (==) e
                                ">"     -> boolOp (>) e
                                "<"     -> boolOp (<) e

                                -- stack ops
                                "drop"  -> dropOp e
                                "swap"  -> swapOp e
                                "dup"   -> dupOp e

                                -- control-flow objects
                                "def"   -> defineCustomFunc e
                                "end"   -> endCustomFunc e
                                "exit"  -> endCustomFunc e
                                "if"    -> startBranch e
                                "endif" -> skip e
                                _       -> callCustomFunc e
                    in interp modEnv



eva :: String -> [Int] -> Env
eva program initStack = interp env
    where --ops = [ ("+", addOp)
          --      , ("-", subOp)
           --     , ("*", mulOp)
            --    , ("")] 
        env = Env { ws = words program, counter = 0, stack = initStack, retAddrs = [], funcs = [] }
    

assert prog expected = (expected == actualRes) && retStackIsEmpty
    where resultEnv = eva prog []
          actualRes = stack resultEnv
          retStackIsEmpty = null (retAddrs resultEnv)

main = do
    print "start interpreting"

    let t1 = assert "1 1 + 1 - 4 * 2 / neg" [-2]
    print t1

    print "============="
    let t2 = assert "def inc 1 1 + * end 3 inc" [6]
    print t2

    print "============="
    let t3 = assert "def cmp 7 5 > if exit endif 20 end cmp" []
    print t3

    print "============="
    let t4 = assert "def -- 1 - end" []
    print t4
    --let res = stack envres !! 0
    --print res
    --return res

