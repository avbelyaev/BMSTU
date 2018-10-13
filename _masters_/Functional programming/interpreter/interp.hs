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



-- stack grows to the left. that means a,b in (a - b) represented as [b a]
doOp :: (Int -> Int -> Int) -> Env -> Env
doOp operation e = e { stack = res : modStack }
    where res = operation (stack e !! 1) (stack e !! 0)
          modStack = tail $ tail (stack e)

negOp e = e { stack = modStack }
    where modStack = -1 * head (stack e) : tail (stack e)


-- stack ops
dropOp e = e { stack = tail (stack e) }

swapOp e = e { stack = x : (y : restOfStack) }
    where x = stack e !! 0
          y = stack e !! 1
          restOfStack = tail $ tail (stack e)

dupOp e = e { stack = head (stack e) : (stack e) }





defineCustomFunc e = 
    let funcName = ws e !! (counter e + 1) -- func name goes after 'define'
        funcBodyStartIndex = counter e + 2 -- skip 'define', skip <name>
        modFuncs = funcs e ++ [(funcName, funcBodyStartIndex)] -- add (name, index) to known funcs

        restOfProgram = slice (counter e) (length $ ws e) (ws e) -- remove everything to the left
        funcBodyEndIndex = indexOfItem "end" restOfProgram
        modCounter = counter e + funcBodyEndIndex
    in trace ("func: " ++ show funcName ++ 
                "[" ++ show funcBodyStartIndex ++ ":" ++ show funcBodyEndIndex ++ "]")
        $ e { counter = modCounter, funcs = modFuncs }


callCustomFunc e =
    let funcNameToCall = ws e !! (counter e)
        funcNames = map (\fn -> fst fn) (funcs e) -- get names of all known functions
        funcToCall = funcs e !! (indexOfItem funcNameToCall funcNames) -- find function by name
        -- get funcAddr from (funcName, funcAddr). sub 1 since next expr in 'interp' is counter++
        funcCallAddress = snd funcToCall - 1 
        
        returnAddr = counter e + 1 -- memorize address where we should return after invocation
    in e { counter = funcCallAddress, retAddrs = returnAddr : retAddrs e }


endCustomFunc e = e { counter = addrToReturn, retAddrs = modRetAddrs }
    -- get retAddr from top of retAddrs stack. sub 1 since next expr in 'interp' is counter++
    where addrToReturn = head (retAddrs e) - 1 
          !modRetAddrs = tail (retAddrs e) -- remove it



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

                                -- stack ops
                                "drop"  -> dropOp e
                                "swap"  -> swapOp e
                                "dup"   -> dupOp e

                                "def"   -> defineCustomFunc e
                                "end"   -> endCustomFunc e
                                _       -> callCustomFunc e
                    in interp modEnv { counter = counter modEnv + 1}



eva :: String -> [Int] -> Env
eva program initStack = interp env
    where --ops = [ ("+", addOp)
          --      , ("-", subOp)
           --     , ("*", mulOp)
            --    , ("")] 
        env = Env { ws = words program, counter = 0, stack = initStack, retAddrs = [], funcs = [] }
    

assert prog expected = (expected == actualRes) && (null actualRetStack)
    where envRes = eva prog []
          actualRes = stack envRes
          actualRetStack = retAddrs envRes

main = do
    print "start interpreting"

    let t1 = assert "1 1 + 1 - 4 * 2 / neg" [-2]
    print t1

    print "============="

    let t2 = assert "def inc 1 + end 2 inc" [3]
    print t2


    --let res = stack envres !! 0
    --print res
    --return res

