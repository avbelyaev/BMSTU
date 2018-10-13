module Interp where

import Text.Read
import Debug.Trace


data Env = Env
    { ws :: [String]        -- program
    , counter :: Int            -- word counter
    , stack :: [Int]            -- data stack. осн вычисления здесь
    , rs :: [Int]           -- return stack. для процедур
    , funcs :: [(String, Int)]         -- dictionary (процедура, номер слова)
    } deriving (Show, Read)


isInteger s = case reads s :: [(Integer, String)] of
    [(_, "")] -> True
    _         -> False

slice from to xs = take (to - from + 1) (drop from xs)


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






{-
defineFunc env = 
    let funcName = wc env !! (counter env + 1) -- func name goes after 'define'
        funcBodyStartIndex = counter env + 2 -- skip 'define', skip <name>
        funcBodyEndIndex = findNearestItemIndex "end" (wc env) (counter env)
        modFuncs = funcs env ++ (funcName, funcBodyStartIndex)
        modCounter = counter env + funcBodyEndIndex
    in env { counter = modCounter, funcs = modFuncs }
-}
--userDefinedOp opName env = env { stack = modStack, counter = modCounter, rs = modRetStack }
 


interp :: Env -> Env
interp e = 
    let cs = "ws:" ++ show (ws e) ++
             " len=" ++ show (length (ws e)) ++
             " exit=" ++ show (counter e == length (ws e)) ++
             " counter:" ++ show (counter e) ++
             " stack:" ++ show (stack e) 
    in 
        if (null $ ws e) || (counter e == length (ws e))
        then e
        else let !current = trace cs $ ws e !! counter e -- force this fucker to evaluate
            in if isInteger current
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

                                --"def"   -> defineFunc e
                                --_       -> userDefinedOp e
                    in interp modEnv { counter = counter e + 1}



eva :: String -> [Int] -> Env
eva program initStack = interp env
    where --ops = [ ("+", addOp)
          --      , ("-", subOp)
           --     , ("*", mulOp)
            --    , ("")] 
        env = Env { ws = words program, counter = 0, stack = initStack, rs = [], funcs = [] }
    

assert prog expected = expected == actualRes
    where envRes = eva prog []
          actualRes = stack envRes

main = do
    print "start interpreting"

    let t1 = assert "1 1 + 1 - 4 * 2 / neg" [-2]
    print t1

    let t2 = assert "def inc 1 + exit" [1]
    print t2


    --let res = stack envres !! 0
    --print res
    --return res

