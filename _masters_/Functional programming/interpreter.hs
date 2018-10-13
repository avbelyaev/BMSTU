module Interp where

import qualified Data.Map as Map
import Text.Read
import Debug.Trace


type FuncMap = Map.Map String Int

data Env = Env
    { ws :: [String]        -- program
    , counter :: Int            -- word counter
    , stack :: [Int]            -- data stack. осн вычисления здесь
    , rs :: [Int]           -- return stack. для процедур
    , as :: FuncMap         -- dictionary (процедура, номер слова)
    } deriving (Show, Read)


isInteger s = 
    let res = case reads s :: [(Integer, String)] of
                [(_, "")] -> True
                _         -> False
    in trace ("isInt=" ++ show res) $ res


-- stack grows to the left. that means a,b in (a - b) represented as [b a]
doOp :: (Int -> Int -> Int) -> [Int] -> [Int]
doOp op stack = res : modStack
    where res = op (stack !! 1) (stack !! 0)
          modStack = tail $ tail stack

addOp env = doOp (+) (stack env)
negOp stack = -1 * head stack : tail stack

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
                    in trace "int" $ interp e { stack = newElem : (stack e)
                                                , counter = counter e + 1 }
                else let !modStack = case current of
                                "+"     -> addOp e
                                "-"     -> doOp (-) (stack e) 
                                "*"     -> doOp (*) (stack e)
                                "/"     -> doOp (div) (stack e)
                                "mod"   -> doOp (mod) (stack e)
                                "neg"   -> negOp $ stack e
                    in trace "str" $ interp e { stack = modStack
                                                , counter = counter e + 1 } 



eva :: String -> [Int] -> Env
eva program initStack = interp env
    where env = Env { ws = words program, counter = 0, stack = initStack, rs = [], as = Map.empty }
    


main = do
    print "start interpreting"

    let prog = ": inc 1 + ; 1 inc ."
    let p2 = "1 1 + 1 - 4 * 2 / neg"
    let envres = eva p2 []
    print $ stack envres

    let res = stack envres !! 0
    print res

