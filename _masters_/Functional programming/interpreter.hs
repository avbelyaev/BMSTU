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


add :: [Int] ->Int
add stack = trace ("add" ++ show (stack)) 
                $ head stack + head (tail stack)


interp :: Env -> Env
interp env = 
    let cs = "curr:" ++ show (ws env !! counter env) ++ 
             " ws:" ++ show (ws env) ++
             " len=" ++ show (length (ws env)) ++
             " exit=" ++ show (counter env == length (ws env) - 1) ++
             " counter:" ++ show (counter env) ++
             " stack:" ++ show (stack env) 
        !current = trace cs $ ws env !! counter env -- force this fucker to evaluate
    in 
        if (null $ ws env) || (counter env == length (ws env) - 1)
        then env
        else if isInteger current
            then let !newElem = (read current :: Int)
                in trace "int" $ interp env { stack = newElem : (stack env)
                              , counter = counter env + 1 }
            else let !res = case current of
                            "+" -> add $ stack env
                            _ -> 228
                in trace "str" $ interp env { stack = res : (tail (tail (stack env)))
                              , counter = counter env + 1 } 



eva :: String -> [Int] -> Env
eva program initStack = interp env
    where env = Env { ws = words program, counter = 0, stack = initStack, rs = [], as = Map.empty }
    


main = do
    print "start interpreting"

    let x = (read "2" :: Int)
    print x

    let prog = ": inc 1 + ; 1 inc ."
    let p2 = "1 2 3 + 4 +"
    let envres = eva p2 []
    print $ stack envres

    let res = stack envres !! 0
    print res

