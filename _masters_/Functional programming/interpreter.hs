module Interp where

import qualified Data.Map as Map
import Text.Read


type FuncMap = Map.Map String Int

data Env = Env
	{ ws :: [String] 		-- program
	, counter :: Int 			-- word counter
	, stack :: [Int]			-- data stack. осн вычисления здесь
	, rs :: [Int]			-- return stack. для процедур
	, as :: FuncMap			-- dictionary (процедура, номер слова)
	} deriving (Show, Read)


isInteger s = case reads s :: [(Integer, String)] of
  [(_, "")] -> True
  _         -> False



interp :: Env -> [Int]
interp env = 
	let current = ws env !! counter env
		--remains = tail $ ws env
		--progLen = length $ ws env
	in 
		if (null $ ws env) || (counter env == length (ws env) - 1)
		then stack env
		else if isInteger current
			then let newElem = (read current :: Int)
				in interp env { stack = newElem : (stack env)
							  , counter = counter env + 1 }
			else [228]




eva :: String -> [Int] -> [Int] 
eva program initStack = interp env
	where env = 
		Env { ws = words program, counter = 0, stack = initStack, rs = [], as = Map.empty }
	


main = do
	print "start interpreting"

	let x = (read "2" :: Int)
	print x

	let prog = ": inc 1 + ; 1 inc ."
	let p2 = "1 2 +"
	let res = eva p2 []
	print res

