

data Env = Env
	{ ws :: [String] 		-- program
	, wc :: Int 			-- word counter
	, xs :: [Int]			-- data stack. осн вычисления здесь
	, rs :: [Int]			-- return stack. для процедур
	, as :: [(String, Int)]	-- dictionary (процедура, номер слова)
	-------------------------- минимум
	, buf :: String			-- output buffer
	, ls :: [Int]			-- ?
	, vs :: [(String, Int)]	-- ?
	} deriving (Show, Read)
