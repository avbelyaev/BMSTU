module Dz where
import Debug.Trace


-- 1.2

data ExactInterval = ExactInterval { from, to :: Float
						 		   } deriving Show

data AvgInterval = AvgInterval { avg :: Float
					 		   , delta :: Int
					 		   } deriving Show


infixl 6 +++
(+++) :: ExactInterval -> ExactInterval -> ExactInterval
a +++ b = ExactInterval { from = from a + from b, to = to a + to b }

subInts :: ExactInterval -> ExactInterval -> ExactInterval
subInts a b = ExactInterval (from a - to b) (to a - from b)

mulInts :: ExactInterval -> ExactInterval -> ExactInterval
mulInts a b = ExactInterval (minimum values) (maximum values)
	where values = [from a * from b, from a * to b, to a * from b, to a * to b]

divInts :: ExactInterval -> ExactInterval -> ExactInterval
divInts a b = ExactInterval (minimum values) (maximum values)
	where values = [from a / from b, from a / to b, to a / from b, to a / to b]

containsZero :: ExactInterval -> Bool
containsZero a = from a * to a < 0


avgToExact :: AvgInterval -> ExactInterval
avgToExact a = ExactInterval (avg a - deltaPart) (avg a + deltaPart)
	where deltaPart = avg a / 100 * (fromIntegral $ delta a)

exactToAvg :: ExactInterval -> AvgInterval
exactToAvg a = AvgInterval mid delta
  where
	deltaPart = (to a - from a) / 2
	mid = from a + deltaPart
	delta = truncate $ 100 * deltaPart / mid


main = do
	let r1 = AvgInterval {avg = 10, delta = 5}
	let r2 = AvgInterval {avg = 20, delta = 5}
	print r1
	print r2

	let r3 = ExactInterval {from = 20, to = 30}
	let r4 = ExactInterval {from = 30, to = 40}
	print $ r3 +++ r4


	let res1 = divInts (mulInts (avgToExact r1) (avgToExact r2)) 
					((avgToExact r1) +++ (avgToExact r2))
	print res1

	print $ containsZero $ ExactInterval {from = -3, to = -1}
	print $ containsZero $ ExactInterval {from = -3, to = 3}

