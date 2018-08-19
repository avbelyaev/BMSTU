//encode-decode e.g: encode [1,2,3,2,3,2]  -> [(1,1),(2,3),(3,2)]; decode [(2,'a'),(3,'b'),(1,'c')] -> "aabbbc"
import Data.List
 
encode xs = map func (group xs) where func xs = (length xs, last xs)
       
decode xs = concat (map (func) xs) where func (y, x) = take y (repeat x)

