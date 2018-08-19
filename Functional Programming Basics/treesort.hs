data Tree a = Leaf | Node (Tree a) a (Tree a)
 
--insert :: Ord a => a -> Tree a -> Tree a
 
insert x Leaf = Node Leaf x Leaf
 
insert x (Node t y t') | x <= y = Node (insert x t) y t'
insert x (Node t y t') | x > y = Node t y (insert x t')
 
--flatten :: Tree a -> [a]
 
flatten Leaf = []
flatten (Node t x t') = flatten t ++ [x] ++ flatten t'
 
 
 
--treesort :: Ord a => [a] -> [a]
 
treesort = flatten . foldr insert Leaf

