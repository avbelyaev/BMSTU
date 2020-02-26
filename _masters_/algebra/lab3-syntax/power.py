def power(x, n):
    res = 1
    while n != 0:
        if n % 2 == 0:
            print "A"
            x, n = x*x, n/2
        else:
            print "B"
            res, n = res * x, n-1
    return res
 
print power(3, 2)