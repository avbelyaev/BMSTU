class Elem(num: Int, mod: Int) {
    val x = num % mod
    val m = mod

    def + (k: Int) = new Elem(x + k, m)

    def + (q: Elem): Elem = {
        if (m != q.m)
            throw new Exception("modules not equal")
        else
            new Elem(x + q.x, m)
    }

    def * (k: Int) = new Elem(x * k, m)

    def * (q: Elem): Elem = {
        if (m != q.m) 
            throw new Exception("modules not equal")
        else 
            new Elem(x * q.x, m)
    }

    def this(mod: Int) = this(1, mod)
}


class ElemFactor(x: Int) {
    def + (q: Elem) = q + x
    def * (q: Elem) = q * x
}

implicit def intToFactor(i: Int) = new ElemFactor(i)

