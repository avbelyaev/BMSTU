public class SetElement<T> {
    private T x;
private int depth;
private SetElement<T> parent, rootx, rooty;

public SetElement (T x) {
	this.x = x;
	depth = 0;
	parent = this;
}

public T x() {
	return x;
}

private SetElement<T> Find (SetElement<T> elem) {
	if (elem.parent == this) {
		return elem.parent;
	} else {
		elem.parent = elem.parent.Find(elem.parent);
		return elem.parent;
	}
}

public Boolean equivalent(SetElement<T> elem) {
	rootx = this.Find(this);
	rooty = elem.Find(elem);
	if (rootx == rooty) {
		return true;
	}
	return false;
}

public void union (SetElement<T> elem) {
	rootx = this.Find(this);	//union for "a" in "a.union(b)"
	rooty = elem.Find(elem);	//union for "b" in "a.union(b)"
	if (rootx.depth < rooty.depth) {
		rootx.parent = rooty;
	} else {
		rooty.parent = rootx;
		if ((rootx.depth == rooty.depth) && (rootx != rooty)) {
			rootx.depth++;
		}
	}
        return;
}
}

//========================================================================================

public class Test 
{ 
        public static void main(String[] args) 
        { 
                SetElement<Integer> a = new SetElement<Integer> (0), 
                        b = new SetElement<Integer> (1), 
                        c = new SetElement<Integer> (2), 
                        d = new SetElement<Integer> (3), 
                        e = new SetElement<Integer> (4), 
                        f = new SetElement<Integer> (5); 
                a.union(b); 
                c.union(a); 
                c.union(d); 
                e.union(f); 
                System.out.println("" + a.x() + "=" + d.x() + ":␣" + a.equivalent(d)); 
                System.out.println("" + a.x() + "=" + f.x() + ":␣" + a.equivalent(f)); 
        } 
}

