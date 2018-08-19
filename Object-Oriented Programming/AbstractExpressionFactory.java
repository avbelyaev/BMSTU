public interface Expr<T> 
{ 
        void accept(ExprVisitor<T> v); 
}

public interface Neg<T> extends Expr<T> 
{ 
        Expr<T> a(); // операнд унарного минуса 
}

public interface Binary<T> extends Expr<T> 
{ 
        Expr<T> a(); // первый операнд 
        Expr<T> b(); // второй операнд 
        char operation(); // ’+’, ’-’, ’*’, ’/’ 
}

public interface Var<T> extends Expr<T> 
{ 
}

public interface Const<T> extends Expr<T> 
{ 
        T value(); // значение константы 
}

//========================================================================================

public interface ExprVisitor<T> 
{ 
        void visitNeg(Neg<T> e); 
        void visitBinary(Binary<T> e); 
        void visitVar(Var<T> e); 
        void visitConst(Const<T> e); 
}

//========================================================================================

public interface AbstractExprFactory<T> 
{ 
        Neg<T> newNeg(Expr<T> a); 
        Binary<T> newBinary(Expr<T> a, Expr<T> b, char operation); 
        Var<T> newVar(); 
        Const<T> newConst(T value); 
}

//========================================================================================

public class ExprFactory<T> implements AbstractExprFactory<T> {
        private class Neg1<T> implements Neg<T> {
		private Expr<T> a1;
		public Neg1 (Expr<T> a) {
			this.a1 = a;
		}
		public Expr<T> a() {
			return a1;
		}
		public void accept (ExprVisitor<T> v) {
			v.visitNeg(this);
		}
	}
	
	private class Binary1<T> implements Binary<T> {
		private Expr<T> a1;
		private Expr<T> b1;
		private char operation1;
		public Binary1 (Expr<T> a, Expr<T> b, char operation) {
			this.a1 = a;
			this.b1 = b;
			this.operation1 = operation;
		}
		public Expr<T> a() {
			return a1;
		}
		public Expr<T> b() {
			return b1;
		}
		public char operation() {
			return operation1;
		}
		public void accept (ExprVisitor<T> v) {
			v.visitBinary(this);
		}
	}
	
	private class Var1<T> implements Var<T> {
		public void accept (ExprVisitor<T> v) {
			v.visitVar(this);
		}
	}
	
	private class Const1<T> implements Const<T> {
		private T value1;
		public Const1 (T value) {
			this.value1 = value;
		}
		public  T value () {
			return value1;
		}
		public void accept (ExprVisitor<T> v) {
			v.visitConst(this);
		}
	}
//==================================================
	public Neg<T> newNeg (Expr<T> a) {
		return new Neg1<T>(a);
	}

	public Binary<T> newBinary (Expr<T> a, Expr<T> b, char operation) {
		return new Binary1<T>(a, b, operation);
	}
	
	public Var<T> newVar () {
		return new Var1<T>();
	}

	public  Const<T> newConst (T value) {
		return new Const1<T>(value);
	}
}

//========================================================================================

// CountVisitor предназначен для подсчёта количества 
// узлов разного типа в дереве арифметических выражений. 
class CountVisitor<T> implements ExprVisitor<T> 
{ 
        private int negs, binaries, vars, consts; 
 
        public CountVisitor(Expr<T> e) 
        { 
                e.accept(this); 
        } 
 
        public int negs() { return negs; } 
        public int binaries() { return binaries; } 
        public int vars() { return vars; } 
        public int consts() { return consts; } 
 
        public void visitNeg(Neg<T> e) 
        { 
                negs++; 
                e.a().accept(this); 
        } 
 
        public void visitBinary(Binary<T> e) 
        { 
                binaries++; 
                e.a().accept(this); 
                e.b().accept(this); 
        } 
 
        public void visitVar(Var<T> e) { vars++; } 
 
        public void visitConst(Const<T> e) { consts++; } 
} 
 
public class Test 
{ 
        public static void main(String[] args) 
        { 
                ExprFactory<Integer> f = new ExprFactory<Integer>(); 
 
                Expr<Integer> e = 
                        f.newBinary( 
                                f.newBinary( 
                                        f.newVar(), 
                                        f.newConst(10), 
                                        ’+’ 
                                ), 
                                f.newVar(), 
                                ’*’ 
                        ); 
 
                CountVisitor<Integer> v = new CountVisitor<Integer>(e); 
                System.out.println(v.negs()); 
                System.out.println(v.binaries()); 
                System.out.println(v.vars()); 
                System.out.println(v.consts()); 
        } 
}

