public interface Expr<T> {
	void accept(ExprVisitor<T> v);
}

public interface Neg<T> extends Expr<T> {
	Expr<T> a(); 	// операнд унарного минуса
}

public interface Binary<T> extends Expr<T> {
	Expr<T> a(); 	// первый операнд 
    	Expr<T> b(); 	// второй операнд 
    	char operation(); // ’+’, ’-’, ’*’, ’/’
}	

public interface Var<T> extends Expr<T> {
}

public interface Const<T> extends Expr<T> {
	T value(); 	// значение константы 
}

//========================================================================================
//========================================================================================
//========================================================================================

public interface AbstractExprFactory<T> {
	Neg<T> newNeg(Expr<T> a); 
    	Binary<T> newBinary(Expr<T> a, Expr<T> b, char operation); 
   	 Var<T> newVar(); 
    	Const<T> newConst(T value); 
}

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
	//==============================
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
//========================================================================================
//========================================================================================

public interface ExprVisitor<T> {
	void visitNeg(Neg<T> e); 
   	 void visitBinary(Binary<T> e); 
   	 void visitVar(Var<T> e); 
   	 void visitConst(Const<T> e);
}

public class EvalVisitor implements ExprVisitor<Integer>{
	private int res;
	private int x1;
	
	public EvalVisitor(int x) {
		this.x1 = x;
	}
	
	public int eval (Expr<Integer> e) {
		e.accept(this);		//recognize subexpression
		return res;			//update result of it on level above
	}
	
	public void visitNeg(Neg<Integer> e) {
		res = -eval(e.a());
	}
	
	public void visitBinary(Binary<Integer> e) {
		if (e.operation() == '+') {
			res = eval(e.a()) + eval(e.b());//dig into operand 
		} else if (e.operation() == '-') {  //to check if its subexpr
				   res = eval(e.a()) - eval(e.b());
			   } else if (e.operation() == '*') {
					      res = eval(e.a()) * eval(e.b());
					  } else if (e.operation() == '/') {
						         res = eval(e.a()) / eval(e.b());
					         }
	}
	
	public void visitVar(Var<Integer> e) {
		res = x1;	//cant dig deeper => get and return variable
	}
	
	public void visitConst(Const<Integer> e) {
		res = e.value();	// -||- => get and return constant
	}
}

//========================================================================================
//========================================================================================
//========================================================================================

public interface StackMachine {
	// Положить константу в стек 
    	void pushConst(int value); 

    // Положить в стек значение переменной x 
    	void pushVar(); 

    // Изменить знак числа на вершине стека 
    	void neg(); 

    // Бинарные арифметические операции. 
    // Каждая операция снимает два операнда со стека 
    // и кладёт на стек результат 
    	void add();  	// сложение 
    	void sub();  	// вычитание 
    	void mul();  	// умножение 
    	void div();  	// деление 
}

public class ExprBuilder implements StackMachine{
        private AbstractExprFactory<Integer> factory1;
	private Expr a, b;
	private Expr[] data;
	private int top;
	
	public ExprBuilder (AbstractExprFactory<Integer> factory) {
		data = new Expr[100];
		this.factory1 = factory;
		top = 0;
	}

	public Expr<Integer> pop () {
		return data[--top];
	}
	
	public void pushConst (int value) {
		data[top++] = factory1.newConst(value);
	}
	
	public void pushVar() {
		data[top++] = factory1.newVar();
	}
 	
	public void neg() {
		this.a = pop();
		data[top++] = factory1.newNeg(a);
	}
	
	public void add() {
		this.a = pop();
		this.b = pop();
		data[top++] = factory1.newBinary(a, b, '+');
	}
	
	public void sub() {
		this.a = pop();
		this.b = pop();
		data[top++] = factory1.newBinary(b, a, '-');
	}
	
	public void mul() {
		this.a = pop();
		this.b = pop();
		Binary<Integer> bin2 = factory1.newBinary(a, b, '*');
		data[top++] = bin2;
	}
	
	public void div() {
		this.a = pop();
		this.b = pop();
		Binary<Integer> bin2 = factory1.newBinary(b, a, '/');
		data[top++] = bin2;
	}
	
	public Expr<Integer> result() {
		return pop();
	}
}

//========================================================================================
//========================================================================================
//========================================================================================

public class Test 
{ 
        public static void main(String[] args) 
        { 
                try { 
                        java.util.Scanner input = 
                                new java.util.Scanner(System.in); 
                        int x = input.nextInt(); 
                        input.nextLine(); 
                        String expr = input.nextLine(); 
 
                        ExprBuilder builder = 
                                new ExprBuilder(new ExprFactory<Integer>()); 
                        ParsingDirector dir = new ParsingDirector(builder); 
                        dir.parse(expr); 
 
                        EvalVisitor v = new EvalVisitor(x); 
                        System.out.println(v.eval(builder.result())); 
                } catch (java.lang.Exception e) { 
                        System.out.println("exception:␣" + e); 
                } 
        } 
}

//========================================================================================
//========================================================================================
//========================================================================================

test1:
=> 50
=> 10-x
<= -40

test2:
=> 50
=> 100-80-x
<= -30

test3:
=> 20
=> -x+2*-(3+4*-(5+6*-(7+8*-(9+10*-(11+12)))))
<= -85186

