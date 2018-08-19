public class EvalVisitor implements ExprVisitor<Integer>{
        private int res;
	private int x1;
	public EvalVisitor(int x) {
		this.x1 = x;
	}
	
	public int eval (Expr<Integer> e) {
		e.accept(this);			//recognize subexpression
		return res;			//update result of it on level above
	}
        
	public void visitNeg(Neg<Integer> e) {
		res = -eval(e.a());
	}
        
	public void visitBinary(Binary<Integer> e) {
		if (e.operation() == '+') {
			res = eval(e.a()) + eval(e.b());	//dig into operand 
		} else if (e.operation() == '-') {      	//to check if its subexpr
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
		res = e.value();	// -||- => get and ret constant
	}
}

//========================================================================================

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
 
                EvalVisitor v = new EvalVisitor(5); 
                System.out.println(v.eval(e)); 
        } 
}

