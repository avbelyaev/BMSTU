public interface StackMachine 
{ 
        // Положить константу в стек 
        void pushConst(int value); 
 
        // Положить в стек значение переменной x 
        void pushVar(); 
 
        // Изменить знак числа на вершине стека 
        void neg(); 
 
        // Бинарные арифметические операции. 
        // Каждая операция снимает два операнда со стека 
        // и кладёт на стек результат 
        void add();  // сложение 
        void sub();  // вычитание 
        void mul();  // умножение 
        void div();  // деление 
 
        // Снимает со стека и возвращает число 
        int result(); 
}

//========================================================================================

public class Evaluator implements StackMachine{
        private int[] data;
	private int top, a, b, c, var;
	public Evaluator (int x) {
		data = new int[100];
		top = 0;
		var = x;
	}
	
	private int pop () {
		return data[--top];
	}
	
	public void pushConst (int value) {
		data[top++] = value;
	}
	
	public void pushVar () {
		data[top++] = var;
	}
	
	public void neg () {
		this.pushConst(pop() * -1);
	}
	
	public void add () {
		this.pushConst(pop() + pop());
	}
	
	public void sub () {
		b = pop();
		a = pop();
		this.pushConst(a-b);
	}	
	
	public void mul () {
		this.pushConst(pop() * pop());
	}
	
	public void div () {
		b = pop();
		a = pop();
		this.pushConst(a/b);
	}
	
	public int result () {
		return pop();
	}
}

//========================================================================================

public class Test 
{ 
        public static void main(String[] args) 
        { 
                try { 
                        Evaluator eval = new Evaluator(10); 
 
                        // x * (x+1) 
                        eval.pushVar(); 
                        eval.pushVar(); 
                        eval.pushConst(1); 
                        eval.add(); 
                        eval.mul(); 
 
                        System.out.println(eval.result()); 
                } catch (java.lang.Exception e) { 
                        System.out.println("exception:␣" + e); 
                } 
        } 
}

