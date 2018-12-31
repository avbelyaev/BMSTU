/**
 * Created by anthony on 02.10.16.
 */
public class PhiStmt extends Stmt{

    //i-th argument corresponds to i-th node of ancestors set

    PhiStmt(Vertex v, Var p) {
        this.lhs = p;
        this.rhs = new Expr();
        this.isPhi = true;
        this.isAss = true;

        v.precs.forEach(prec -> this.rhs.vars.add(p));
        System.out.println(this.rhs.vars.size() + " vars were added into phi");
    }

    @Override
    public boolean updateRhsVarVersion(int version, int indexInRhs) {
        if (this.rhs.vars.size() <= indexInRhs || -1 == indexInRhs)
            return false;

        Var t = new Var(this.rhs.vars.get(indexInRhs).name);
        t.version = version;
        this.rhs.vars.set(indexInRhs, t);
        return true;
    }

    @Override
    public String toString() {
        String s =  lhs.toString() + "Ф( ";
        for (Var v : rhs.vars)
            s += v.name + "_" + v.version + " | ";
        s = s.substring(0, s.length() - 2) + ")";
        return s;
    }
}
