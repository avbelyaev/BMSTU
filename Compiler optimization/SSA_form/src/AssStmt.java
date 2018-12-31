/**
 * Created by anthony on 02.10.16.
 */
public class AssStmt extends Stmt{


    AssStmt(Var lhs, Expr rhs) {
        this.lhs = lhs;
        this.rhs = rhs;
        this.isPhi = false;
        this.isAss = true;
    }

    @Override
    public String toString() {
        String s = lhs.toString();
        for (Var v : rhs.vars) s += v.toString();
        return s;
    }

    @Override
    public boolean updateRhsVarVersion(int version, int indexInRhs) {
        return false;
    }
}
