/**
 * Created by anthony on 02.10.16.
 */
public abstract class Stmt {
    Var lhs;
    Expr rhs;

    boolean isPhi;
    boolean isAss;



    public void renameRhsVar(String varName, int varVersion) {
        int i = 0;
        for (Var v : this.rhs.vars) {
            if (varName.equals(v.name)) {
                Var t = new Var(varName, varVersion);
                t.sign = v.sign;
                this.rhs.vars.set(i, t);
            }
            i++;
        }
    }

    public boolean renameLhsVar(String varName, int varVersion) {
        if (varName.equals(this.lhs.name)) {
            Var t = new Var(varName, varVersion);
            t.sign = "=";
            this.lhs = t;
            return true;
        }
        return false;
    }

    public abstract boolean updateRhsVarVersion(int version, int indexInRhs);

}
