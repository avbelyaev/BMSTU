import java.util.ArrayList;

/**
 * Created by anthony on 02.10.16.
 */
public class BranchStmt extends Stmt{

    private String condition;
    private Vertex positive;
    private Vertex negative;

    @Override
    public void renameRhsVar(String varName, int varVersion) {
        //this.lhs = new Var(varName, varVersion);
        this.lhs.version = varVersion;
    }

    @Override
    public boolean updateRhsVarVersion(int version, int indexInRhs) {
        return false;
    }

    BranchStmt(Var left, String condition, Var right, Vertex positive, Vertex negative) {
        this.lhs = new Var(left.name);
        this.rhs = new Expr();
        this.rhs.vars.add(right);

        this.negative = negative;
        this.positive = positive;

        this.isPhi = false;
        this.isAss = false;
        this.condition = condition;
    }

    @Override
    public String toString() {
        String s = this.lhs.name + "_" + this.lhs.version + " " + this.condition + " " + this.rhs.vars.get(0).name +
            " ? " + positive.name.toUpperCase() + " : " + negative.name.toUpperCase() + " ";
        return s;
    }
}
