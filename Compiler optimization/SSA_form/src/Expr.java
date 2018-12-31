import java.util.ArrayList;
import java.util.List;

/**
 * Created by anthony on 02.10.16.
 */
public class Expr {

    public List<Var> vars;

    Expr() {
        this.vars = new ArrayList<Var>();
        //this.vars.add(new Var("", ";"));
    }
    /*
    void append(Var vs) {
        Var enderVS = vars.get(vars.size() - 1);
        vars.set(vars.size() - 1, vs);
        vars.add(enderVS);
    }

    public Var getLast() {
        return vars.size() > 1 ? vars.get(vars.size() - 1) : null;
    }

    List<Var> getVars() {
        List<Var> clearList = new ArrayList<>(this.vars);
        clearList.remove(this.vars.size() - 1);
        return clearList;
    }

    Var getVarByIndex(int version) {
        return vars.size() > 1 ? vars.get(version) : null;
    }

    void setVarByIndex(int version, Var v) {
        this.vars.set(version, v);
    }
    */
    @Override
    public String toString() {
        String e = "";
        for (Var v : vars)
            e += v.toString();
        return e;
    }
}
