import java.util.*;

/**
 * Created by anthony on 02.10.16.
 */
class Vertex implements Comparable<Vertex>{
    String name;

    List<Vertex> precs;
    Set<Vertex> succs;          //potomki v CFG
    Set<Vertex> children;       //potomki v DOM tree

    Vertex immediateDom;        //the one and only and blizhaishiy

    List<Stmt> statements;
    List<Stmt> phis;

    public void prependStmt(Stmt s) {
        List<Stmt> tmp = new ArrayList<>(statements);
        this.statements = new ArrayList<>();
        this.statements.add(s);
        this.statements.addAll(tmp);
    }

    Vertex(String name) {
        this.precs = new ArrayList<Vertex>();
        this.succs = new HashSet<Vertex>();
        this.children = new HashSet<Vertex>();

        this.statements = new ArrayList<Stmt>();
        this.phis = new ArrayList<Stmt>();

        this.name = name;

        this.immediateDom = null;
    }

    @Override
    public int compareTo(Vertex o) {
        return this.name.compareTo(o.name);
    }

    @Override
    public int hashCode() {
        return name.hashCode() < 0 ? -name.hashCode() : name.hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        return ((Vertex)o).name.equals(this.name);
    }

    @Override
    public String toString() {
        String stmts = "\n ";
        for (Stmt s : statements)
            stmts += "\t" + s.toString() + "\n";
        stmts = stmts.substring(0, stmts.length() - 2);
        return name.toUpperCase() + ": {" + stmts + "}\n";
    }
}
