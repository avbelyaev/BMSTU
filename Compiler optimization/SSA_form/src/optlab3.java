import java.util.*;

public class optlab3 {

    Map<Vertex, Set<Vertex>> dominanceFrontier;
    Map<Var, Set<Vertex>> varToUsageVertices;
    Set<Var> uniqueVaribales;
    Set<Vertex> postOrder;
    Vertex root;
    List<Integer> stack;
    int counter;

    private boolean compareVertexSets(Set<Vertex> a, Set<Vertex> b) {
        return a.size() == b.size();
    }

    private void postOrder(Vertex v) {
        v.succs.forEach(this::postOrder);
        postOrder.add(v);
    }

    private void generatePostOrder(Vertex root) {
        System.out.println("Post order:");

        postOrder(root);

        postOrder.forEach(v -> System.out.print(v.name.toUpperCase() + " "));
        System.out.println();
    }

    private void generateVarToUsageVertices() {

        this.postOrder.forEach(v -> {
            v.statements.forEach(stmt -> {
                Set<Vertex> additionalSet = varToUsageVertices.containsKey(stmt.lhs) ?
                        varToUsageVertices.get(stmt.lhs) : new HashSet<Vertex>();

                additionalSet.add(v);
                varToUsageVertices.put(stmt.lhs, additionalSet);
            });
        });


        System.out.println("Var to UsageVertices:");
        for (Var v : varToUsageVertices.keySet()) {
            Set<Vertex> vSet = varToUsageVertices.get(v);
            String s = "";
            for (Vertex x : vSet)
                s += x.name.toUpperCase() + " ";
            System.out.println("var: " + v.name + ": [" + s + "]");
        }
        System.out.println();
    }

    private void generateDFGlobal() {

        postOrder.forEach(x -> {
            dominanceFrontier.put(x, new HashSet<Vertex>());
            x.succs.forEach(succ -> {
                if (x != succ.immediateDom) {
                    dominanceFrontier.get(x).add(succ);

                    System.out.println("added succ " + succ.name);
                }
            });
            x.children.forEach(child -> {
                if (null != dominanceFrontier.get(child)) {
                    dominanceFrontier.get(child).forEach(y -> {
                        if (x != y.immediateDom) {
                            dominanceFrontier.get(x).add(y);

                            System.out.println("added child " + y.name);
                        }
                    });
                }
            });
        });


        System.out.println("Dominance Frontier:");
        for (Vertex v : postOrder) {
            Set<Vertex> vSet = dominanceFrontier.get(v);
            String s = "";
            for (Vertex x : vSet)
                s += x.name.toUpperCase() + " ";
            System.out.println("vert: " + v.name.toUpperCase() + ": [" + s + "]");
        }
        System.out.println();
    }

    private Set<Vertex> generateDFSet(Set<Vertex> vertexSet) {

        Set<Vertex> res = new HashSet<>();
        vertexSet.forEach(x -> res.addAll(dominanceFrontier.get(x)));


        String src = "";
        for (Vertex x : vertexSet)
            src += x.name.toUpperCase() + " ";

        System.out.println("DF-Set for [" + src + "]:");
        String s = "";
        for (Vertex x : res)
            s += x.name.toUpperCase() + " ";
        System.out.println("[" + s + "]");

        return res;
    }

    private Set<Vertex> generateDFPSet(Set<Vertex> vertexSet) {

        Set<Vertex> res = new HashSet<Vertex>();
        Set<Vertex> DFP = generateDFSet(vertexSet);
        boolean change;

        do {
            change = false;

            DFP.addAll(vertexSet);
            DFP = generateDFSet(DFP);

            if (!compareVertexSets(DFP, res)) {
                res = new HashSet<>(DFP);
                change = true;
            }
        } while (change);

        return res;
    }

    private void placePhi() {

        for (Var k : varToUsageVertices.keySet()) {

            System.out.println(k.name + ":");
            Set<Vertex> phiSet = generateDFPSet(varToUsageVertices.get(k));

            for (Vertex x : phiSet) {

                System.out.println("placing PHI");
                PhiStmt phi = new PhiStmt(x, k);

                x.phis.add(phi);
                x.prependStmt(phi);
            }

        }
    }

    private int whichPred(Vertex childVertexToBeSearchedWithin,
                          Vertex ancestorVertexToBeSearhedFor) {

        return childVertexToBeSearchedWithin.precs.indexOf(ancestorVertexToBeSearhedFor);
    }

    String currTrav, prevTrav;

    private void traverse(Vertex v, Var p) {
        currTrav = v.name;
        //if (currTrav.equals(prevTrav)) return;

        System.out.println("\nTRV " + v.name.toUpperCase() +
                " [" + p.name + ", s: " + stack + ", ctr: " + counter + "]:");
        //System.out.println("curr: " + currTrav + " prev: " + prevTrav);

        for (Stmt stmt : v.statements) {
            System.out.println("  stmt: " + stmt);

            if (!stmt.isPhi) {
                System.out.print("    Rhs: " + stmt + " \t-> ");

                stmt.renameRhsVar(p.name, stack.get(stack.size() - 1));

                System.out.println(stmt);
            }

            if (stmt.isAss && stmt.lhs.equals(p)) {
                System.out.print("    Lhs: " + stmt + " \t-> ");

                stmt.renameLhsVar(p.name, counter);
                stack.add(counter);
                counter++;

                System.out.println(stmt + " ctr++");
            }
        }

        System.out.println("s: " + stack + ", ctr: " + counter);

        v.succs.forEach(succ -> {

            int j = whichPred(succ, v);

            if (-1 != j) {
                System.out.println("  " + j + " = whichPred(" + succ.name + ", " + v.name + ")");
                System.out.println("  phis:");

                succ.phis.forEach(phiStmt -> {
                    if (p.equals(phiStmt.lhs)) {

                        System.out.print("    " + phiStmt + " \t-> ");

                        phiStmt.updateRhsVarVersion(stack.get(stack.size() - 1), j);

                        System.out.println(phiStmt);
                    }

                });
            }
        });

        prevTrav = v.name;

        v.children.forEach(child -> traverse(child, p));

        v.statements.forEach(stmt -> {
            if (stmt.lhs.equals(p))
                stack.remove(stack.size() - 1);
        });
    }

    private void renameSingleVar(Var p) {

        stack.clear();
        stack.add(0);
        counter = 0;

        traverse(root, p);
    }

    private void renameVars() {
        this.varToUsageVertices.keySet().forEach(this::renameSingleVar);
    }

    /*
    //   [A]
    //  /   \
    // [B]  [C]
     */
    private void buildCFGSample1() {



        Vertex a = new Vertex("a");
        Vertex b = new Vertex("b");
        Vertex c = new Vertex("c");
        Vertex d = new Vertex("d");
        Vertex e = new Vertex("e");

        this.postOrder.add(b);
        this.postOrder.add(c);
        this.postOrder.add(a);

        //=================== A ===================
        Var leftPartA1 = new Var("x", "=");
        Expr rightPartA1 = new Expr();
        rightPartA1.vars.add(new Var("5", ";"));
        Expr rightPartA2 = new Expr();
        rightPartA2.vars.add(new Var("x", "-"));
        rightPartA2.vars.add(new Var("3", ";"));
        Expr rightPartA3 = new Expr();
        rightPartA3.vars.add(new Var("1", ";"));

        a.statements.add(new AssStmt(leftPartA1, rightPartA1));
        a.statements.add(new AssStmt(new Var("x", "="), rightPartA2));
        a.statements.add(new AssStmt(new Var("y", "="), rightPartA3));
        a.children.add(b);
        a.children.add(c);
        a.succs.add(b);
        a.succs.add(c);
        a.immediateDom = null;

        //=================== B ===================
        Var leftPartB = new Var("t", "=");
        Expr rightPartB1 = new Expr();
        rightPartB1.vars.add(new Var("x", "*"));
        rightPartB1.vars.add(new Var("2", ";"));
        Expr rightPartB2 = new Expr();
        rightPartB2.vars.add(new Var("y", ";"));

        b.statements.add(new AssStmt(new Var("y", "="), rightPartB1));
        b.statements.add(new AssStmt(new Var("w", "="), rightPartB2));
        b.precs.add(a);
        b.immediateDom = a;

        //=================== C ===================
        Expr rightPartC1 = new Expr();
        rightPartC1.vars.add(new Var("x", "-"));
        rightPartC1.vars.add(new Var("3", ";"));

        c.statements.add(new AssStmt(new Var("y", "="), rightPartC1));
        c.precs.add(b);
        c.immediateDom = a;


        //renameVar(a, new Var("x", ""));
        //renameVar(a, new Var("y", ""));
        root = a;
    }

    /*
    //   [A]
    //  /   \
    // [B]  [C]
    //  \   /
    //   [D]
     */
    private void buildCFGSample2() {



        Vertex a = new Vertex("a");
        Vertex b = new Vertex("b");
        Vertex c = new Vertex("c");
        Vertex d = new Vertex("d");

        root = a;

        //=================== A ===================
        Var leftPartA1 = new Var("x", "=");
        Expr rightPartA1 = new Expr();
        rightPartA1.vars.add(new Var("5"));
        Expr rightPartA2 = new Expr();
        rightPartA2.vars.add(new Var("x", "-"));
        rightPartA2.vars.add(new Var("3"));
        Expr rightPartA3 = new Expr();
        rightPartA3.vars.add(new Var("1"));

        a.statements.add(new AssStmt(leftPartA1, rightPartA1));
        a.statements.add(new AssStmt(new Var("x", "="), rightPartA2));
        a.statements.add(new AssStmt(new Var("y", "="), rightPartA3));
        a.statements.add(new BranchStmt(new Var("x"), "<", new Var("3"), b, c));
        a.children.add(b);
        a.children.add(d);
        a.children.add(c);
        a.succs.add(b);
        a.succs.add(c);
        a.succs.add(d);
        a.immediateDom = null;


        //=================== B ===================
        Var leftPartB = new Var("t", "=");
        Expr rightPartB1 = new Expr();
        rightPartB1.vars.add(new Var("x", "*"));
        rightPartB1.vars.add(new Var("2"));
        Expr rightPartB2 = new Expr();
        rightPartB2.vars.add(new Var("y"));

        b.statements.add(new AssStmt(new Var("y", "="), rightPartB1));
        b.statements.add(new AssStmt(new Var("w", "="), rightPartB2));
        b.precs.add(a);
        b.succs.add(d);
        b.immediateDom = a;


        //=================== C ===================
        Expr rightPartC1 = new Expr();
        rightPartC1.vars.add(new Var("x", "-"));
        rightPartC1.vars.add(new Var("3", "+"));
        rightPartC1.vars.add(new Var("y", "*"));
        rightPartC1.vars.add(new Var("y"));
        Expr rightPartC2 = new Expr();
        rightPartC2.vars.add(new Var("z", "+"));
        rightPartC2.vars.add(new Var("9"));
        Expr rightPartC3 = new Expr();
        rightPartC3.vars.add(new Var("y"));

        c.statements.add(new AssStmt(new Var("y", "="), rightPartC1));
        c.statements.add(new AssStmt(new Var("x", "="), rightPartC2));
        c.statements.add(new AssStmt(new Var("x", "="), rightPartC3));
        c.precs.add(b);
        c.succs.add(d);
        c.immediateDom = a;


        //=================== D ===================
        Expr rightPartD1 = new Expr();
        rightPartD1.vars.add(new Var("x", "-"));
        rightPartD1.vars.add(new Var("y"));
        Expr rightPartD2 = new Expr();
        rightPartD2.vars.add(new Var("x", "+"));
        rightPartD2.vars.add(new Var("y"));
        Expr rightPartD1Phi1 = new Expr();
        rightPartD1Phi1.vars.add(new Var("y"));

        d.statements.add(new AssStmt(new Var("w", "="), rightPartD1));
        d.statements.add(new AssStmt(new Var("z", "="), rightPartD2));
        d.precs.add(b);
        d.precs.add(c);
        d.immediateDom = a;

        generatePostOrder(a);
    }

    /*
    //    /[ A ]\
    //   /   |   \
    //  /    |    \
    // [B]  [C]   [D]
    //  \    |    /
    //   \   |   /
    //    \[ E ]/
     */
    private void buildCFGSample3() {

        Vertex a = new Vertex("a");
        Vertex b = new Vertex("b");
        Vertex c = new Vertex("c");
        Vertex d = new Vertex("d");
        Vertex e = new Vertex("e");

        root = a;

        //=================== A ===================
        Expr rightPartA1 = new Expr();
        rightPartA1.vars.add(new Var("5"));
        Expr rightPartA2 = new Expr();
        rightPartA2.vars.add(new Var("x", "-"));
        rightPartA2.vars.add(new Var("3"));
        Expr rightPartA3 = new Expr();
        rightPartA3.vars.add(new Var("1"));

        a.statements.add(new AssStmt(new Var("x", "="), rightPartA1));
        a.statements.add(new AssStmt(new Var("w", "="), rightPartA2));
        a.statements.add(new AssStmt(new Var("y", "="), rightPartA3));
        a.statements.add(new BranchStmt(new Var("y"), "<", new Var("3"), b, c));
        a.children.add(b);
        a.children.add(c);
        a.children.add(d);
        a.children.add(e);
        a.succs.add(b);
        a.succs.add(c);
        a.succs.add(d);
        a.succs.add(e);


        //=================== B ===================
        Expr rightPartB1 = new Expr();
        rightPartB1.vars.add(new Var("x", "*"));
        rightPartB1.vars.add(new Var("(2", "-"));
        rightPartB1.vars.add(new Var("y", "+"));
        rightPartB1.vars.add(new Var("x", ")"));
        Expr rightPartB2 = new Expr();
        rightPartB2.vars.add(new Var("y"));

        b.statements.add(new AssStmt(new Var("y", "="), rightPartB1));
        b.statements.add(new AssStmt(new Var("w", "="), rightPartB2));
        b.precs.add(a);
        b.succs.add(e);
        b.immediateDom = a;


        //=================== C ===================
        Expr rightPartC1 = new Expr();
        rightPartC1.vars.add(new Var("x", "-"));
        rightPartC1.vars.add(new Var("3", "+"));
        rightPartC1.vars.add(new Var("y", "*"));
        rightPartC1.vars.add(new Var("y"));
        Expr rightPartC2 = new Expr();
        rightPartC2.vars.add(new Var("z", "+"));
        rightPartC2.vars.add(new Var("9"));
        Expr rightPartC3 = new Expr();
        rightPartC3.vars.add(new Var("y"));

        //c.statements.add(new AssStmt(new Var("y", "="), rightPartC1));
        c.statements.add(new AssStmt(new Var("w", "="), rightPartC2));
        //c.statements.add(new AssStmt(new Var("x", "="), rightPartC3));
        c.precs.add(a);
        c.succs.add(e);
        c.immediateDom = a;


        //=================== D ===================
        Expr rightPartD1 = new Expr();
        rightPartD1.vars.add(new Var("x", "-"));
        rightPartD1.vars.add(new Var("3"));
        Expr rightPartD2 = new Expr();
        rightPartD2.vars.add(new Var("x", "+"));
        rightPartD2.vars.add(new Var("y))"));
        Expr rightPartD3 = new Expr();
        rightPartD3.vars.add(new Var("z", "-"));
        rightPartD3.vars.add(new Var("3", "*"));
        rightPartD3.vars.add(new Var("(5", "+"));
        rightPartD3.vars.add(new Var("w", "+"));
        rightPartD3.vars.add(new Var("y", ")"));

        d.statements.add(new AssStmt(new Var("w", "="), rightPartD1));
        d.statements.add(new AssStmt(new Var("y", "="), rightPartD2));
        d.statements.add(new AssStmt(new Var("w", "="), rightPartD3));
        d.precs.add(a);
        d.succs.add(e);
        d.immediateDom = a;

        //=================== E ===================
        Expr rightPartE1 = new Expr();
        rightPartE1.vars.add(new Var("y"));
        Expr rightPartE2 = new Expr();
        rightPartE2.vars.add(new Var("x", "*"));
        rightPartE2.vars.add(new Var("w"));

        e.statements.add(new AssStmt(new Var("x", "="), rightPartE1));
        e.statements.add(new AssStmt(new Var("z", "="), rightPartE2));
        e.precs.add(b);
        e.precs.add(c);
        e.precs.add(d);
        e.immediateDom = a;

        generatePostOrder(a);
    }

    public static void main(String args[]) {
        optlab3 CFG = new optlab3();

        CFG.buildCFGSample3();

        CFG.generateDFGlobal();
        CFG.generateVarToUsageVertices();

        CFG.printBaseBlocks();

        CFG.placePhi();

        CFG.renameVars();
        //CFG.renameSingleVar(new Var("y"));

        CFG.printGraph(3);
        CFG.printBaseBlocks();
    }



    private optlab3() {
        this.postOrder = new LinkedHashSet<>();
        this.dominanceFrontier = new HashMap<Vertex, Set<Vertex>>();
        this.uniqueVaribales = new HashSet<Var>();
        this.varToUsageVertices = new HashMap<Var, Set<Vertex>>();
        this.stack = new ArrayList<Integer>();
    }

    private void printGraph(int n) {
        String str = "No graph representation";
        switch (n) {
            case 1:
                str = "\n" +
                "   [A]\n" +
                "  /   \\\n" +
                " [B]  [C]";
                break;
            case 2:
                str = "\n" +
                    "   [A]\n" +
                    "  /   \\\n" +
                    " [B]  [C]\n" +
                    "  \\   /\n" +
                    "   [D]";
                break;
            case 3:
                str = "\n" +
                    "    /[ A ]\\\n" +
                    "   /   |   \\\n" +
                    "  /    |    \\\n" +
                    " [B]  [C]   [D]\n" +
                    "  \\    |    /\n" +
                    "   \\   |   /\n" +
                    "    \\[ E ]/";
                break;
            default:
                break;
        }
        System.out.println(str);
    }

    private void printBaseBlocks() {
        System.out.println("\n->");
        dominanceFrontier.keySet().forEach(System.out::println);
    }

    private void buildCFGSample4() {

        //   [A]
        //  /   \
        // [B]  [C]
        //  \   /
        //   [D]
        //    |
        //   [E]

        Vertex a = new Vertex("a");
        Vertex b = new Vertex("b");
        Vertex c = new Vertex("c");
        Vertex d = new Vertex("d");
        Vertex e = new Vertex("e");

        root = a;

        //=================== A ===================
        Var leftPartA1 = new Var("x", "=");
        Expr rightPartA1 = new Expr();
        rightPartA1.vars.add(new Var("5"));
        Expr rightPartA2 = new Expr();
        rightPartA2.vars.add(new Var("x", "-"));
        rightPartA2.vars.add(new Var("3"));
        Expr rightPartA3 = new Expr();
        rightPartA3.vars.add(new Var("1"));

        a.statements.add(new AssStmt(leftPartA1, rightPartA1));
        a.statements.add(new AssStmt(new Var("x", "="), rightPartA2));
        a.statements.add(new AssStmt(new Var("y", "="), rightPartA3));
        a.statements.add(new BranchStmt(new Var("x"), "<", new Var("3"), b, c));
        a.children.add(b);
        a.children.add(c);
        a.children.add(d);
        a.children.add(e);
        a.succs.add(b);
        a.succs.add(c);
        a.succs.add(d);
        a.succs.add(e);
        a.immediateDom = null;


        //=================== B ===================
        Var leftPartB = new Var("t", "=");
        Expr rightPartB1 = new Expr();
        rightPartB1.vars.add(new Var("x", "*"));
        rightPartB1.vars.add(new Var("2"));
        Expr rightPartB2 = new Expr();
        rightPartB2.vars.add(new Var("y"));

        b.statements.add(new AssStmt(new Var("y", "="), rightPartB1));
        b.statements.add(new AssStmt(new Var("w", "="), rightPartB2));
        b.precs.add(a);
        b.succs.add(d);
        b.succs.add(e);
        b.immediateDom = a;


        //=================== C ===================
        Expr rightPartC1 = new Expr();
        rightPartC1.vars.add(new Var("x", "-"));
        rightPartC1.vars.add(new Var("3"));

        c.statements.add(new AssStmt(new Var("y", "="), rightPartC1));
        c.precs.add(b);
        c.succs.add(d);
        c.succs.add(e);
        c.immediateDom = a;


        //=================== D ===================
        Expr rightPartD1 = new Expr();
        rightPartD1.vars.add(new Var("x", "-"));
        rightPartD1.vars.add(new Var("y"));
        Expr rightPartD2 = new Expr();
        rightPartD2.vars.add(new Var("x", "+"));
        rightPartD2.vars.add(new Var("y"));

        d.statements.add(new AssStmt(new Var("w", "="), rightPartD1));
        d.statements.add(new AssStmt(new Var("z", "="), rightPartD2));
        d.precs.add(b);
        d.precs.add(c);
        d.succs.add(e);
        d.children.add(e);
        d.immediateDom = a;

        //=================== E ===================
        Expr rightPartE1 = new Expr();
        rightPartE1.vars.add(new Var("y"));
        Expr rightPartE2 = new Expr();
        rightPartE2.vars.add(new Var("2", "*"));
        rightPartE2.vars.add(new Var("x"));

        e.statements.add(new AssStmt(new Var("x", "="), rightPartE1));
        //e.statements.add(new AssStmt(new Var("z", "="), rightPartE2));
        e.precs.add(d);
        e.precs.add(b);
        e.immediateDom = d;

        generatePostOrder(a);
    }
}
