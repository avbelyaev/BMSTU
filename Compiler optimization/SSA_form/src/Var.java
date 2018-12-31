/**
 * Created by anthony on 02.10.16.
 */

public class Var {

    String name;
    int version;
    String sign;

    String rightPar;

    Var(String name, int version) {
        this.rightPar = name.contains(")") ? ")+" : "";
        this.name = name.contains(")") ? name.substring(0, 1) : name;
        this.version = version;
        this.sign = "";
    }
    Var(String name, String sign) {
        this.rightPar = name.contains(")") ? ")+" : "";
        this.name = name.contains(")") ? name.substring(0, 1) : name;
        this.version = !"".equals(name) ? 0 : -1;
        this.sign = sign;
    }
    Var(String name) {
        this.rightPar = name.contains(")") ? ")+" : "";
        this.name = name.contains(")") ? name.substring(0, 1) : name;
        this.version = !"".equals(name) ? 0 : -1;
        this.sign = "";
    }

    private boolean isNumeric(String s) {
        return s.matches("[-+()]?\\d*\\.?\\d+");
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        return ((Var)o).name.equals(this.name);
    }

    @Override
    public int hashCode() {
        return (name.hashCode() < 0) ? -name.hashCode() : name.hashCode();
    }

    @Override
    public String toString() {
        String nameAndIndex = this.isNumeric(this.name) ? this.name : this.name + "_" + this.version + rightPar;

        return nameAndIndex + " " + this.sign + " ";
    }
}
