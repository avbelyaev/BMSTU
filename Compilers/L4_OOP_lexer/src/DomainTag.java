public enum DomainTag {
    IDENT(0),
    KEYWORD(1),
    EOP(2);

    private final int val;

    DomainTag(int val) {
        this.val = val;
    }

    public int getVal() {
        return val;
    }
}
