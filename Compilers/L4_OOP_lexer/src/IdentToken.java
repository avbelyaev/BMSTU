public class IdentToken extends Token {
    public int code;

    public IdentToken(int code, String image, Position starting, Position following) {
        super(DomainTag.IDENT, image, starting, following);
        this.code = code;
    }

    @Override
    public String toString() {
        return super.toString();
    }
}
