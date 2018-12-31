public class Fragment {
    public final Position starting, following;

    Fragment(Position starting, Position following) {
        this.starting = starting;
        this.following = following;
    }

    @Override
    public String toString() {
        return starting.toString() + "-" + following.toString();
    }
}
