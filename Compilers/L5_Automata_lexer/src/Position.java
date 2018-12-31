/**
 * Created by anthony on 17.04.16.
 */
public class Position implements Comparable<Position> {

    private String text;
    private int line, pos, index;
    private int Cp;

    Position(String text) {
        this.text = text;
        line = pos = 1;
        index = 0;
    }

    public int getLine() {
        return line;
    }

    public void setLine(int line) {
        this.line = line;
    }

    public int getPos() {
        return pos;
    }

    public void setPos(int pos) {
        this.pos = pos;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public int getCp() {
        return (index == text.length()) ?
                -1 : Character.codePointAt(text.toCharArray(), index);
    }

    public Boolean isWhiteSpace() {
        return index != text.length() &&
                Character.isWhitespace(this.getCp());
    }

    public Boolean isLetter() {
        return index != text.length() &&
                Character.isLetter(this.getCp());
    }

    public Boolean isLetterOrDigit() {
        return index != text.length() &&
                Character.isLetterOrDigit(this.getCp());
    }

    public Boolean isNewLine() {
        if (index == text.length())
            return true;

        if ('\r' == text.charAt(index) && index + 1 < text.length())
            return ('\n' == text.charAt(index + 1));

        return ('\n' == text.charAt(index));
    }

    // position.nxt == position++
    public Position nxt() {
        if (index < text.length()) {
            if (isNewLine()) {
                if ('\r' == text.charAt(index)) {
                    index++;
                }
                line++;
                pos = 1;
            } else {
                if (Character.isHighSurrogate(text.charAt(index))) {
                    index++;
                }
                pos++;
            }
            index++;
        }
        return this;
    }

    @Override
    public int compareTo(Position o) {
        return Integer.compare(this.index, o.index);
    }

    @Override
    public String toString() {
        return "(" + line + ", " + pos + ")";
    }

    public Position copy() {
        Position tmp = new Position(this.text);
        tmp.setPos(this.pos);
        tmp.setLine(this.line);
        tmp.setIndex(this.index);

        return tmp;
    }
}