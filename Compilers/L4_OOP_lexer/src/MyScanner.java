import java.util.ArrayList;
import java.util.List;

public class MyScanner {
    public final String program;

    private Compiler compiler;
    private Position cur;
    private List<Fragment> comments;

    public void outputComments() {
        System.out.println();
        System.out.println("Comments:");

        for (Fragment f : comments) {
            int commStart = f.starting.getIndex();
            int commEnd = f.following.getIndex();

            System.out.println(f.toString() + ": " + program.substring(commStart, commEnd).replaceAll("\n", " "));
        }
    }

    public MyScanner(String program, Compiler compiler) {
        this.compiler = compiler;
        this.cur = new Position(program);
        this.program = program;
        this.comments = new ArrayList<>();
    }

    private boolean isValidIdent(String ident) {
        Position identPos = new Position(ident);
        boolean letter = false, digit = false;

        if (identPos.isLetter()) digit = true;  //assume "previous char" was a digit
        if (identPos.isDigit()) letter = true;  //-||- letter

        while (identPos.isLetterOrDigit()) {
            if (identPos.isLetter()) {  //letter found
                if (digit) {            //if "prev char" on current iteration is digit then
                    letter = true;      //"prev char" on next iteration will be a letter
                    digit = false;

                    identPos.nxt();
                    continue;
                }
                if (letter)
                    return false;
            }
            if (identPos.isDigit()) {
                if (letter) {
                    letter = false;
                    digit = true;

                    identPos.nxt();
                    continue;
                }
                if (digit)
                    return false;
            }
        }
        return true;
    }

    //Var 2
    //Комментарии: начинаются с «/∗», заканчиваются на «∗/» и могут пересекать границы строк текста.
    //Идентификаторы: последовательности латинских букв и десятичных цифр,
    // в которых буквы и цифры чередуются.
    //Ключевые слова: «for», «if», «m1».

    public Token nextToken() {
        while (-1 != cur.getCp()) {

            while (cur.isWhiteSpace() || cur.isNewLine())
                cur.nxt();

            Position start = cur.copy();   //avoid referencing object (otherwise changing 'cur' also changes 'start')

            switch (cur.getCp()) {
                case '/':

                    cur.nxt();
                    if ('*' != cur.getCp()) {

                        Position curCopy = cur.copy();
                        compiler.addMessage(true, curCopy, "unexpected character in comment: '" +
                                (char)cur.getCp() + "' found, '*' expected");
                    }

                    do {
                        do {
                            cur.nxt();
                        } while ('*' != cur.getCp() && -1 != cur.getCp());
                        cur.nxt();

                    } while ('/' != cur.getCp() && -1 != cur.getCp());

                    if (-1 == cur.getCp())
                        compiler.addMessage(true, cur, "end of program found, '*/' expected");

                    cur.nxt();

                    Position curCopy = cur.copy();
                    comments.add(new Fragment(start, curCopy));


                    break;

                default:
                    if (cur.isLetterOrDigit()) {

                        String currWord = "";
                        while (cur.isLetterOrDigit()) {
                            currWord += (char)(cur.getCp());
                            cur.nxt();
                            if (currWord.equals("for") || currWord.equals("if") || currWord.equals("m1"))
                                return new KeyToken(currWord, start, cur);
                        }

                        if (isValidIdent(currWord))
                            return new IdentToken(compiler.addName(currWord), currWord, start, cur);

                        compiler.addMessage(true, start, "unexpected character in ident");
                    }
                    break;
            }
            cur.nxt();
        }
        return null;
    }
}