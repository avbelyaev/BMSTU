import java.io.File;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class complab3 {
    //var 11
    public static void main(String args[]) {

        String text = "";
        Scanner scanner;

        try {
            scanner = new Scanner(new File("src/in.txt"));
        } catch (java.io.FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }

        System.out.println("cols - [123456789]");
        int linenum = 1;
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine();
            text += line + "\n";
            System.out.println("line " + linenum + " [" + line + "]");
            linenum++;
        }

        //Комментарии: начинаются с «(∗» или «{», заканчиваются на «∗)» или «}»
        // и могут пересекать границы строк текста.
        String comment_1 = "(\\{[^\\}]*\\})";
        String comment_2 = "(\\(\\*([^\\*]|\\*($|[^\\)]))*\\*\\))";
        // exclude "not|this": ^([^nt]|n($|[^o]|o($|[^t]))|t($|[^h]|h($|[^i]|i($|[^s]))))*$

        //Идентификаторы: последовательности латинских букв,
        // представляющие собой конкатенации двух одинаковых слов («zz», «abab»).
        String ident = "(([A-Za-z]+)\\10)";

        //Ключевые слова: « ifif », «do», «dodo».
        String keyword = "(ifif|dodo|do)";

        String whitespace = "(\r| |\t)";

        String linebreak = "\n";



        String pattern = "(^" + comment_1 + "|^" + comment_2 +       // 1 2 3 4 5
                        ")|(^" + keyword +                           // 6 7
                        ")|(^" + ident +                             // 8 9 10
                        ")|(^" + whitespace +                        // 11 12
                        ")|(^" + linebreak + ")";                    // 13

        Pattern p = Pattern.compile(pattern);

        int line = 1, linepos = 1;
        boolean err = false;

        while (!text.equals("")) {

            Matcher m = p.matcher(text);
            if (m.find()) {

                err = false;

                if (null != m.group(1)) {
                    String tmp = m.group(1);

                    System.out.println("COMMENT (" + line + ", " + linepos + "): " + tmp.replaceAll("\n", " "));

                    while (tmp.contains("\n")) {
                        line++;
                        linepos = 1;
                        tmp = tmp.substring(tmp.indexOf("\n") + 1);
                    }
                    linepos += tmp.length();
                    text = text.substring(m.end());
                    continue;
                }

                if (null != m.group(6))
                    System.out.println("KEYWORD (" + line + ", " + linepos + "): " + m.group(6));

                if (null != m.group(8))
                    System.out.println("IDENT (" + line + ", " + linepos + "): " + m.group(8));



                if (null != m.group(13)) {
                    //System.out.println("lb (" + line + ", " + linepos + ")");
                    line++;
                    linepos = 0;
                }

                //if (null != m.group(11))
                //    System.out.println("ws (" + line + ", " + linepos + ")");


                text = text.substring(m.end());
                linepos += m.end();

            } else {

                if (true != err) {
                    System.out.println("SYNTAX ERR (" + line + ", " + linepos + ")");
                    err = true;
                }

                text = text.substring(1);
                linepos++;
            }


        }

        System.out.println("\ndone");
    }

}