import java.io.File;
import java.util.*;

/**
 * Created by anthony on 05.04.16.
 */
/*
digraph {
        rankdir = LR
        node [shape = doublecircle]; "5_KEY", "6_ID", "7_NUM", "8_OP", "9_WS";
        node [shape = circle]; "1", "2", "3", "4";
        0 [shape = circle]

        0 -> "1"                 [label = "d"]
        "1" -> "2"               [label = "e"]
        "2" -> "5_KEY"           [label = "f"]
        0 -> "3"                 [label = "v"]
        "3" -> "4"               [label = "a"]
        "4" -> "5_KEY"           [label = "r"]
        "4" -> "5_KEY"           [label = "l"]

        0 -> "6_ID"              [label = "a-z"]
        "6_ID" -> "6_ID"         [label = "a-z0-9"]

        0 -> "7_NUM"             [label = "0-9"]
        "7_NUM" -> "7_NUM"       [label = "0-9"]

        0 -> "8_OP"              [label = "["]
        0 -> "8_OP"              [label = "]"]

        0 -> "9_WS"              [label = "\\n\\r\\t' '"]
        "9_WS" -> "9_WS"         [label = "\\n\\r\\t' '"]
}

digraph {
        rankdir = LR
        node [shape = doublecircle]; "1_ID" "2_ID" "3_KEY" "4_ID" "5_ID"
                                     "6_ID" "7_NUM" "8_OP" "9_WS";
        0 [shape = circle]

        0 -> "1_ID"            [label = "d"]
        "1_ID" -> "2_ID"       [label = "e"]
        "2_ID" -> "3_KEY"      [label = "f"]
        0 -> "4_ID"            [label = "v"]
        "4_ID" -> "5_ID"       [label = "a"]
        "5_ID" -> "3_KEY"      [label = "r"]
        "5_ID" -> "3_KEY"      [label = "l"]
        0 -> "6_ID"            [label = "a-z\\{d,v}"]

        "1_ID" -> "6_ID"       [label = "a-z0-9\\{e}"]
        "2_ID" -> "6_ID"       [label = "a-z0-9\\{f}"]
        "4_ID" -> "6_ID"       [label = "a-z0-9\\{a}"]
        "5_ID" -> "6_ID"       [label = "a-z0-9\\{r,l}"]
        "3_KEY" -> "6_ID"      [label = "a-z0-9"]

        "6_ID" -> "6_ID"       [label = "a-z0-9"]

        0 -> "7_NUM"           [label = "0-9"]
        "7_NUM" -> "7_NUM"     [label = "0-9"]

        0 -> "8_OP"            [label = "["]
        0 -> "8_OP"            [label = "]"]

        0 -> "9_WS"            [label = "\\n\\r\\t' '"]
        "9_WS" -> "9_WS"        [label = "\\n\\r\\t' '"]
}
*/

public class complab5 {

    private class Automata {
        private SortedMap<Position, String> messages;
        public String program;
        private Position pos;
        private int state;

        public Automata(String program) {
            this.program = program;
            this.pos = new Position(program);
            this.state = 0;
            this.messages = new TreeMap<>();
        }

        private int get_code(char c) {
            if (c >= '0' && c <= '9')
                return 7;
            if (']' == c || '[' == c)
                return 8;
            if (' ' == c || '\n' == c || '\r' == c || '\t' == c)
                return 9;

            switch (c) {
                case 'd':
                    return 0;
                case 'e':
                    return 1;
                case 'f':
                    return 2;
                case 'v':
                    return 3;
                case 'a':
                    return 4;
                case 'r':
                    return 5;
                case 'l':
                    return 6;
            }

            if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))
                return 10;

            return -1;
        }

        private String get_state_name(int state) {
            switch (state) {
                case 1:
                case 2:
                    return "IDENT";
                case 3:
                    return "KEYWORD";
                case 4:
                case 5:
                case 6:
                    return "IDENT";
                case 7:
                    return "NUMBER";
                case 8:
                    return "OPERATION";
                case 9:
                    return "WHITESPACE";
                default:
                    return "ERROR";
            }
        }

        boolean err = false;

        public void run() {
            System.out.println("\nTokens:");
            while (-1 != pos.getCp()) {
                String word = "";
                state = 0;
                boolean final_state = false;
                Position start = pos.copy();

                while (-1 != pos.getCp()) {

                    char curr_char = program.charAt(pos.getIndex());
                    int jump_code = get_code(curr_char);

                    if (-1 == jump_code) {
                        if (!err) {
                            messages.put(pos.copy(), "Unexpected characters");
                            err = true;
                        }
                        break;
                    }
                    err = false;

                    System.out.print("(" + state + ")->");
                    System.out.print("[" + curr_char + "]->");

                    int next_state = table[state][jump_code];

                    if (-1 == next_state) {
                        final_state = true;
                        System.out.print("(-1)\n");
                        break;
                    }

                    word += curr_char;
                    state = next_state;
                    pos.nxt();
                }
                if (final_state) {

                    Fragment frag = new Fragment(start, pos);

                    System.out.println(get_state_name(state) + " " +
                            frag.toString() + ": " + word.replaceAll("\n"," "));

                    continue;
                }

                pos.nxt();
            }
        }

        public void output_messages() {
            System.out.println("\nMessages:");
            for (Map.Entry<Position, String> entry : messages.entrySet()) {
                System.out.print("ERROR ");
                System.out.print("(" + entry.getKey().getLine() + ", " +
                        entry.getKey().getPos() + "): ");
                System.out.println(entry.getValue());
            }
        }

    }



    final static int[][] table = {
                /*  d   e   f   v   a   r   l  num  par  ws  oth*/
    /*  START   */{ 1,  6,  6,  4,  6,  6,  6,  7,   8,   9,  6},
    /*  ID_1    */{ 6,  2,  6,  6,  6,  6,  6,  6,  -1,  -1,  6},
    /*  ID_2    */{ 6,  6,  3,  6,  6,  6,  6,  6,  -1,  -1,  6},
    /*  KEY_3   */{ 6,  6,  6,  6,  6,  6,  6,  6,  -1,  -1,  6},
    /*  ID_4    */{ 6,  6,  6,  6,  5,  6,  6,  6,  -1,  -1,  6},
    /*  ID_5    */{ 6,  6,  6,  6,  6,  3,  3,  6,  -1,  -1,  6},
    /*  ID_6    */{ 6,  6,  6,  6,  6,  6,  6,  6,  -1,  -1,  6},
    /*  NUM_7   */{-1, -1, -1, -1, -1, -1, -1,  7,  -1,  -1, -1},
    /*  OP_8    */{-1, -1, -1, -1, -1, -1, -1, -1,   8,  -1, -1},
    /*  WS_9    */{-1, -1, -1, -1, -1, -1, -1, -1,  -1,   9, -1}
    };

    //Var 10
    //Лексические домены:
    // Пробелы - \n, \r, \t, ' '
    // Идентификаторы - непустые последовательности латинских букв и десятичных цифр, начинающиеся с буквы
    // Целочисленные литералы - непустые последовательности десятичных цифр
    // Ключевые слова: def, val, var
    // Знаки операций: [, ]
    //P.S. чтобы не усложнять лексический анализатор,
    // разрешим идентификаторам примыкать справа к целочисленным литералам

    public static void main(String []args) {

        String text = "";
        Scanner scanner;

        try {
            scanner = new Scanner(new File("in.txt"));
        } catch (java.io.FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }

        System.out.println("   1234567890123456789");
        int i = 1;
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine();
            System.out.println(i + " [" + line + "]");
            text += line + "\n";
            i++;
        }

        Automata auto = new complab5().new Automata(text);

        auto.run();
        auto.output_messages();
    }
}