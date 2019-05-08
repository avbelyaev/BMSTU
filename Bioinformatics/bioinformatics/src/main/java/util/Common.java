package util;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Created by anthony on 07.01.17.
 */
public class Common {

    public static int max(int a, int b, int c, int d) {
        return max(a, max(b, c, d));
    }

    public static int max(int a, int b, int c) {
        return max(a, max(b, c));
    }

    public static int max(int a, int b) {
        return Math.max(a, b);
    }

    public static int similarityIndex(char c, char[] aminoBases) {
        List<Character> charList = new ArrayList<Character>();

        if ("-".equals(String.valueOf(c))) {
            return -2;
        }
        for (char e : aminoBases)
            charList.add(e);
        return -1 == charList.indexOf(String.valueOf(c).toUpperCase().charAt(0)) ?
                -1 : charList.indexOf(String.valueOf(c).toUpperCase().charAt(0));
    }

    public static int similarityIndex(char c) {
        switch (c) {
            case 'A':
                return 0;
            case 'G':
                return 1;
            case 'C':
                return 2;
            case 'T':
                return 3;
            case '-':
                return Const.starFlag;
            default:
                return -1;
        }
    }


    public static void printSequences(String a, String b, List<String> res) {
        System.out.println("\nin:");
        System.out.println(a);
        System.out.println(b);

        System.out.println("out:");
        System.out.println(res.get(0));
        System.out.println(res.get(1));
        System.out.println();
    }

    public static String printCell(int x) {
        String sx = String.valueOf(x);
        switch (sx.length()) {
            case 1:
                return "  " + sx;
            case 2:
                return " " + sx;
            case 3:
                return sx;
            default:
                return sx;
        }
    }
}
