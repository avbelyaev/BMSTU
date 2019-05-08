package labs;

import algo.needlemanWunsch.INeedlemanWunsch;

import java.util.List;

import static util.Const.basicAmeno;
import static util.Const.basicSimilarity;

/**
 * Created by anthony on 27.10.16.
 */
public class NWBasicLab {

    public static void main(String[] args) {

        String A ="GCATGCAT";
        String B = "GATTACA";
        int penalty = -5;

        INeedlemanWunsch nw = new algo.needlemanWunsch.NWBasic(A, B, penalty, basicSimilarity, basicAmeno);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();

        System.out.println("in:");
        System.out.println(A + "\n" + B);
        System.out.println("out:");
        System.out.println(r.get(0) + "\n" + r.get(1));
    }
}
