package labs;

import algo.needlemanWunsch.INeedlemanWunsch;
import algo.needlemanWunsch.NWGapPenalty;

import java.util.List;

import static util.Const.basicAmeno;
import static util.Const.basicSimilarity;

/**
 * Created by anthony on 08.01.17.
 */
public class NWGapPenaltyLab {

    public static void main(String[] args) {

        String A ="GCATGTCGCAT";
        String B = "GATTACA";
        int penalty = -5;

        INeedlemanWunsch nw = new NWGapPenalty(A, B, penalty, basicSimilarity, basicAmeno);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();

        System.out.println("in:");
        System.out.println(A + "\n" + B);
        System.out.println("out:");
        System.out.println(r.get(0) + "\n" + r.get(1));
    }
}
