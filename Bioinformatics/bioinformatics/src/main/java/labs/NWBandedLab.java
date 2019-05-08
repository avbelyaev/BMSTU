package labs;

import algo.needlemanWunsch.INeedlemanWunsch;
import util.Common;
import util.Const;

import java.util.List;

/**
 * Created by anthony on 07.01.17.
 */
public class NWBandedLab {

    public static void main(String[] args) {

        String A = "GGGAAA";
        String B = "TGAC";
        int penalty = -5;
        int bandWidth = 3;

        INeedlemanWunsch nwb = new algo.needlemanWunsch.NWBanded(
                A, B, penalty, bandWidth, Const.basicSimilarity, Const.basicAmeno);
        System.out.println(nwb);

        List<String> r = nwb.computeAlignment();

        System.out.println("in:");
        System.out.println(A + "\n" + B);
        System.out.println("out:");
        System.out.println(r.get(0) + "\n" + r.get(1));
    }
}
