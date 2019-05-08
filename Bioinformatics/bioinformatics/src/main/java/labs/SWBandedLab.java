package labs;

import algo.smithWaterman.ISmithWaterman;
import algo.smithWaterman.SWBanded;

import java.util.List;

import static util.Const.basicAmeno;
import static util.Const.basicSimilarity;

/**
 * Created by anthony on 12.02.17.
 */
public class SWBandedLab {

    public static void main(String[] args) {

        String A ="TTGACACCCTCCCAATT";
        String B = "ACCCCAGGCTTTACACAT";
        int penalty = -5;
        int bandWidth = 2;

        ISmithWaterman nw = new SWBanded(A, B, penalty, bandWidth, basicSimilarity, basicAmeno);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();

        System.out.println("in:");
        System.out.println(A + "\n" + B);
        System.out.println("out:");
        System.out.println(r.get(0) + "\n" + r.get(1));
    }
}
