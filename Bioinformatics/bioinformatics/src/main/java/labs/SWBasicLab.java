package labs;

import algo.smithWaterman.ISmithWaterman;
import algo.smithWaterman.SWBasic;

import java.util.List;

import static util.Const.*;

/**
 * Created by anthony on 12.02.17.
 */
public class SWBasicLab {

    public static void main(String[] args) {

        String A ="GATTCA";
        String B = "GTCTGA";
        int penalty = -5;

        ISmithWaterman nw = new algo.smithWaterman.SWBasic(A, B, penalty, basicSimilarity, basicAmeno);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();

        System.out.println("in:");
        System.out.println(A + "\n" + B);
        System.out.println("out:");
        System.out.println(r.get(0) + "\n" + r.get(1));
    }
}
