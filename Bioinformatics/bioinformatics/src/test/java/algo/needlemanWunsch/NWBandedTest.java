package algo.needlemanWunsch;

import algo.needlemanWunsch.INeedlemanWunsch;
import algo.needlemanWunsch.NWBanded;
import org.junit.Test;
import util.Common;
import util.Const;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static util.Const.basicSimilarity;

/**
 * Created by anthony on 12.02.17.
 */
public class NWBandedTest {

    @Test
    public void testNormalBand() {

        String A = "AGCTATC";
        String B = "CACTGAA";
        int d = -5;
        int k = 3;

        INeedlemanWunsch nw = new NWBanded(A, B, d, k, basicSimilarity, Const.basicAmeno);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("-AGCT-ATC", r.get(0));
        assertEquals("CA-CTGA-A", r.get(1));
    }

    @Test
    public void testWideBand() {

        String A = "ACAATCC";
        String B = "AGCATGC";
        int d = -5;
        int k = 6;

        INeedlemanWunsch nw = new NWBanded(A, B, d, k, basicSimilarity, Const.basicAmeno);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("A-CAATCC", r.get(0));
        assertEquals("AGC-ATGC", r.get(1));
    }
}
