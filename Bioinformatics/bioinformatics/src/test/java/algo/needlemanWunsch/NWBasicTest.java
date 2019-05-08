package algo.needlemanWunsch;

import algo.needlemanWunsch.INeedlemanWunsch;
import org.junit.Test;
import algo.needlemanWunsch.NWBasic;
import util.Common;
import util.Const;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static util.Const.basicSimilarity;

/**
 * Created by anthony on 28.10.16.
 */
public class NWBasicTest {
    private List<String> runBasic(String A, String B, int d) {
        INeedlemanWunsch nw = new NWBasic(A, B, d, basicSimilarity, Const.basicAmeno);
        System.out.println(nw);

        List<String> res = nw.computeAlignment();

        Common.printSequences(A, B, res);

        return res;
    }

    @Test
    public void test1() {
        List<String> r = runBasic("GCATGCAT", "GATTACA", -5);

        assertEquals("GCA-TGCAT", r.get(0));
        assertEquals("G-ATTACA-", r.get(1));
    }

    @Test
    public void test3() {
        List<String> r = runBasic("AGGA", "GA", -5);

        assertEquals("GCA-TGCAT", r.get(0));
        assertEquals("G-ATTACA-", r.get(1));
    }

    @Test
    public void test2() {
        List<String> r = runBasic("GCATGGCATACGT", "GATCATCGGACACAA", -5);

        assertEquals("G--CAT-GGCATACGT", r.get(0));
        assertEquals("GATCATCGG-ACACAA", r.get(1));
    }
}
