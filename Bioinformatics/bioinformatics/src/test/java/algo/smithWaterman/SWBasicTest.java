package algo.smithWaterman;

import org.junit.Test;
import util.Common;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static util.Const.*;

/**
 * Created by anthony on 12.02.17.
 */
public class SWBasicTest {

    @Test
    public void testBLOSUM45() {
        //http://docencia.ac.upc.edu/master/AMPP/slides/ampp_sw_presentation.pdf
        String A = "PAWHEAE";
        String B = "HEAGAWGHEE";
        int d = -5;

        ISmithWaterman sw = new SWBasic(A, B, d, BLOSUM45, blosumAminoBases);
        System.out.println(sw);

        List<String> r = sw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("----PAW-HEAE", r.get(0));
        assertEquals("HEAG-AWGHE-E", r.get(1));
    }

    @Test
    public void testBLOSUM62() {
        String A = "VSPAGMASGYD";
        String B = "IPGKASYD";
        int penalty = -5;

        ISmithWaterman nw = new SWBasic(A, B, penalty, BLOSUM62, blosumAminoBases);
        //System.out.println(nw);

        List<String> r = nw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("-VSPAGMASGYD", r.get(0));
        assertEquals("I--P-GKAS-YD", r.get(1));
    }

    @Test
    public void testBasicSimilarityBig() {
        String A = "TTGACACCCTCCCAATTGTA";
        String B = "ACCCCAGGCTTTACACAT";
        int penalty = -5;

        ISmithWaterman nw = new SWBasic(A, B, penalty, basicSimilarity, basicAmeno);
        //System.out.println(nw);

        List<String> r = nw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("TTGACACC---CTCC-CA-ATTGTA", r.get(0));
        assertEquals("---ACCCCAGGCTTTACACAT----", r.get(1));
    }
}
