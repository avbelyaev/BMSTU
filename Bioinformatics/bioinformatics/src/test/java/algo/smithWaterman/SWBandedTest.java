package algo.smithWaterman;

import org.junit.Test;
import util.Common;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static util.Const.*;

/**
 * Created by anthony on 12.02.17.
 */
public class SWBandedTest {

    @Test
    public void testBasicNormalBand() {
        //http://docencia.ac.upc.edu/master/AMPP/slides/ampp_sw_presentation.pdf
        String A = "ACAATCAG";
        String B = "CTCATCCA";
        int d = -5;
        int k = 4;

        ISmithWaterman sw = new SWBanded(A, B, d, k, basicSimilarity, basicAmeno);
        System.out.println(sw);

        List<String> r = sw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("--ACAAT-CAG", r.get(0));
        assertEquals("CT-C-ATCCA-", r.get(1));
    }

    @Test
    public void testBLOSUM45WideBand() {
        //http://docencia.ac.upc.edu/master/AMPP/slides/ampp_sw_presentation.pdf
        String A = "PAWHEAE";
        String B = "HEAGAWGHEE";
        int d = -5;
        int k = 4;

        ISmithWaterman sw = new SWBanded(A, B, d, k, BLOSUM45, blosumAminoBases);
        //System.out.println(sw);

        List<String> r = sw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("----PAW-HEAE", r.get(0));
        assertEquals("HEAG-AWGHE-E", r.get(1));
    }

    @Test
    public void testBLOSUM62TightBand() {
        String A = "VSPAGMASGYD";
        String B = "IPGKASYD";
        int penalty = -5;
        int k = 2;

        ISmithWaterman nw = new SWBanded(A, B, penalty, k, BLOSUM62, blosumAminoBases);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("-VSPAGMAS--GYD", r.get(0));
        assertEquals("I--P-GKASYD---", r.get(1));
    }

    @Test
    public void testBasicBigTightBand() {
        String A = "TTGACACCCTCCCAATTGTA";
        String B = "ACCCCAGGCTTTACACAT";
        int penalty = -5;
        int k = 3;

        ISmithWaterman nw = new SWBanded(A, B, penalty, k, basicSimilarity, basicAmeno);
        //System.out.println(nw);

        List<String> r = nw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("TTGACACC---CTCC-CA-ATTGTA", r.get(0));
        assertEquals("---ACCCCAGGCTTTACACAT----", r.get(1));

    }

    @Test
    public void testBasicBigNormalBand() {
        String A = "TTGACACCCTCCCAATTGTA";
        String B = "ACCCCAGGCTTTACACAT";
        int penalty = -5;
        int k = 2;

        ISmithWaterman nw = new SWBanded(A, B, penalty, k, basicSimilarity, basicAmeno);
        //System.out.println(nw);

        List<String> r = nw.computeAlignment();

        Common.printSequences(A, B, r);

        assertEquals("----TTGACA-CCCTCCCA-ATTGTA", r.get(0));
        assertEquals("ACCC----CAGGCTTTACACAT----", r.get(1));
    }
}
