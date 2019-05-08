package algo.needlemanWunsch;

import algo.needlemanWunsch.INeedlemanWunsch;
import algo.needlemanWunsch.NWGapPenalty;
import org.junit.Test;
import util.Common;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static util.Const.*;

/**
 * Created by anthony on 16.02.17.
 */
public class NWGapPenaltyTest {

    @Test
    public void testBasicSimilarity() {

        String A = "GCATGTCGCAT";
        String B = "GATTACA";
        int penalty = -5;

        INeedlemanWunsch nw = new NWGapPenalty(A, B, penalty, basicSimilarity, basicAmeno);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();
        Common.printSequences(A, B, r);

        assertEquals("GCATGTCGCAT", r.get(0));
        assertEquals("G-AT-T-ACA-", r.get(1));
    }

    @Test
    public void testBlosum45() {

        String A = "GCATGTCGCAT";
        String B = "GATTACA";
        int penalty = -5;

        INeedlemanWunsch nw = new NWGapPenalty(A, B, penalty, basicSimilarity, basicAmeno);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();
        Common.printSequences(A, B, r);

        assertEquals("GCATGTCGCAT", r.get(0));
        assertEquals("G-AT-T-ACA-", r.get(1));
    }


    @Test
    public void testDnaFull() {

        String A = "TATATAAATAATAATTAAAATA";
        String B = "ATATATTATATAT";
        int penalty = -5;

        INeedlemanWunsch nw = new NWGapPenalty(A, B, penalty, DNAFull, dnaFullAminoBases);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();
        Common.printSequences(A, B, r);

        assertEquals("TATATAAATAATAATTAAAATA", r.get(0));
        assertEquals("-ATAT--AT--T-A-TAT-AT-", r.get(1));
    }

    @Test
    public void testBlosum62Big() {
        String A = "AGTWAGTCGYKKRKAKSKWKABKHKKHWKWKKAK";
        String B = "AACHKABAVSWTAHKTVBA";
        int penalty = -5;

        INeedlemanWunsch nw = new NWGapPenalty(A, B, penalty, BLOSUM62, blosumAminoBases);
        System.out.println(nw);

        List<String> r = nw.computeAlignment();
        Common.printSequences(A, B, r);

        assertEquals("AGTWAGTCGYKKRKAKSKWKABKHKKHWKWKKAK", r.get(0));
        assertEquals("A---A--C--HKABAVS-WTA--H-K--T-VBA-", r.get(1));
    }
}
